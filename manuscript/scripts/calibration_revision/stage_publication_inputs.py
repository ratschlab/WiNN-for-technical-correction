#!/usr/bin/env python3
"""Stage complete recalculated WiNN panel inputs without changing comparator rows.

This performs deterministic aggregation only, never model fitting or selection.
All inputs are checked before the first output is written. --dry-run validates
and reports the pending files without modifying the manuscript checkout.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

WI = {"winn_fixed_no_qc", "winn_auto_qc", "winn_auto_batch_qc"}
LABEL_FIELDS = {"method", "method_label", "dataset", "dataset_label"}


def numeric(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Expected a finite value, got {value!r}")
    return number


def encoded(rows):
    return sorted(json.dumps(r, sort_keys=True, separators=(",", ":")) for r in rows)


def digest_rows(rows):
    return hashlib.sha256("\n".join(encoded(rows)).encode()).hexdigest()


def is_winn(row):
    if "method_id" in row:
        return row["method_id"] in WI
    return row.get("method", "").lower().startswith("winn")


def key(row, fields):
    return tuple(str(numeric(row[f])) if f == "level" else row[f] for f in fields)


def index(rows, fields):
    result = {}
    for row in rows:
        k = key(row, fields)
        if k in result:
            raise ValueError(f"Duplicate key {fields}: {k}")
        result[k] = row
    return result


def summary(values, expected=10):
    values = [numeric(v) for v in values]
    if len(values) != expected:
        raise ValueError(f"Expected {expected} paired observations, found {len(values)}")
    avg = statistics.mean(values)
    sd = statistics.stdev(values)
    se = sd / math.sqrt(len(values))
    return dict(mean=avg, sd=sd, se=se, n_seeds=len(values), lower=avg-se, upper=avg+se)


class Staging:
    def __init__(self, root, final, derivative):
        self.root = root
        self.final = final
        self.derivative = derivative
        self.panel = root / "generated/revision_panel_inputs"
        self.pending = {}
        self.sources = {}
        self.audit = []

    def read(self, path, source=True):
        path = Path(path)
        if not path.exists() and Path(str(path) + ".gz").exists():
            path = Path(str(path) + ".gz")
        data = path.read_bytes()
        plain = gzip.decompress(data) if path.suffix == ".gz" else data
        reader = csv.DictReader(io.StringIO(plain.decode()))
        rows = list(reader)
        if not reader.fieldnames or not rows:
            raise ValueError(f"Empty input: {path}")
        if source:
            self.sources[str(path)] = dict(sha256=hashlib.sha256(plain).hexdigest(), rows=len(rows))
        return list(reader.fieldnames), rows

    def table(self, stem, derivative=False):
        base = self.derivative if derivative else self.final
        path = base / "tables" / (stem + ".csv")
        if stem == "ablation_stage_impacts_table" and not path.exists() and not Path(str(path) + ".gz").exists():
            path = base / "tables/ablation_stage_impacts.csv"
        return self.read(path)[1]

    def stage(self, name, transform, all_rows=False, outside=False):
        path = self.root / name if outside else self.panel / name
        columns, original = self.read(path, source=False)
        updated = []
        changed = 0
        for old in original:
            if all_rows or is_winn(old):
                new = transform(old.copy())
                if set(new) != set(old):
                    raise ValueError(f"Output schema changed: {path}")
                updated.append(new)
                changed += 1
            else:
                updated.append(old)
        if not changed:
            raise ValueError(f"No target rows in {path}")
        before = [r for r in original if not is_winn(r)] if not all_rows else []
        after = [r for r in updated if not is_winn(r)] if not all_rows else []
        if encoded(before) != encoded(after):
            raise ValueError(f"Comparator changed: {path}")
        # A second equality check protects method-only S5 inputs, which the
        # historical method_id-only comparator audit does not cover.
        self.audit.append(dict(file=str(path.relative_to(self.root)), rows=len(updated),
                               refreshed_rows=changed, comparator_rows=len(before),
                               comparator_sha256=digest_rows(before)))
        self.pending[path] = (columns, updated)

    @staticmethod
    def merge(old, new):
        for field in old:
            if field in new and field not in LABEL_FIELDS:
                old[field] = new[field]
        return old

    def matching(self, filename, rows, fields, all_rows=False, outside=False):
        lookup = index(rows, fields)
        self.stage(filename, lambda old: self.merge(old, lookup[key(old, fields)]),
                   all_rows=all_rows, outside=outside)

    def bootstrap(self, dataset):
        return self.read(self.derivative / dataset / "bootstrap_replicates.csv")[1]


def prepare(s, expected_calibration="permutation"):
    params = s.table("primary_winn_selected_parameters")
    if len([r for r in params if r["method_id"] in WI]) != 18:
        raise ValueError("Expected 18 dataset/mode parameter records")
    for r in params:
        if r["method_id"] in WI:
            if (r.get("drift_calibration", "asymptotic") != expected_calibration or
                (expected_calibration == "permutation" and numeric(r.get("n_permutations", "nan")) != 999)):
                raise ValueError("Release does not match the declared calibration and permutation count")

    configuration = s.table("configuration_selection_manifest")
    winn_configuration = [r for r in configuration if r["method_id"] in WI]
    if len(winn_configuration) != 18:
        raise ValueError("Expected 18 WiNN configuration records")
    if expected_calibration == "permutation":
        for r in winn_configuration:
            expected_candidates = 1 if r["method_id"] == "winn_fixed_no_qc" else 36
            if numeric(r["candidate_count"]) != expected_candidates or r["source_version"] != "0.1.7.9999":
                raise ValueError("WiNN configuration metadata has stale version or candidate count")
    s.matching("generated/supplementary_table_data/endpoint_free_selection_manifest.csv",
               configuration, ["dataset_key", "method_id"], outside=True)

    all_safeguards = s.table("intensity_domain_safeguard")
    safeguard_key = ["analysis_family", "task_id", "dataset_key", "seed_id", "split_id", "method_id"]
    s.matching("generated/supplementary_table_data/intensity_domain_safeguard.csv",
               all_safeguards, safeguard_key, outside=True)
    safeguards = [r for r in all_safeguards if r["analysis_family"] == "primary"]
    s.matching("generated/review_validation/primary_intensity_domain_safeguard.csv",
               safeguards, safeguard_key,
               outside=True)

    # Long metric sources retain exact denominators needed for descriptive SEs.
    primary = index(s.table("primary_method_metrics_long"), ["dataset_key", "method_id", "metric"])
    boots = s.table("bootstrap_estimates", derivative=True)
    boot = index(boots, ["dataset_key", "method_id", "metric"])
    bootstrap_values = defaultdict(list)
    for ds in ("batchcorr_set1", "clinical_fiams", "mtbls79", "waveica", "sacurine"):
        for r in s.bootstrap(ds):
            if r["method_id"] in WI:
                bootstrap_values[(ds, r["method_id"], r["metric"])].append(numeric(r["value"]))

    def boot_stats(ds, method, metric):
        k = (ds, method, metric)
        point = numeric(boot[k]["point_estimate"])
        values = bootstrap_values[k]
        if len(values) != 2000 or numeric(boot[k]["valid_bootstrap"]) != 2000:
            raise ValueError(f"Expected 2,000 valid bootstrap draws: {k}")
        return point, statistics.stdev(values), len(values)

    simulation = s.table("simulation_seed_method_metrics")
    s.matching("simulation_seed_method_metrics.csv", simulation,
               ["dataset_key", "seed_id", "method_id", "metric"])
    ablation = s.table("simulation_seed_ablation_metrics")
    s.matching("simulation_seed_ablation_metrics.csv", ablation,
               ["dataset_key", "seed_id", "method_id", "metric"], all_rows=True)
    gate = [r for r in ablation if r["method_id"] in {"G_SS", "G_AS", "G_SA", "G_AA"}
            and r["metric"] == "truth_log1p_rmse"]
    gate_idx = index(gate, ["dataset_key", "seed_id", "method_id", "metric"])
    gate_stats = {m: summary(r["value"] for r in gate if r["method_id"] == m)
                  for m in {r["method_id"] for r in gate}}
    def gate_row(old):
        s.merge(old, gate_idx[key(old, ["dataset_key", "seed_id", "method_id", "metric"])])
        stat = gate_stats[old["method_id"]]
        old.update(method_mean=stat["mean"], method_se=stat["se"])
        return old
    s.stage("selective_forced_truth_rmse.csv", gate_row, all_rows=True)
    s.matching("selective_vs_forced_metrics_long.csv", s.table("selective_vs_forced_metrics_long"),
               ["dataset_key", "method_id", "metric"], all_rows=True)

    # Cumulative stage experiments are entirely WiNN-derived.
    stages = s.table("ablation_stage_impacts_table")
    stage_idx = index(stages, ["dataset_key", "metric", "stage"])
    max_benefit = defaultdict(float)
    for r in stages:
        if r["metric"] not in {"heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance", "batch_weighted_pc_r2_categorical"}:
            continue
        max_benefit[r["metric"]] = max(max_benefit[r["metric"]], abs(numeric(r["direction_normalized_benefit"])))
    def stage_row(old):
        s.merge(old, stage_idx[key(old, ["dataset_key", "metric", "stage"])])
        scale = max_benefit[old["metric"]]
        old["within_metric_score"] = 100 * numeric(old["direction_normalized_benefit"]) / scale if scale else 0
        return old
    s.stage("ablation_stage_impacts.csv", stage_row, all_rows=True)
    stage_fields = {"Raw": "cum_raw", "Outlier shrinkage": "cum_outlier", "Selective drift": "cum_drift",
                    "Selective batch": "cum_batch", "PQN": "cum_pqn"}
    def trajectory_row(old):
        r = stage_idx[(old["dataset_key"], "batch_weighted_pc_r2_categorical", "Outlier_shrinkage")]
        old["value"] = r[stage_fields[old["stage"]]]
        return old
    s.stage("categorical_batch_r2_trajectory.csv", trajectory_row, all_rows=True)

    confounding = [r for r in s.table("partial_confounding_method_metrics") if is_winn(r)]
    conf_groups = defaultdict(list)
    for r in confounding:
        conf_groups[(r["confounding_type"], numeric(r["level"]), r["metric"])].append(r["value"])
    def conf_row(old):
        level = numeric(old["level"])
        kind = "none" if level == 0 else old["confounding_type"]
        stat = summary(conf_groups[(kind, level, old["metric"])])
        for field in old:
            if field in stat:
                old[field] = stat[field]
        return old
    s.stage("confounding_fixed_winn.csv", conf_row, all_rows=True)
    for p in sorted(s.panel.glob("S5*_partial_*.csv")):
        s.stage(p.name, conf_row)

    reference = s.table("reference_split_method_metrics")
    reference_idx = index(reference, ["dataset_key", "method_id", "split_id", "metric"])
    reference_values = defaultdict(list)
    for r in reference:
        if r["metric"] == "heldout_qc_cv_mean":
            reference_values[(r["dataset_key"], r["method_id"])].append(r["value"])
    def ref_row(old):
        s.merge(old, reference_idx[key(old, ["dataset_key", "method_id", "split_id", "metric"])])
        old["method_mean"] = summary(reference_values[(old["dataset_key"], old["method_id"])])["mean"]
        return old
    for p in sorted(s.panel.glob("S4*_reference_*.csv")):
        s.stage(p.name, ref_row)
    for p in sorted(s.panel.glob("S8*.csv")):
        s.matching(p.name, boots, ["dataset_key", "method_id", "metric"])

    def scatter_row(old, identity=False):
        ds, method = old["dataset_key"], old["method_id"]
        metrics = ["batch_weighted_pc_r2_categorical",
                   "genuine_replicate_icc_median" if identity else "biological_weighted_pc_r2"]
        for metric in metrics:
            point, se, count = boot_stats(ds, method, metric)
            old["point_estimate_" + metric] = point
            old["bootstrap_se_" + metric] = se
            old["bootstrap_n_" + metric] = count
        old["batch_pc_variance_percent"] = 100 * old["point_estimate_batch_weighted_pc_r2_categorical"]
        old["batch_pc_variance_percent_se"] = 100 * old["bootstrap_se_batch_weighted_pc_r2_categorical"]
        if not identity:
            old["biological_pc_variance_percent"] = 100 * old["point_estimate_biological_weighted_pc_r2"]
            old["biological_pc_variance_percent_se"] = 100 * old["bootstrap_se_biological_weighted_pc_r2"]
        else:
            def metric(name):
                return primary[(ds, method, name)]
            cv = metric("heldout_qc_cv_mean")
            cv_sd = metric("heldout_qc_cv_sd")
            gam = metric("residual_run_order_gam_mean_deviance")
            gam_sd = metric("residual_run_order_gam_sd_deviance")
            if numeric(cv["denominator"]) != numeric(cv_sd["denominator"]):
                raise ValueError("QC CV mean/SD denominators disagree")
            if numeric(gam["denominator"]) != numeric(gam_sd["denominator"]):
                raise ValueError("GAM mean/SD denominators disagree")
            old.update(qc_cv=numeric(cv["value"]),
                       qc_cv_se=numeric(cv_sd["value"]) / math.sqrt(numeric(cv["denominator"])),
                       residual_run_order_gam_mean_deviance=numeric(gam["value"]),
                       residual_run_order_gam_sd_deviance=numeric(gam_sd["value"]),
                       retained_features=numeric(metric("retained_features")["value"]),
                       residual_gam_profiles=numeric(gam["denominator"]),
                       residual_run_order_gam_deviance_percent=100 * numeric(gam["value"]),
                       residual_run_order_gam_deviance_percent_se=100 * numeric(gam_sd["value"]) /
                           math.sqrt(numeric(gam["denominator"])))
        return old
    s.stage("eccb2026_identity_tradeoffs_v2.csv", lambda r: scatter_row(r, identity=True))
    s.stage("eccb2026_public_effect_tradeoffs_v2.csv", scatter_row)
    for letter, covariate in [("A", "age"), ("B", "bmi"), ("C", "gender")]:
        def sacurine_row(old, covariate=covariate):
            x, x_se, _ = boot_stats("sacurine", old["method_id"], "batch_weighted_pc_r2_categorical")
            y, y_se, _ = boot_stats("sacurine", old["method_id"], covariate + "_weighted_pc_r2")
            old.update(x=100*x, x_se=100*x_se, y=100*y, y_se=100*y_se)
            return old
        s.stage(f"S10{letter}_sacurine_{covariate}_weighted_pc_r2.csv", sacurine_row)

    associations = s.table("biological_association_summary", derivative=True)
    association_idx = index(associations, ["dataset_key", "method_id", "endpoint"])
    permutations = s.table("permutation_summary", derivative=True)
    s.matching("biological_association_summary.csv", associations, ["dataset_key", "method_id", "endpoint"])
    s.matching("permutation_summary.csv", permutations, ["dataset_key", "method_id", "endpoint"])
    nulls = []
    for ds in ("mtbls79", "waveica", "sacurine"):
        for row in s.read(s.derivative / ds / "permutation_null.csv")[1]:
            if row["method_id"] in WI:
                row["dataset_key"] = ds
                nulls.append(row)
    grouped_nulls = defaultdict(list)
    for r in nulls:
        grouped_nulls[(r["dataset_key"], r["method_id"], r["endpoint"])].append(numeric(r["significant_features"]))
    for r in permutations:
        if r["method_id"] not in WI:
            continue
        k = (r["dataset_key"], r["method_id"], r["endpoint"])
        observed = numeric(r["observed_count"])
        values = grouped_nulls[k]
        if len(values) != 5000 or observed != numeric(association_idx[k]["significant_features"]):
            raise ValueError(f"Incomplete or inconsistent permutation family: {k}")
        p = (1 + sum(v >= observed for v in values)) / 5001
        if abs(p - numeric(r["empirical_p_value"])) > 1e-12:
            raise ValueError(f"Permutation p-value mismatch: {k}")
    s.matching("permutation_null.csv.gz", nulls, ["dataset_key", "method_id", "endpoint", "permutation"])
    s.matching("generated/supplementary_table_data/leave_one_batch_out_bootstrap_summary.csv",
               s.table("leave_one_batch_out_bootstrap_summary", derivative=True),
               ["dataset_key", "method_id", "heldout_batch", "endpoint", "metric"], outside=True)
    s.matching("generated/supplementary_table_data/table_s12_lobo_fold_level.csv",
               s.table("leave_one_batch_out_bootstrap_summary", derivative=True),
               ["dataset_key", "method_id", "heldout_batch", "endpoint", "metric"], outside=True)



def prepare_review_validation(s, source_root):
    """Stage WiNN-only effects and all paired contrasts involving changed WiNN."""
    mappings = {
        "frozen_identity_retrieval_estimates.csv": (["dataset", "method_id", "metric"], False),
        "frozen_identity_retrieval_paired_differences.csv":
            (["dataset", "winn_method_id", "comparator_method_id", "metric"], True),
        "retrieval_query_weighted_bootstrap_replicates.csv.gz":
            (["dataset", "bootstrap", "method_id", "metric"], False),
        "baseline_simulation_summary.csv": (["method_id", "metric"], False),
        "baseline_simulation_seed_metrics.csv": (["seed_id", "method_id"], False),
        "baseline_simulation_paired_differences.csv": (["comparator", "metric"], True),
        "baseline_empirical_feature_effects.csv": (["dataset", "method_id", "endpoint"], False),
        "baseline_empirical_lobo_effects.csv": (["dataset", "method_id", "endpoint", "heldout_batch"], False),
        "paired_feature_effects_descriptive.csv": (["dataset_key", "method_id", "endpoint"], False),
        "clinical_batch_r2_bootstrap_bias_audit.csv": (["method_id", "metric"], False),
        "batchcorr_identity_retrieval_paired.csv": (["method_id", "endpoint"], False),
    }
    for name, (fields, all_rows) in mappings.items():
        rows = s.read(source_root / name)[1]
        s.matching("generated/review_validation/" + name, rows, fields,
                   all_rows=all_rows, outside=True)
    # This complete proof has repeated check labels; validate and replace the
    # ordered records directly instead of assuming labels are unique keys.
    name = "baseline_biology_validation_checks.csv"
    columns, rows = s.read(source_root / name)
    if len(rows) != 91 or not all(str(r["passed"]).lower() == "true" for r in rows):
        raise ValueError("Expected all 91 successful baseline biology validation checks")
    target = s.root / "generated/review_validation" / name
    s.pending[target] = (columns, rows)
    s.audit.append(dict(file=str(target.relative_to(s.root)), rows=len(rows),
                        refreshed_rows=len(rows), comparator_rows=0,
                        comparator_sha256=digest_rows([])))

    # These manifests identify corrected-matrix exports and must come from the
    # new analysis. Their rows have no method keys for safe numerical splicing.
    # biology_source_manifest belongs to unchanged preparation-only experiments.
    for name in ("frozen_derivative_input_manifest.csv",
                 "frozen_identity_retrieval_input_manifest.csv"):
        path = source_root / name
        columns, rows = s.read(path)
        target = s.root / "generated/review_validation" / name
        s.pending[target] = (columns, rows)
        s.audit.append(dict(file=str(target.relative_to(s.root)), rows=len(rows),
                            refreshed_rows=len(rows), comparator_rows=0,
                            comparator_sha256=digest_rows([])))


def serialize(columns, rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-root", type=Path, required=True)
    parser.add_argument("--derivative-root", type=Path,
                        help="Contains tables/ plus per-dataset bootstrap/permutation exports")
    parser.add_argument("--manuscript-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--review-validation-root", type=Path,
                        help="Fresh retrieval, matched-effect and baseline comparison exports")
    parser.add_argument("--package-commit", required=True)
    parser.add_argument("--analysis-commit", required=True)
    parser.add_argument("--expected-calibration", choices=["permutation", "asymptotic"], default="permutation",
                        help="Use asymptotic only to replay the old release for a no-op mapping audit")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    root = args.manuscript_root.resolve()
    final = args.final_root.resolve()
    derivative = (args.derivative_root or final / "derivative").resolve()
    staging = Staging(root, final, derivative)
    prepare(staging, args.expected_calibration)
    if len(staging.pending) != 38:
        raise ValueError("Expected all 38 main/supplementary numerical, configuration and safeguard input files")
    if args.review_validation_root:
        prepare_review_validation(staging, args.review_validation_root.resolve())
    session_info = None
    provenance_files = []
    if args.review_validation_root:
        source = args.review_validation_root.resolve() / "baseline_biology_sessionInfo.txt"
        session_info = source.read_bytes()
        if not session_info.strip():
            raise ValueError(f"Empty session provenance: {source}")
        checksum = hashlib.sha256(session_info).hexdigest()
        staging.sources[str(source)] = dict(sha256=checksum, bytes=len(session_info))
        provenance_files.append(dict(file="generated/review_validation/baseline_biology_sessionInfo.txt",
                                     sha256=checksum, bytes=len(session_info)))
    report = dict(package_commit=args.package_commit, analysis_commit=args.analysis_commit,
                  expected_calibration=args.expected_calibration,
                  review_validation_staged=bool(args.review_validation_root),
                  final_root=str(final), derivative_root=str(derivative),
                  sources=staging.sources, outputs=staging.audit,
                  provenance_files=provenance_files)
    if args.dry_run:
        print(json.dumps(report, indent=2))
        return
    for path, (columns, rows) in staging.pending.items():
        plain = serialize(columns, rows)
        path.write_bytes(gzip.compress(plain, mtime=0) if path.suffix == ".gz" else plain)
    if session_info is not None:
        (root / provenance_files[0]["file"]).write_bytes(session_info)
    # Update the existing uncertainty source map only after numerical staging.
    manifest = staging.panel / "supplementary_uncertainty_source_manifest.csv"
    columns, rows = staging.read(manifest, source=False)
    for row in rows:
        path = root / row["file"]
        row["original_source"] = f"permutation-calibrated WiNN release {args.analysis_commit}; comparator rows retained"
        row["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest.write_bytes(serialize(columns, rows))
    out = root / "generated/calibration_revision/publication_input_staging.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Staged {len(staging.pending)} files; all comparator dictionaries unchanged.")
    if not args.review_validation_root:
        print("Retrieval, matched-effect and baseline comparison exports were not staged; "
              "provide --review-validation-root before regenerating their tables.")


if __name__ == "__main__":
    main()
