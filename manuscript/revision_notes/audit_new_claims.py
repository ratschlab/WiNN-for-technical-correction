#!/usr/bin/env python3
"""Check new manuscript numbers directly against bundled exact CSV evidence."""
import csv
import gzip
import hashlib
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NOTES = ROOT / "revision_notes"
VALIDATION = ROOT / "generated/review_validation"
checks = []
sources = set()
texts = {}

def read_csv(path):
    path = Path(path); sources.add(path)
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))

def text(path):
    path = Path(path); sources.add(path)
    if path not in texts:
        texts[path] = path.read_text()
    return texts[path]

def one(rows, **keys):
    matches = [r for r in rows if all(r[k] == str(v) for k, v in keys.items())]
    assert len(matches) == 1, (keys, len(matches))
    return matches[0]

def check_token(claim, target, anchor, token, source, expected):
    blocks = [p for p in text(ROOT / target).split("\n\n") if anchor in p]
    assert len(blocks) == 1, (target, anchor, len(blocks))
    checks.append(dict(claim=claim, target=target, expected_display=token,
        source=str(source.relative_to(ROOT)), source_value=expected,
        passed=token in blocks[0]))

null = read_csv(VALIDATION / "null_gate_verification_summary.csv")
for n in (20, 50, 100, 200):
    for stage in ("unprepared", "mad_prepared"):
        r = one(null, n=n, stage=stage)
        if stage == "mad_prepared":
            check_token(f"Prepared gate global-null rejection n={n}", "supplementary.tex",
                "We tested whether WiNN selects features", f"{float(r['any_rejection'])*100:.1f}\\%",
                VALIDATION / "null_gate_verification_summary.csv", r["any_rejection"])
        else:
            for field in ("any_rejection", "any_rejection_lower", "any_rejection_upper", "rejected_features"):
                check_token(f"Unprepared gate n={n} {field}", "supplementary.tex",
                    "An independent null experiment", f"{float(r[field]):.3f}",
                    VALIDATION / "null_gate_verification_summary.csv", r[field])

# Standalone calibration comparison is supplementary only; benchmarks stay at 0.1.6.
calibration_path = ROOT / "generated/calibration_revision/gate_validation_summary.csv"
calibration = read_csv(calibration_path)
for n in (20, 50, 100, 200):
    r = one(calibration, experiment="original_null", n=n, stage="mad_prepared", algorithm="permutation_maxT")
    check_token(f"Alternative maxT null rate n={n}", "supplementary.tex",
        "We tested whether WiNN selects features", f"{float(r['fwer'])*100:.1f}\\%", calibration_path, r['fwer'])
    for field in ("fwer", "fwer_lower", "fwer_upper"):
        digits = 3 if field == "fwer" else 4
        check_token(f"Alternative null n={n} {field}", "supplementary.tex", "The paired null comparison reuses",
            f"{float(r[field]):.{digits}f}", calibration_path, r[field])
for rho in (0, 0.6):
    for algorithm in ("asymptotic_BH", "permutation_maxT"):
        r = one(calibration, experiment="robustness", n=20, rho=rho, distribution="gaussian",
            scenario="mixed_linear_drift", algorithm=algorithm)
        check_token(f"Alternative power rho={rho} {algorithm}", "supplementary.tex", "For Gaussian linear drift",
            f"{float(r['sensitivity']):.4f}", calibration_path, r['sensitivity'])

biology = read_csv(VALIDATION / "preparation_biology_summary.csv")
for scenario in ("balanced_sparse", "rare_subgroup", "broad_coordinated"):
    r = one(biology, scenario=scenario, stage="mad_then_pqn")
    value = float(r["median_recovery_ratio"])
    check_token(f"MAD+PQN recovery {scenario}", "supplementary.tex", "We tested whether WiNN selects features",
        f"{value:.3f}", VALIDATION / "preparation_biology_summary.csv", value)
    check_token(f"Supplement MAD+PQN recovery {scenario}", "supplementary.tex", "MAD followed by PQN retained",
        f"{value:.3f}", VALIDATION / "preparation_biology_summary.csv", value)
for scenario in ("rare_subgroup", "broad_coordinated"):
    for stage in ("mad_only", "pqn_only"):
        r = one(biology, scenario=scenario, stage=stage)
        check_token(f"Stage recovery {scenario} {stage}", "supplementary.tex", "MAD followed by PQN retained",
            f"{float(r['median_recovery_ratio']):.3f}", VALIDATION / "preparation_biology_summary.csv", r["median_recovery_ratio"])

retrieval = read_csv(VALIDATION / "frozen_identity_retrieval_estimates.csv")
for dataset, anchor in (("batchcorr_set1", "BatchCorr tested whether"), ("clinical_fiams", "The internal clinical FIA-MS cohort provided")):
    for method in ("raw", "winn_fixed_no_qc"):
        r = one(retrieval, dataset=dataset, metric="cross_batch_top1_accuracy", method_id=method)
        check_token(f"Query-weighted retrieval {dataset} {method}", "main.tex", anchor,
            f"{float(r['point_estimate'])*100:.1f}\\%", VALIDATION / "frozen_identity_retrieval_estimates.csv", r["point_estimate"])
        check_token(f"Supplement retrieval {dataset} {method}", "supplementary.tex", "Repeated samples in BatchCorr",
            f"{float(r['point_estimate']):.4f}", VALIDATION / "frozen_identity_retrieval_estimates.csv", r["point_estimate"])

differences = read_csv(VALIDATION / "frozen_identity_retrieval_paired_differences.csv")
for dataset in ("batchcorr_set1", "clinical_fiams"):
    r = one(differences, dataset=dataset, metric="cross_batch_top1_accuracy", comparator_method_id="raw")
    for field in ("point_direction_normalized_difference", "ci_lower", "ci_upper"):
        check_token(f"Frozen paired retrieval delta {dataset} {field}", "supplementary.tex", "Repeated samples in BatchCorr",
            f"{float(r[field]):.4f}", VALIDATION / "frozen_identity_retrieval_paired_differences.csv", r[field])

effects = read_csv(VALIDATION / "paired_feature_effects_descriptive.csv")
table_path = ROOT / "generated/supplementary_tables/table_s24_paired_feature_effects.tex"
table = text(table_path)
method_names = dict(combat="ComBat", qc_rlsc="QC-RLSC", qc_rfsc="QC-RFSC (legacy)", tiger="TIGER", serrf="SERRF",
    winn_auto_qc="WiNN (QC-tuned)", winn_auto_batch_qc="WiNN (inferred segments)", winn_fixed_no_qc="WiNN")
for r in effects:
    if r["method_id"] in ("raw", "qc_rfsc"):
        # The corrected RFSC fit is reported separately; legacy effects are not published.
        continue
    endpoint = {"mtbls79":"MTBLS79 cow--sheep", "waveica":"WaveICA groups"}.get(r["dataset_key"], "Sacurine " + r["endpoint"])
    rows = [line for line in table.splitlines() if line.startswith(endpoint + " & " + method_names[r["method_id"]] + " & ")]
    assert len(rows) == 1, (endpoint, r["method_id"])
    expected = " & ".join(f"{float(r[field]):.3f}" for field in ("slope_through_origin_vs_raw", "sign_concordance_vs_raw", "pearson_vs_raw"))
    checks.append(dict(claim=f"S24 matched-feature row {r['dataset_key']} {r['endpoint']} {r['method_id']}",
        target=str(table_path.relative_to(ROOT)), expected_display=expected,
        source="generated/review_validation/paired_feature_effects_descriptive.csv", source_value="three matched-feature summaries",
        passed=expected in rows[0]))
for dataset, endpoint in (("mtbls79","cow_vs_sheep"), ("waveica","biological_group"), ("sacurine","age"), ("sacurine","bmi"), ("sacurine","gender")):
    r = one(effects, dataset_key=dataset, endpoint=endpoint, method_id="winn_fixed_no_qc")
    for field in ("slope_through_origin_vs_raw", "sign_concordance_vs_raw"):
        check_token(f"Main fixed-WiNN feature effects {dataset} {endpoint} {field}", "main.tex",
            "We compared biological effects for the same features", f"{float(r[field]):.3f}", VALIDATION / "paired_feature_effects_descriptive.csv", r[field])

# Independently recompute S6 null tails from all 225,000 frozen counts.
distribution_path = ROOT / "generated/revision_panel_inputs/permutation_null.csv.gz"
sources.add(distribution_path)
distributions = defaultdict(list)
with gzip.open(distribution_path, "rt", newline="") as handle:
    for r in csv.DictReader(handle):
        distributions[(r["dataset_key"], r["endpoint"], r["method_id"])].append(int(r["significant_features"]))
def quantile(values, probability):
    index = (len(values)-1) * probability
    lo = math.floor(index); hi = math.ceil(index)
    return values[lo] + (index-lo) * (values[hi]-values[lo])
specs = [("A","mtbls79","cow_vs_sheep","mtbls79_cow_vs_sheep"),
    ("C","waveica","biological_group","waveica_biological_group"),
    ("E","sacurine","age","sacurine_age"), ("G","sacurine","bmi","sacurine_bmi"),
    ("I","sacurine","gender","sacurine_gender")]
for letter, dataset, endpoint, stem in specs:
    panel_path = ROOT / f"figures/source_data/S6{letter}_permutation_{stem}.csv"
    for r in read_csv(panel_path):
        values = sorted(distributions[(dataset,endpoint,r["method_id"])])
        observed = int(r["observed"])
        expect = dict(null_q99=quantile(values,.99), null_q999=quantile(values,.999), null_max=max(values),
            empirical_p=(1+sum(v >= observed for v in values))/(len(values)+1), permutations=len(values))
        for field, value in expect.items():
            checks.append(dict(claim=f"S6{letter} {r['method_id']} {field}", target=str(panel_path.relative_to(ROOT)),
                expected_display=str(value), source=str(distribution_path.relative_to(ROOT)), source_value=value,
                passed=abs(float(r[field])-value) < 1e-9))

# Verify the newly integrated QC-RFSC results against validated Euler summaries.
rfsc_summary_path = VALIDATION / "qcrfsc_completed_simulation_summary.csv"
rfsc_summary = read_csv(rfsc_summary_path)
# The main text no longer quotes QC-RFSC simulation values; the supplement checks follow.
for metric, digits in (("heldout_qc_cv_mean", 2), ("truth_log1p_rmse", 4), ("truth_metabolite_profile_icc_mean", 4), ("truth_sample_profile_pearson_mean", 4)):
    r = one(rfsc_summary, metric=metric)
    check_token("Supplement completed RFSC simulation " + metric, "supplementary.tex", "QC-RFSC yielded mean",
        f"{float(r['mean']):.{digits}f}", rfsc_summary_path, r["mean"])
rfsc_empirical_path = VALIDATION / "qcrfsc_completed_empirical_summary.csv"
for r in read_csv(rfsc_empirical_path):
    for field in ("raw_cv_percent", "qcrfsc_cv_percent"):
        check_token("Completed RFSC precision " + r["dataset_key"] + " " + field, "supplementary.tex", "QC-RFSC reduced held-out QC CV",
            f"{float(r[field]):.2f}\\%", rfsc_empirical_path, r[field])
    table_path = ROOT / "generated/supplementary_tables/table_s26_qcrfsc_empirical.tex"
    row = [line for line in text(table_path).splitlines() if line.startswith(r["dataset"] + " & ") and "multicolumn" not in line]
    assert len(row) == 1
    expected = f"{float(r['raw_cv_percent']):.3f} & {float(r['qcrfsc_cv_percent']):.3f} & {r['fallback_features']} & {r['retained_features']}"
    checks.append(dict(claim="S26 completed RFSC row " + r["dataset_key"], target=str(table_path.relative_to(ROOT)),
        expected_display=expected, source=str(rfsc_empirical_path.relative_to(ROOT)), source_value="precision and retained-feature summaries", passed=expected in row[0]))

with (NOTES / "new_claims_numerical_audit.csv").open("w", newline="") as handle:
    writer=csv.DictWriter(handle, fieldnames=list(checks[0])); writer.writeheader(); writer.writerows(checks)

timestamp = datetime.now(timezone.utc).isoformat()
sources.update(ROOT.glob("scripts/targeted_validation/**/*.R"))
sources.add(ROOT / "scripts/generate_revision_s6_panels.R")
sources.add(Path(__file__).resolve())
manifest=[]
for path in sorted(sources):
    raw=path.read_bytes()
    manifest.append(dict(path=str(path.relative_to(ROOT)), sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw), audit_utc=timestamp))
with (NOTES / "new_claims_source_hashes.csv").open("w", newline="") as handle:
    writer=csv.DictWriter(handle, fieldnames=list(manifest[0])); writer.writeheader(); writer.writerows(manifest)
failed=[r for r in checks if not r["passed"]]
for r in failed:
    print("FAILED:", r["claim"], "expected", r["expected_display"])
print(f"Audited {len(checks)} new numerical claims/source computations; {len(failed)} failed.")
raise SystemExit(bool(failed))
