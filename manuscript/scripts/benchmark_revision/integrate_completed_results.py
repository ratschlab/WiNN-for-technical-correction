#!/usr/bin/env python3
"""Validate retrieved Euler aggregates and update publication tables/plot inputs.

Inputs under generated/review_validation/euler_completed contain aggregate
outputs only. Private matrices and identity maps remain on Euler.
"""
import csv
import hashlib
import math
import re
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VALIDATION = ROOT / "generated/review_validation"
EULER = VALIDATION / "euler_completed"


def read(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows, columns=None):
    columns = columns or list(rows[0])
    with open(path, "w", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=columns)
        out.writeheader()
        out.writerows(rows)


checks = []


def check(name, passed, detail=""):
    checks.append(dict(check=name, passed=bool(passed), detail=detail))
    if not passed:
        raise ValueError(name + ": " + detail)


metrics = ["heldout_qc_cv_mean", "truth_log1p_rmse",
           "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean"]
source = ROOT / "generated/revision_panel_inputs/simulation_seed_method_metrics.csv"
old = read(source)
archive = VALIDATION / "qcrfsc_previous_simulation_metrics.csv"
if not archive.exists():
    write(archive, [r for r in old if r["method_id"] == "qc_rfsc"])
legacy = read(archive)
corrected = []
for n in range(1, 11):
    seed = f"SIM{n:02d}"
    directory = EULER / "simulation_results" / seed
    check(seed + " completed", "status=complete" in (directory / "COMPLETE.txt").read_text())
    identity = read(directory / "input_identity_checks.csv")
    check(seed + " exact raw and clean hashes", len(identity) == 2 and
          all(r["passed"] == "TRUE" and r["observed_sha256"] == r["expected_sha256"] for r in identity))
    assertions = read(directory / "REGfit_training_index_assertions.csv")
    check(seed + " training-reference isolation", bool(assertions) and all(
        int(r["unauthorized"]) == 0 and r["used"] == r["authorized"] for r in assertions))
    rows = read(directory / "simulation_truth_metrics.csv")
    for method in ["raw", "legacy"]:
        row = next(r for r in rows if r["method"] == method)
        original = old if method == "raw" else legacy
        for metric in metrics:
            expected = next(float(r["value"]) for r in original if r["seed_id"] == seed and
                            r["metric"] == metric and r["method_id"] == ("raw" if method == "raw" else "qc_rfsc"))
            error = abs(float(row[metric]) - expected)
            check(seed + " " + method + " " + metric + " reproduced", error < 1e-11, str(error))
    row = next(r for r in rows if r["method"] == "corrected")
    full = EULER / "simulation_evaluation" / seed / "metrics.csv"
    if full.exists():
        complete = (full.parent / "COMPLETE.txt").read_text()
        check(seed + " full endpoints complete", "status=complete" in complete)
        replacement = read(full)
        for metric in metrics:
            observed = next(float(r["value"]) for r in replacement if r["metric"] == metric)
            check(seed + " independent endpoint agreement " + metric,
                  abs(observed - float(row[metric])) < 1e-11)
    else:
        replacement = []
        for metric in metrics:
            original = next(r for r in legacy if r["seed_id"] == seed and r["metric"] == metric).copy()
            original["value"] = row[metric]
            replacement.append(original)
    for item in replacement:
        item = {key: item.get(key, "") for key in old[0]}
        item.update(method_id="qc_rfsc", method="QC-RFSC", dataset_key="simulation", seed_id=seed)
        corrected.append(item)
    # Validate completeness, never require a particular ranking of methods.
    # Calibration revisions may change WiNN estimates in either direction.
    for metric in metrics[1:]:
        matched = [r for r in old if r["seed_id"] == seed and
                   r["method_id"] == "winn_fixed_no_qc" and r["metric"] == metric]
        check(seed + " unique finite WiNN truth endpoint " + metric,
              len(matched) == 1 and math.isfinite(float(matched[0]["value"])))
        check(seed + " finite corrected RFSC truth endpoint " + metric,
              math.isfinite(float(row[metric])))
write(source, [r for r in old if r["method_id"] != "qc_rfsc"] + corrected)

# Every published RFSC simulation cell comes from the same corrected matrices.
table_path = ROOT / "generated/supplementary_table_data/table_s4_simulation.csv"
table = read(table_path)
summary = []
names = {"heldout_qc_cv_mean": "QC CV", "residual_run_order_gam_mean_deviance": "GAM",
         "batch_weighted_pc_r2_categorical": "Batch R2", "truth_log1p_rmse": "RMSE",
         "truth_metabolite_profile_icc_mean": "ICC", "truth_sample_profile_pearson_mean": "Pearson"}
target = next(r for r in table if r["method_id"] == "qc_rfsc")
for metric, label in names.items():
    values = [float(r["value"]) for r in corrected if r["metric"] == metric]
    if len(values) == 10:
        mean = statistics.mean(values)
        se = statistics.stdev(values) / math.sqrt(len(values))
        digits = 2 if label == "QC CV" else 4
        target[label] = f"{mean:.{digits}f} $\\pm$ {se:.{digits}f}"
        summary.append(dict(metric=metric, n=10, mean=mean, se=se))
    else:
        target[label] = "--"
write(table_path, table)
write(VALIDATION / "qcrfsc_completed_simulation_summary.csv", summary)
tex_path = ROOT / "generated/supplementary_tables/table_s4_simulation.tex"
lines = tex_path.read_text().splitlines()
replacement = " & ".join(["QC-RFSC"] + [target[k] for k in ["QC CV", "GAM", "Batch R2", "RMSE", "ICC", "Pearson"]]) + r" \\"
tex_path.write_text("\n".join(replacement if line.startswith("QC-RFSC &") else line for line in lines) + "\n")

# Remove invalid previous-adapter estimates from publication tables. The
# completed corrected empirical endpoints are reported separately (Table S26).
for number in [5, 6, 7, 8, 9, 10, 11, 12, 14, 23, 24]:
    for path in (ROOT / "generated/supplementary_tables").glob(f"table_s{number}_*.tex"):
        text = path.read_text()
        text = "\n".join(line for line in text.splitlines() if not ("QC-RFSC" in line and " & " in line)) + "\n"
        text = text.replace("across all nine outputs", "across eight outputs")
        text = text.replace("QC-RFSC legacy estimates remain provisional.", "")
        text = text.replace("QC-RFSC legacy comparisons remain provisional because its output-scale adapter is being corrected.", "")
        path.write_text(text)
    for path in (ROOT / "generated/supplementary_table_data").glob(f"table_s{number}_*.csv"):
        rows = read(path)
        if rows:
            kept = [r for r in rows if not any(v in ["qc_rfsc", "QC-RFSC", "QC-RFSC (legacy)"] for v in r.values())]
            write(path, kept, list(rows[0]))

labels = {"mtbls79": "MTBLS79", "batchcorr_set1": "BatchCorr Set 1",
          "clinical_fiams": "Clinical FIA-MS", "sacurine": "Sacurine", "waveica": "WaveICA"}
feature_counts = dict(mtbls79=2488, batchcorr_set1=199, clinical_fiams=1594, sacurine=113, waveica=6461)
empirical = []
all_empirical_endpoints = []
for dataset, label in labels.items():
    directory = EULER / "results" / dataset
    check(dataset + " fit complete", "status=complete" in (directory / "COMPLETE.txt").read_text())
    assertions = read(directory / "REGfit_training_index_assertions.csv")
    check(dataset + " training-reference isolation", bool(assertions) and all(
        int(r["unauthorized"]) == 0 and r["used"] == r["authorized"] for r in assertions))
    rows = read(directory / "revision_diagnostic_metrics.csv")
    raw = next(r for r in rows if r["method"] == "raw")
    result = next(r for r in rows if r["method"] == "corrected")
    selection = read(directory / "selected_parameters.csv")[0]
    check(dataset + " finite output", int(result["nonfinite"]) == 0)
    selected_scores = [r for r in read(directory / "selection_summary.csv") if r["selected"] == "TRUE"]
    check(dataset + " unique recorded candidate", len(selected_scores) == 1 and
          selected_scores[0]["candidate_id"] == selection["candidate_id"])
    endpoint_dir = EULER / "empirical_evaluation" / dataset
    check(dataset + " full point endpoints complete", "status=complete" in (endpoint_dir / "COMPLETE.txt").read_text())
    endpoints = read(endpoint_dir / "corrected/metrics.csv")
    endpoint_values = {r["metric"]: float(r["value"]) for r in endpoints if r["value"] not in ("NA", "", "NaN")}
    check(dataset + " evaluated reference precision agrees", abs(endpoint_values["heldout_qc_cv_mean"] - float(result["heldout_cv_percent"])) < 1e-11)
    check(dataset + " full endpoint dimensions", endpoint_values["retained_features"] == feature_counts[dataset])
    all_empirical_endpoints.extend(endpoints)
    empirical.append(dict(dataset_key=dataset, dataset=label,
                          raw_cv_percent=float(raw["heldout_cv_percent"]),
                          qcrfsc_cv_percent=float(result["heldout_cv_percent"]),
                          fallback_features=int(result["uncorrected_features"]),
                          retained_features=feature_counts[dataset],
                          selected_parameters=selection["params"],
                          residual_gam_percent=100 * endpoint_values["residual_run_order_gam_mean_deviance"],
                          batch_weighted_pc_r2=endpoint_values["batch_weighted_pc_r2_categorical"],
                          replicate_profile_icc=endpoint_values.get("genuine_replicate_icc_median", ""),
                          replicate_profile_pearson=endpoint_values.get("sample_replicate_pearson_median", "")))
write(VALIDATION / "qcrfsc_completed_empirical_summary.csv", empirical)
write(VALIDATION / "qcrfsc_completed_empirical_all_endpoints.csv", all_empirical_endpoints)
write(ROOT / "generated/supplementary_table_data/table_s26_qcrfsc_empirical.csv", empirical)
body = [r"\footnotesize", r"\begin{longtable}{lrrrrrr}",
        r"\caption{QC-RFSC empirical point estimates and feature retention. Mean CV is calculated across features in the held-out reference set. GAM is residual acquisition-order deviance explained; batch $R^2$ is categorical weighted-PC $R^2$ among study injections. Features without a returned correction retain their input abundance. Replicate-profile agreement is summarized across repeated biological identities. Clinical results use the reconstructed pair-level input described in the Methods.}\label{tab:supp_qcrfsc_empirical}\\",
        r"\toprule", r"Dataset & Raw CV (\%) & QC-RFSC CV (\%) & Raw-restored & Retained & GAM (\%) & Batch $R^2$ \\",
        r"\midrule", r"\endfirsthead", r"\toprule",
        r"Dataset & Raw CV (\%) & QC-RFSC CV (\%) & Raw-restored & Retained & GAM (\%) & Batch $R^2$ \\",
        r"\midrule", r"\endhead", r"\bottomrule", r"\endlastfoot"]
for row in empirical:
    body.append(f"{row['dataset']} & {row['raw_cv_percent']:.3f} & {row['qcrfsc_cv_percent']:.3f} & {row['fallback_features']} & {row['retained_features']} & {row['residual_gam_percent']:.3f} & {row['batch_weighted_pc_r2']:.5f}" + r" \\")
body += [r"\midrule", r"\multicolumn{7}{l}{Replicate-profile agreement in the repeated-identity datasets} \\"]
for row in empirical:
    if row['dataset_key'] in ('batchcorr_set1', 'clinical_fiams'):
        body.append(row['dataset'] + r" & \multicolumn{3}{l}{ICC(A,1): " + f"{row['replicate_profile_icc']:.4f}" + r"} & \multicolumn{3}{l}{Pearson correlation: " + f"{row['replicate_profile_pearson']:.4f}" + r"} \\")
body.append(r"\end{longtable}")
(ROOT / "generated/supplementary_tables/table_s26_qcrfsc_empirical.tex").write_text("\n".join(body) + "\n")

# Selected configurations now refer to the validated corrected RFSC runs.
path = ROOT / "generated/supplementary_tables/table_s21_selected_parameters.tex"
lines = path.read_text().splitlines()
for i, line in enumerate(lines):
    if line.startswith(r"\caption{"):
        lines[i] = r"\caption{Comparator settings by dataset. Training-QC selection uses the prespecified candidate grids; single candidates are fixed configurations.}\label{tab:supp_selected_parameters}\\"
    for row in empirical:
        if line.startswith(row["dataset"] + " & QC-RFSC &"):
            lines[i] = row["dataset"] + " & QC-RFSC & " + row["selected_parameters"] + r" & 4 & Training-QC score \\"
path.write_text("\n".join(lines) + "\n")

# Keep the machine-readable configuration manifest aligned with Table S21.
# Legacy configuration and seeds remain in actual_comparator_parameters.csv.
manifest_path = ROOT / "generated/supplementary_table_data/endpoint_free_selection_manifest.csv"
manifest = read(manifest_path)
for row in manifest:
    if row["method_id"] != "qc_rfsc" or row["dataset_key"] == "simulation":
        continue
    dataset = row["dataset_key"]
    selected = read(EULER / "results" / dataset / "selected_parameters.csv")[0]
    row["selected_parameters"] = selected["params"]
    key = "::".join((dataset, "QC-RFSC", selected["candidate_id"], "full_training", "post_selection"))
    row["selection_seed"] = str(int(hashlib.sha256(key.encode()).hexdigest()[:7], 16) + 1)
    row["selection_criterion"] = "recomputed training-QC score with neutral identifiers and abundance restoration"
    row["information_available_during_selection"] = "intensity matrix, acquisition order, supplied batches and authorized training-reference identities; withheld-reference identities and biological labels masked"
    row["departure_justification"] = "Selected from the four-candidate grid using the corrected adapter and training-QC score only."
write(manifest_path, manifest)

captions = {
    "table_s22_primary_safeguards.tex": r"Numerical safeguards in the benchmark fits. The table records negative intensities floored at zero and features restored to input abundance. QC-RFSC fallback counts are reported in Table~\ref{tab:supp_qcrfsc_empirical}.",
    "table_s23_identity_retrieval.tex": r"Cross-batch identity discrimination. Values are query-weighted estimates and paired biological-identity bootstrap percentile intervals (2,000 draws). Top-1 accuracy requires the most similar profile in another batch to have the correct identity; MRR is mean reciprocal rank. The last column is WiNN minus the indicated comparator.",
    "table_s24_paired_feature_effects.tex": r"Same-feature effects relative to input abundance. Slope is the through-origin least-squares slope of corrected on input feature effects. Sign concordance and Pearson correlation use matched finite features. Effects use the log1p scale and the dataset-specific osmolality adjustment.",
    "table_s25_nearest_baselines.tex": r"Reference precision and batch association after study-sample regression. Wehrens ANCOVA and mean-only fits use all injections, acquisition order and batch labels. Study-only ANCOVA additionally uses reference identities. Batch $R^2$ is categorical weighted-PC $R^2$ on study injections."
}
for name, caption in captions.items():
    path = ROOT / "generated/supplementary_tables" / name
    text = path.read_text()
    text = re.sub(r"\\caption\{.*?\}\\label", lambda _: r"\caption{" + caption + r"}\label", text, count=1)
    path.write_text(text)

caption_replacements = {
    "frozen point estimates": "point estimates",
    "A dash marks a correction reused from a validated frozen run, for which matrix-loading time was not reported as algorithm runtime.":
        "A dash indicates an unavailable algorithm runtime for a reused correction matrix.",
    "Recorded numerical safeguards. The trigger was retained as part of the benchmark outcome and no affected run was silently dropped.":
        "Numerical safeguard triggers and affected analyses.",
    " No relaxed model was substituted.": "",
    "same frozen settings": "same prespecified settings",
    "Columns are immutable matrices after each ordered stage.": "Columns describe the matrix after each ordered stage.",
    "Complete reference-split stability summary across all datasets and methods.": "Reference-split stability across all datasets and eight outputs.",
    "Complete all-method summary of leave-one-batch-out biological-effect replication.": "Leave-one-batch-out biological-effect replication across eight outputs.",
    "Complete release-validation checklist for the frozen benchmark.": "Validation checklist for the current benchmark release. QC-RFSC updates are documented by the separate integration checks in the reproducibility package."
}
for path in (ROOT / "generated/supplementary_tables").glob("*.tex"):
    text = path.read_text()
    for old_text, new_text in caption_replacements.items():
        text = text.replace(old_text, new_text)
    path.write_text(text)

# Publication labels are distinct from stable analysis method IDs.
display_names = {
    "WiNN auto-batch (QC)": "WiNN (inferred segments)",
    "WiNN auto (QC)": "WiNN (QC-tuned)",
    "WiNN auto-batch": "WiNN (inferred segments)",
    "WiNN auto": "WiNN (QC-tuned)",
    "WiNN fixed default (no QC)": "WiNN",
    "WiNN fixed": "WiNN",
    "Fixed-WiNN": "WiNN", "fixed-WiNN": "WiNN",
    "Fixed WiNN": "WiNN", "fixed WiNN": "WiNN",
}
for directory, extension in (("supplementary_tables", "tex"), ("supplementary_table_data", "csv")):
    for path in (ROOT / "generated" / directory).glob("*." + extension):
        content = path.read_bytes()
        for previous, current in display_names.items():
            content = content.replace(previous.encode(), current.encode())
        path.write_bytes(content)

write(VALIDATION / "qcrfsc_completed_integration_checks.csv", checks)
print(f"Validated {len(checks)} checks; integrated {len(corrected)} simulation metrics.")
