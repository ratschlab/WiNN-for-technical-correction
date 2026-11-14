#!/usr/bin/env python3
"""Assemble the condensed supplementary tables.

The detailed generators write one table per analysis into generated/supplementary_tables/
and generated/supplementary_table_data/. Bioinformatics asks that supplementary material
contain only directly relevant content, so the supplement merges related tables:

  condensed_selected_settings.tex   comparator settings (S21) + WiNN parameters (S17)
  condensed_complete_benchmark.tex  simulation and five empirical benchmarks (S4-S9),
                                    with the corrected QC-RFSC empirical points (S26)
  condensed_selective_forced.tex    compact pivot of the gate-policy experiment (S15)
  condensed_baselines.tex           same-information baselines (S27, S25, S28)

Tables that were dropped from the supplement (reference splits, cumulative stages,
leave-one-batch-out summaries, safeguards, validation checklist, reference placement)
remain available as CSV source data. Run from the manuscript directory after the
table generators:

  python3 scripts/condense_supplementary_tables.py
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEX = ROOT / "generated/supplementary_tables"
DATA = ROOT / "generated/supplementary_table_data"
SOURCE = ROOT / "figures/source_data"

DATASETS = ["Simulation", "MTBLS79", "Clinical FIA-MS", "BatchCorr Set 1", "Sacurine", "WaveICA"]


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def body_rows(name):
    """Data rows of a generated longtable, split into cells."""
    text = (TEX / name).read_text()
    body = text.split("\\endlastfoot", 1)[1].split("\\end{longtable}", 1)[0]
    rows = []
    for line in body.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("\\midrule") or line.startswith("\\multicolumn"):
            rows.append(line)
            continue
        rows.append([cell.strip() for cell in re.sub(r"\\\\\s*$", "", line).split(" & ")])
    return rows


def row(cells):
    return " & ".join(cells) + " \\\\"


def longtable(spec, caption, label, header, lines, size="\\scriptsize", landscape=True):
    """Landscape or portrait longtable; header is a list of cells, a preformatted row or None."""
    ncol = len(re.sub(r"\{[^}]*\}", "", spec))
    if header is None:
        head = []
    else:
        head = [header if isinstance(header, str) else row(header), "\\midrule"]
    out = []
    if landscape:
        out.append("\\begin{landscape}")
    out += [size, f"\\begin{{longtable}}{{{spec}}}",
            f"\\caption{{{caption}}}\\label{{{label}}}\\\\",
            "\\toprule", *head, "\\endfirsthead",
            f"\\multicolumn{{{ncol}}}{{l}}{{\\textit{{Table \\thetable{{}} continued}}}}\\\\",
            "\\toprule", *head, "\\endhead",
            "\\midrule", f"\\multicolumn{{{ncol}}}{{r}}{{\\textit{{Continued on next page}}}}\\\\", "\\endfoot",
            "\\bottomrule", "\\endlastfoot"]
    out += lines
    out.append("\\end{longtable}")
    if landscape:
        out.append("\\end{landscape}")
    return "\n".join(out) + "\n"


def selected_settings():
    comparators = [r for r in body_rows("table_s21_selected_parameters.tex") if isinstance(r, list)]
    winn = read_csv(DATA / "table_s17_winn_parameters.csv")
    pqn = {"shrink": "conservative", "normalize": "full"}
    lines = []
    for dataset in DATASETS:
        for cells in comparators:
            if cells[0] == dataset:
                lines.append(row(cells))
        for r in winn:
            if r["Dataset"] != dataset:
                continue
            source = "inferred" if r["BatchSource"].startswith("auto") else "supplied"
            drift_error = "FWER" if "maxT" in r["Test"] else "FDR"
            setting = (f"batches {source} ({r['Segments']}); {r['Test']}; {r['Smoother']} smoother; "
                       f"drift {drift_error} {r['DriftFDR']}; batch FDR {r['BatchFDR']}; PQN {pqn.get(r['PQN'], r['PQN'])}")
            tuned = r["Mode"] != "WiNN"
            lines.append(row([dataset, r["Mode"], setting, "grid" if tuned else "1",
                              "Package training-QC score" if tuned else "Prespecified"]))
        lines.append("\\midrule")
    lines = lines[:-1]
    caption = ("Selected correction settings for every dataset and method. Comparator candidates were "
               "prespecified or scored on labelled training QCs only. WiNN (QC-tuned) and WiNN (inferred "
               "segments) selected settings from the package grid (36 nominal candidates; Section~S1) using "
               "training QCs only; segment counts are shown in parentheses. WiNN uses one prespecified "
               "configuration in every dataset and receives no QC identities. Raw is the identity operation, "
               "and TIGER uses its native $4\\times4$ ensemble with injection order and 500 trees in every dataset.")
    return longtable("p{2.6cm}p{3.4cm}p{9.6cm}rp{3.6cm}", caption, "tab:supp_selected_parameters",
                     ["Dataset", "Method", "Selected settings", "Candidates", "Selection"], lines)


def decimals(value):
    match = re.match(r"-?\d+\.(\d+)", value.strip())
    return len(match.group(1)) if match else 0


def complete_benchmark():
    rfsc = {r["dataset"]: r for r in read_csv(DATA / "table_s26_qcrfsc_empirical.csv")}
    rfsc_mtbls79 = next(r for r in read_csv(ROOT / "generated/review_validation/qcrfsc_mtbls79_sample_means.csv")
                       if r["method_id"] == "qc_rfsc")
    assert int(rfsc_mtbls79["biological_n"]) == 20
    rfsc_waveica = next(r for r in read_csv(ROOT / "generated/review_validation/qcrfsc_completed_empirical_all_endpoints.csv")
                       if r["dataset_key"] == "waveica" and r["metric"] == "biological_weighted_pc_r2")
    assert int(rfsc_waveica["denominator"]) == 568
    blocks = [
        ("Simulation", "table_s4_simulation.csv", "mean $\\pm$ SE, ten seeds",
         ["Batch R2", "QC CV", "GAM", "RMSE", "ICC", "Pearson"], ["RMSE", "ICC(A,1)", "Sample $r$"]),
        ("MTBLS79", "table_s5_mtbls79.csv", None,
         ["QC CV", "GAM", "Batch R2", "Cow-sheep R2 [95% CI]", "Replicate Pearson", "Feature ICC"],
         ["Cow-sheep $R^2$ [95\\% CI]", "Replicate $r$", "Batch-pair ICC"]),
        ("BatchCorr Set 1", "table_s6_batchcorr.csv", None,
         ["QC CV", "GAM", "Batch R2", "Replicate Pearson [95% CI]", "Replicate ICC [95% CI]"],
         ["Replicate $r$ [95\\% CI]", "Replicate ICC [95\\% CI]", ""]),
        ("Clinical FIA-MS", "table_s7_clinical.csv", "plate $R^2$",
         ["QC CV", "GAM", "Batch R2", "Replicate Pearson [95% CI]", "Replicate ICC [95% CI]"],
         ["Replicate $r$ [95\\% CI]", "Replicate ICC [95\\% CI]", ""]),
        ("Sacurine", "table_s8_sacurine.csv", None,
         ["QC CV", "GAM", "Batch R2", "Age R2 [95% CI]", "BMI R2 [95% CI]", "Gender R2 [95% CI]"],
         ["Age $R^2$ [95\\% CI]", "BMI $R^2$ [95\\% CI]", "Gender $R^2$ [95\\% CI]"]),
        ("WaveICA", "table_s9_waveica.csv", None,
         ["QC CV", "GAM", "Batch R2", "Group R2 [95% CI]", "BH features", "Standardized effect"],
         ["Group $R^2$ [95\\% CI]", "BH features", "Std.\\ effect"]),
    ]
    lines = []
    for dataset, filename, note, columns, labels in blocks:
        rows = read_csv(DATA / filename)
        title = f"\\textbf{{{dataset}}}" + (f" ({note})" if note else "")
        lines += ["\\midrule", f"\\multicolumn{{4}}{{l}}{{{title}}} & " + " & ".join(f"\\textit{{{x}}}" if x else "" for x in labels) + " \\\\", "\\midrule"]
        if dataset == "Simulation":
            order = ["QC CV", "GAM", "Batch R2", "RMSE", "ICC", "Pearson"]
        else:
            order = columns + [""] * (6 - len(columns))
        for r in rows:
            method = r.get("Method") or r.get("method")
            cells = [method] + [r[c] if c else "" for c in order]
            lines.append(row(cells))
            if dataset != "Simulation" and method == "QC-RLSC" and dataset in rfsc:
                q = rfsc[dataset]
                raw = rows[0]
                cv = f"{float(q['qcrfsc_cv_percent']):.2f}"
                gam = f"{float(q['residual_gam_percent']) / 100:.{decimals(raw['GAM'])}f}"
                batch = f"{float(q['batch_weighted_pc_r2']):.{decimals(raw['Batch R2'])}f}"
                extra = ["--", "--", "--"]
                if dataset in ("BatchCorr Set 1", "Clinical FIA-MS"):
                    extra = [f"{float(q['replicate_profile_pearson']):.4f}", f"{float(q['replicate_profile_icc']):.4f}", ""]
                elif dataset == "MTBLS79":
                    extra = [f"{float(rfsc_mtbls79['biological_pc_variance_percent']) / 100:.3f}",
                             f"{float(q['replicate_profile_pearson']):.3f}", "--"]
                elif dataset == "WaveICA":
                    extra = [f"{float(rfsc_waveica['value']):.4f}", "--", "--"]
                lines.append(row(["QC-RFSC$^{a}$", cv, gam, batch] + extra))
    lines = lines[1:]
    caption = ("Complete benchmark results for all outputs. QC CV is mean held-out reference CV (\\%); GAM is "
               "residual run-order deviance explained (proportion); batch $R^2$ is categorical weighted-PC $R^2$ "
               "among study injections. Simulation cells are mean $\\pm$ seed-level standard error and include "
               "clean-truth endpoints. Empirical cells are point estimates; brackets are 95\\% percentile intervals "
               "from 2,000 paired bootstrap resamples of the independent biological unit, conditional on the fitted "
               "corrections. Replicate $r$ and ICC are replicate-profile Pearson correlation and within-identity "
               "feature-profile ICC(A,1); BH features counts group-associated features after Benjamini--Hochberg "
               "adjustment. $^{a}$QC-RFSC empirical rows are point estimates from the corrected adapter, whose "
               "resampling analyses were not repeated; 14, 37, 1, 2 and 459 features in MTBLS79, BatchCorr, the "
               "clinical cohort, Sacurine and WaveICA kept their input values after nonfinite fits.")
    header = "Method & QC CV (\\%) & GAM & Batch $R^2$ & \\multicolumn{3}{c}{Preservation endpoints} \\\\"
    return longtable("lrrrrrr", caption, "tab:supp_complete_benchmark", header, lines)


def selective_forced():
    rows = read_csv(DATA / "table_s15_selective_forced.csv")
    value = {(r["Dataset"], r["Policy"], r["Metric"]): r["Value"] for r in rows}
    policies = ["selective drift / selective batch", "forced drift / selective batch",
                "selective drift / forced batch", "forced drift / forced batch"]
    metrics = [("held-out QC CV", "Held-out QC CV (\\%)"), ("residual GAM deviance", "Residual GAM deviance"),
               ("categorical batch R2", "Batch $R^2$")]
    lines = []
    for dataset in DATASETS:
        label = "Simulation (seed 1)" if dataset == "Simulation" else dataset
        for key, name in metrics:
            lines.append(row([label, name] + [value[(dataset, p, key)] for p in policies]))
            label = ""
        if dataset != DATASETS[-1]:
            lines.append("\\midrule")
    # Ten-seed clean-truth endpoints from the Figure 2A / S2D-F source data.
    lines.append("\\midrule")
    names = {"Selective / selective": 0, "Forced drift": 1, "Forced batch": 2, "Forced drift + batch": 3}
    first = True
    for stem, name, digits in (("S2D_truth_log1p_rmse", "Clean-truth RMSE", 4),
                               ("S2E_truth_metabolite_profile_icc_mean", "Feature ICC(A,1)", 4),
                               ("S2F_truth_sample_profile_pearson_mean", "Sample correlation", 4)):
        means = [""] * 4
        for r in read_csv(SOURCE / f"{stem}.csv"):
            means[names[r["gate_label"]]] = f"{float(r['method_mean']):.{digits}f}"
        lines.append(row(["Simulation (ten seeds)" if first else "", name] + means))
        first = False
    caption = ("Selective versus forced correction. Each policy changes one gate of the WiNN workflow: forced "
               "drift detrends every testable feature-segment profile, forced batch adjusts every eligible "
               "feature, and MAD preparation and PQN remain active. Technical rows use one prespecified "
               "representative simulation seed and the empirical point estimates; clean-truth rows are means over "
               "ten simulation seeds. "
               "Lower values are favorable for CV, GAM, batch $R^2$ and RMSE.")
    return longtable("llrrrr", caption, "tab:supp_selective_forced",
                     ["Dataset", "Endpoint", "Selective", "Forced drift", "Forced batch", "Forced both"],
                     lines, size="\\footnotesize", landscape=False)


def baselines():
    truth = body_rows("table_s27_baseline_simulation.tex")
    precision = [r for r in body_rows("table_s25_nearest_baselines.tex") if isinstance(r, list)]
    effects = [r for r in body_rows("table_s28_baseline_empirical_effects.tex") if isinstance(r, list)]
    lines = ["\\multicolumn{7}{l}{\\textbf{A. Simulation truth (ten seeds; mean $\\pm$ SE)}} \\\\",
             row(["", "\\textit{Method}", "\\textit{Hidden-ref.\\ CV (\\%)}", "\\textit{RMSE}",
                  "\\textit{Feature ICC}", "\\textit{Sample $r$}", ""]), "\\midrule"]
    advantage = False
    for r in truth:
        if isinstance(r, list):
            if advantage:
                lines.append(row(["WiNN advantage", "over " + r[0]] + r[1:] + [""]))
            else:
                lines.append(row(["Simulation", r[0]] + r[1:] + [""]))
        elif r.startswith("\\multicolumn"):
            advantage = True
    lines += ["\\midrule", "\\multicolumn{7}{l}{\\textbf{B. Empirical reference precision and batch association}} \\\\",
              row(["\\textit{Dataset}", "\\textit{Fit}", "\\textit{Held-out QC CV (\\%)}", "\\textit{Batch $R^2$}",
                   "\\textit{Floored cells}", "", ""]), "\\midrule"]
    lines += [row(r + ["", ""]) for r in precision]
    lines += ["\\midrule", "\\multicolumn{7}{l}{\\textbf{C. Matched-feature effects relative to raw input}} \\\\",
              row(["\\textit{Dataset / endpoint}", "\\textit{Method}", "\\textit{Slope}", "\\textit{Sign agreement}",
                   "\\textit{Pearson}", "\\textit{LOBO Pearson}", "\\textit{BH count}"]), "\\midrule"]
    lines += [row(r) for r in effects]
    caption = ("Same-information baselines. Mean-only batch centering and the Wehrens study-sample ANCOVA use all "
               "injections with the same intensity, acquisition-order and batch information as WiNN; the study-only "
               "ANCOVA additionally uses reference identities to choose fitting observations. (A) Ten-seed clean-truth "
               "endpoints; WiNN advantages are comparator minus WiNN for RMSE and WiNN minus comparator otherwise, "
               "with 95\\% $t_9$ intervals, and WiNN improved every endpoint in every seed. (B) Held-out reference CV, "
               "categorical batch weighted-PC $R^2$ among study injections and negative extrapolations floored at zero. "
               "(C) Through-origin slope, sign agreement and Pearson correlation of corrected against raw feature "
               "effects, median leave-one-batch-out effect correlation and BH-significant feature counts. Raw effects "
               "are an observed reference, not biological truth.")
    return longtable("p{3.0cm}p{3.6cm}rrrrr", caption, "tab:supp_baselines", None, lines, size="\\footnotesize")


def main():
    outputs = {
        "condensed_selected_settings.tex": selected_settings(),
        "condensed_complete_benchmark.tex": complete_benchmark(),
        "condensed_selective_forced.tex": selective_forced(),
        "condensed_baselines.tex": baselines(),
    }
    for name, text in outputs.items():
        (TEX / name).write_text(text)
        print("wrote", (TEX / name).relative_to(ROOT))


if __name__ == "__main__":
    main()
    from normalize_table_terminology import main as normalize_tables
    normalize_tables()
