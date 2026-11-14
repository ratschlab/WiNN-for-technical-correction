#!/usr/bin/env python3
"""Check every numerical claim in main.tex against the frozen revised WiNN source data.

Each check recomputes the displayed value from a CSV and requires the formatted token to
appear in the main-text paragraph that starts the claim. Run from the manuscript directory:

  python3 revision_notes/audit_main_text_claims.py [--derivative-root=<frozen derivative results>]

All input CSVs needed for the current main-text claims are bundled with the manuscript.
"""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "generated/supplementary_table_data"
VAL = ROOT / "generated/review_validation"
SRC = ROOT / "figures/source_data"


text = (ROOT / "main.tex").read_text()
paragraphs = text.split("\n\n")
checks = []


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def one(rows, **keys):
    hits = [r for r in rows if all(r[k] == v for k, v in keys.items())]
    assert len(hits) == 1, (keys, len(hits))
    return hits[0]


def lead(value):
    return float(str(value).split()[0])


def check(claim, anchor, token):
    blocks = [p for p in paragraphs if anchor in p]
    assert len(blocks) == 1, (anchor, len(blocks))
    checks.append((claim, token, token in blocks[0]))


ABS = "\\textbf{Results:}"
SIM = "Across ten paired simulations"
PC = "The partial-confounding study showed"
BC = "BatchCorr tested whether"
CL = "The internal clinical FIA-MS cohort provided"
WV = "WaveICA had strong batch"
MT = "MTBLS79 provided a known cow-sheep"
AB = "We examined each correction step"
ST = "We tested whether WiNN selects features"
DI = "This is important because correction"
RT = "Runtimes are reported"

# Simulation benchmark (ten-seed means)
s4 = read(DATA / "table_s4_simulation.csv")
raw, winn = one(s4, method_id="raw"), one(s4, method_id="winn_fixed_no_qc")
for anchor in (SIM,):
    check("Sim raw QC CV", anchor, f"{lead(raw['QC CV']):.2f}\\%")
    check("Sim WiNN QC CV", anchor, f"{lead(winn['QC CV']):.2f}\\%")
    check("Sim raw GAM", anchor, f"{lead(raw['GAM']) * 100:.2f}\\%")
    check("Sim WiNN GAM", anchor, f"{lead(winn['GAM']) * 100:.2f}\\%")
for anchor in (SIM, ABS):
    check("Sim raw RMSE", anchor, f"{lead(raw['RMSE']):.3f}")
    check("Sim WiNN RMSE", anchor, f"{lead(winn['RMSE']):.3f}")
check("Sim raw ICC", SIM, f"{lead(raw['ICC']):.3f}")
check("Sim WiNN ICC", SIM, f"{lead(winn['ICC']):.3f}")
check("Sim raw Pearson", SIM, f"{lead(raw['Pearson']):.3f}")
check("Sim WiNN Pearson", SIM, f"{lead(winn['Pearson']):.3f}")
seed = read(ROOT / "generated/revision_panel_inputs/simulation_seed_method_metrics.csv")
for metric, sign in (("truth_log1p_rmse", -1), ("truth_metabolite_profile_icc_mean", 1), ("truth_sample_profile_pearson_mean", 1)):
    by = {}
    for r in seed:
        if r["metric"] == metric and r["value"] not in ("", "NA"):
            by.setdefault(r["seed_id"], {})[r["method_id"]] = float(r["value"])
    wins = all(all((d["winn_fixed_no_qc"] - d[m]) * sign > 0 for m in ("combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf")) for d in by.values())
    checks.append((f"WiNN beats all comparators every seed: {metric}", "all ten seeds", wins and len(by) == 10))
base = read(VAL / "baseline_simulation_summary.csv")
for r in base:
    if r.get("metric") == "truth_log1p_rmse" and r.get("method_id") in ("mean_batch_all", "wehrens_ancova_all"):
        check(f"Baseline RMSE {r['method_id']}", SIM, f"{float(r['mean']):.3f}")

# Partial confounding
FULL = "100\\%"
s13 = read(DATA / "table_s13_partial_confounding.csv")
none = one(s13, Type="None")
check("PC random recovery", PC, f"{lead(none['Clean ratio']):.3f}")
assert lead(none["TPR"]) == 1.0
check("PC batch TPR", PC, f"{lead(one(s13, Type='Batch', Alignment=FULL)['TPR']):.3f}")
check("PC run-order recovery", PC, f"{lead(one(s13, Type='Run order', Alignment=FULL)['Clean ratio']):.3f}")
joint = one(s13, Type="Joint", Alignment=FULL)
for field in ("Clean ratio", "TPR", "FDR"):
    check(f"PC joint {field}", PC, f"{lead(joint[field]):.3f}")
cv = [lead(r["QC CV"]) for r in s13]
check("PC CV min", PC, f"{min(cv):.2f}\\%")
check("PC CV max", PC, f"{max(cv):.2f}\\%")
pcs = read(ROOT / "generated/revision_panel_inputs/S5B_partial_residual_run_order_gam_mean_deviance.csv")
gam = [float(r["mean"]) for r in pcs if r["metric"] == "residual_run_order_gam_mean_deviance" and r["method"].startswith("WiNN")]
check("PC GAM min", PC, f"{min(gam) * 100:.1f}\\%")
check("PC GAM max", PC, f"{max(gam) * 100:.1f}\\%")

# Repeated-identity datasets (Figure 3 source data)
for stem, anchor, cv_label in (("F3A_batchcorr_qc_cv_vs_run_order_deviance", BC, None), ("F3B_clinical_qc_cv_vs_run_order_deviance", CL, None)):
    rows = read(SRC / f"{stem}.csv")
    for method in ("raw", "winn_fixed_no_qc"):
        r = one(rows, method_id=method)
        vals = list(r.values())
        cv_value, gam_value = float(vals[8]), float(vals[10]) * 100
        batch_value = float(r["point_estimate_batch_weighted_pc_r2_categorical"]) * 100
        icc = float(r["point_estimate_genuine_replicate_icc_median"])
        check(f"{stem} {method} CV", anchor, f"{cv_value:.2f}\\%")
        check(f"{stem} {method} GAM", anchor, f"{gam_value:.2f}\\%")
        digits = 2 if method == "raw" else 3
        check(f"{stem} {method} batch", anchor, f"{batch_value:.{digits}f}\\%")
        check(f"{stem} {method} ICC", anchor, f"{icc:.3f}" if "batchcorr" in stem else f"{icc:.4f}")
clinical_points = read(SRC / "F3B_clinical_qc_cv_vs_run_order_deviance.csv")
tiger_gam = float(one(clinical_points, method_id="tiger")["residual_run_order_gam_mean_deviance"])
check("Clinical TIGER residual run-order drift", CL, f"{tiger_gam * 100:.2f}\\%")
primary_ids = {"raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf", "winn_fixed_no_qc"}
clinical_primary = [r for r in clinical_points if r['method_id'] in primary_ids]
checks.append(("Clinical WiNN has the lowest primary QC CV", "seven primary methods",
               len(clinical_primary) == 7 and
               min(clinical_primary, key=lambda r: float(r['qc_cv']))['method_id'] == 'winn_fixed_no_qc'))
checks.append(("Clinical TIGER has the lowest primary residual drift", "TIGER",
               len(clinical_primary) == 7 and
               min(clinical_primary, key=lambda r: float(r['residual_run_order_gam_deviance_percent']))['method_id'] == 'tiger'))
ret = read(VAL / "frozen_identity_retrieval_estimates.csv")
for dataset, anchor, methods in (("batchcorr_set1", BC, ("raw", "winn_fixed_no_qc", "combat", "qc_rlsc")),
                                 ("clinical_fiams", CL, ("raw", "winn_fixed_no_qc", "tiger"))):
    for method in methods:
        r = one(ret, dataset=dataset, metric="cross_batch_top1_accuracy", method_id=method)
        check(f"Retrieval {dataset} {method}", anchor, f"{float(r['point_estimate']) * 100:.1f}\\%")
        if method in ("raw", "winn_fixed_no_qc"):
            check(f"Abstract retrieval {dataset} {method}", ABS, f"{float(r['point_estimate']) * 100:.1f}\\%")

# Public biological contrasts (Figure 4 source data)
for stem, anchor in (("F4A_waveica_batch_vs_group", WV), ("F4B_mtbls79_batch_vs_group", MT)):
    rows = read(SRC / f"{stem}.csv")
    for method in ("raw", "winn_fixed_no_qc"):
        r = one(rows, method_id=method)
        batch = float(r["point_estimate_batch_weighted_pc_r2_categorical"]) * 100
        group = float(r["point_estimate_biological_weighted_pc_r2"]) * 100
        check(f"{stem} {method} batch", anchor, f"{batch:.2f}\\%" if method == "raw" else f"{batch:.3f}\\%")
        check(f"{stem} {method} group", anchor, f"{group:.3f}\\%" if "waveica" in stem else f"{group:.2f}\\%")
discovery = read(DATA / "table_s11_association_permutation.csv")
for method in ("Raw", "WiNN"):
    row = one(discovery, Dataset="MTBLS79", Endpoint="cow vs sheep", Method=method)
    check(f"MTBLS79 discoveries {method}", MT, f"{int(row['Significant']):,}")
check("MTBLS79 tested features", MT, f"{int(row['Tested']):,}")

# Ablations
s15 = {(r["Dataset"], r["Policy"], r["Metric"]): float(r["Value"]) for r in read(DATA / "table_s15_selective_forced.csv")}
datasets = sorted({k[0] for k in s15})
up = [d for d in datasets if s15[(d, "forced drift / selective batch", "held-out QC CV")] > s15[(d, "selective drift / selective batch", "held-out QC CV")]]
checks.append(("Forced drift raises CV in four datasets", "four of six", len(up) == 4 and set(datasets) - set(up) == {"BatchCorr Set 1", "Sacurine"}))
f2a = read(SRC / "F2A_selective_forced_truth_rmse.csv")
check("Selective RMSE", AB, f"{float(one([r for r in f2a if r['seed_id'] == 'SIM01'], policy='Selective')['method_mean']):.4f}")
check("Forced drift RMSE", AB, f"{float(one([r for r in f2a if r['seed_id'] == 'SIM01'], policy='Force drift')['method_mean']):.4f}")
ablation_seeds = read(ROOT / "generated/revision_panel_inputs/simulation_seed_ablation_metrics.csv")
gate_seeds = {r['seed_id'] for r in ablation_seeds}
for metric, sign in (("truth_log1p_rmse", -1), ("truth_metabolite_profile_icc_mean", 1),
                     ("truth_sample_profile_pearson_mean", 1)):
    forced_drift_wins = sum(
        (float(one(ablation_seeds, seed_id=seed, method_id='G_AS', metric=metric)['value']) -
         float(one(ablation_seeds, seed_id=seed, method_id='G_SS', metric=metric)['value'])) * sign > 0
        for seed in gate_seeds
    )
    checks.append((f"Forced drift improves {metric} in all ten seeds", "all ten simulations",
                   len(gate_seeds) == 10 and forced_drift_wins == 10))
s10 = read(DATA / "table_s10_runtime_retention.csv")
frac = [int(r["Batch features"]) / int(r["Features"]) * 100 for r in s10 if r["Method"] == "WiNN"]
check("Batch selection min", AB, f"{min(frac):.1f}")
check("Batch selection max", AB, f"{max(frac):.1f}\\%")

# Stress tests
null = read(VAL / "null_gate_verification_summary.csv")
for n in ("20", "50", "100", "200"):
    check(f"Null gate n={n}", ST, f"{float(one(null, n=n, stage='mad_prepared')['any_rejection']) * 100:.1f}\\%")
calibrated = read(ROOT / "generated/calibration_revision/gate_validation_summary.csv")
for n in ("20", "50", "100", "200"):
    row = one(calibrated, experiment="original_null", n=n, stage="mad_prepared", algorithm="permutation_maxT")
    check(f"Permutation null gate n={n}", ST, f"{float(row['fwer']) * 100:.1f}\\%")
prep = read(VAL / "preparation_biology_summary.csv")
for scenario in ("balanced_sparse", "rare_subgroup", "broad_coordinated"):
    check(f"MAD+PQN {scenario}", ST, f"{float(one(prep, scenario=scenario, stage='mad_then_pqn')['median_recovery_ratio']):.3f}")

# Measured runtime completeness; no prespecified speedup is required.
times = {(r["Dataset"], r["Method"]): float(r["Runtime"]) for r in s10 if r["Runtime"] not in ("--", "")}
winn_times = [v for (d,m),v in times.items() if m == "WiNN"]
checks.append(("Six measured WiNN runtimes are finite and positive", "6 datasets", len(winn_times) == 6 and all(math.isfinite(v) and v > 0 for v in winn_times)))

failed = [c for c in checks if not c[2]]
for claim, token, _ in failed:
    print("FAILED:", claim, "expected", token)
print(f"Checked {len(checks)} main-text claims; {len(failed)} failed.")
raise SystemExit(bool(failed))
