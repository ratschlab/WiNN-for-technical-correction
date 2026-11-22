#!/usr/bin/env python3
"""Rebuild current RFSC table overlays without loading superseded fit results."""
from pathlib import Path
import csv, re, subprocess, sys
ROOT=Path(__file__).resolve().parents[2]
V=ROOT/'generated/review_validation'
def read(path):
    with path.open(newline='') as f: return list(csv.DictReader(f))
def write(path, rows, fields=None):
    if fields is None: fields=list(rows[0])
    with path.open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=fields);out.writeheader();out.writerows(rows)
# Fail before changing tables if the current release cannot be reconstructed.
for ds in ['mtbls79','batchcorr_set1','clinical_fiams','sacurine','waveica']:
    for folder in ['results','empirical_evaluation']:
        p=V/'chronology_empirical_20261122'/folder/ds/'COMPLETE.txt'
        if not p.exists() or 'status=complete' not in p.read_text(): raise RuntimeError('Missing current chronology output: '+str(p))
    p=V/'qcrfsc_bootstrap_20261122'/ds/'COMPLETE.txt'
    if not p.exists(): raise RuntimeError('Missing current uncertainty output: '+str(p))
    provenance=read(p.with_name('input_manifest.csv'))
    matrix=next(r for r in provenance if r['object']=='corrected_matrix')
    current=read(V/'chronology_empirical_20261122'/'results'/ds/'chronology_matrix_audit.csv')[0]
    if matrix['sha256']!=current['new_matrix_file_sha256'] or int(matrix['nboot'])!=2000:
        raise RuntimeError('Uncertainty does not belong to the current corrected matrix: '+ds)
for i in range(1,11):
    for folder in ['fits','evaluation']:
        p=V/'qcrfsc_chronology_simulation'/folder/f'SIM{i:02d}'/'COMPLETE.txt'
        if not p.exists() or 'status=complete' not in p.read_text(): raise RuntimeError('Missing current simulation output: '+str(p))
if '--check-only' in sys.argv:
    print('PASS: all current chronology fits and matching 2000-draw uncertainty inputs are complete.')
    sys.exit(0)
# Archive-derived RFSC rows are excluded where that analysis has not been rerun.
# Current matching and other completed endpoint rows are restored below.
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

captions = {
    "table_s22_primary_safeguards.tex": r"Numerical safeguards in the benchmark fits. The table records negative intensities floored at zero and features restored to input abundance. QC-RFSC fallback counts are reported in Table~\ref{tab:supp_qcrfsc_empirical}.",
    "table_s23_identity_retrieval.tex": r"Cross-batch identity discrimination. Values are query-weighted estimates and paired biological-identity bootstrap percentile intervals (2,000 draws). Top-1 accuracy requires the most similar profile in another batch to have the correct identity; MRR is mean reciprocal rank. The last column is WiNN minus the indicated comparator.",
    "table_s24_paired_feature_effects.tex": r"Same-feature effects relative to input abundance. Slope is the through-origin least-squares slope of corrected on input feature effects. Sign concordance and Pearson correlation use matched finite features. Effects use the log1p scale and the dataset-specific osmolality adjustment.",
    "table_s25_nearest_baselines.tex": r"Reference precision and batch association. All-injection ANCOVA uses all injections, acquisition order and batch labels; mean-only batch adjustment uses all injections and batch labels without acquisition order. Study-only ANCOVA additionally uses reference identities. Batch $R^2$ is categorical weighted-PC $R^2$ on study injections."
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
    "Complete release-validation checklist for the frozen benchmark.": "Validation checklist for the original benchmark release. QC-RFSC updates are documented by the separate integration checks in the reproducibility package."
}
for path in (ROOT / "generated/supplementary_tables").glob("*.tex"):
    text = path.read_text()
    for old_text, new_text in caption_replacements.items():
        text = text.replace(old_text, new_text)
    path.write_text(text)

# Publication labels are distinct from stable analysis method IDs.
display_names = {
    "ANCOVA, all": "All-injection ANCOVA", "ANCOVA, study only": "Study-only ANCOVA",
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

for name in ['integrate_chronological_empirical.py','integrate_chronology_simulation.py','integrate_qcrfsc_bootstrap.py']:
    subprocess.run([sys.executable,str(Path(__file__).with_name(name))],check=True)
print('PASS: complete chronological fits and current uncertainty applied; no legacy fit integration executed.')
