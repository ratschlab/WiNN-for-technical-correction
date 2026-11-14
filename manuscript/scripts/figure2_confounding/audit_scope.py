#!/usr/bin/env python3
"""Verify immutable non-confounding results against the pre-revision baseline."""
import csv,hashlib,json,sys,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
path=ROOT/'generated/figure2_confounding/protected_baseline.json'
ref=json.loads(path.read_text())
# The later supplement finalization corrected only stale S7 rendering inputs.
# Explicit opt-in retains the original confounding-only audit.
allowed = {}
if '--allow-finalized-s7' in sys.argv:
    subprocess.run([sys.executable, str(ROOT/'scripts/audit_supplement_finalization.py')], check=True)
    revision = json.loads((ROOT/'generated/review_validation/supplement_finalization_provenance.json').read_text())
    allowed = {r['path']: r['sha256'] for r in revision['updated_s7_files']}
    assert len(allowed) == 13
    assert all('S7' in p or p in ('generated/supplementary_table_data/leave_one_batch_out_bootstrap_summary.csv', 'generated/supplementary_table_data/README.md') for p in allowed)
changed = []
for r in ref['files']:
    current = hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()
    if current != r['sha256']:
        assert allowed.get(r['path']) == current, r['path']
        changed.append(r['path'])
blocks=(ROOT/'main.tex').read_text().split('\n\n')
for r in ref['main_text_blocks']:
    hits=[b for b in blocks if r['anchor'] in b]
    assert len(hits)==1 and hashlib.sha256(hits[0].encode()).hexdigest()==r['sha256'],r['anchor']
new=list(csv.DictReader((ROOT/'generated/main_figure_assembly_manifest.csv').open()))
for row in ref['assembly']:
    match=next(r for r in new if r['panel']==row['panel'])
    for key in row:
        assert row[key]==match[key] or (key in ('width_mm','height_mm','x_mm','y_mm') and abs(float(row[key])-float(match[key]))<1e-12),(key,row,match)
report=dict(baseline_commit=ref['baseline_commit'],protected_files=len(ref['files']),all_identical=not changed,
    unchanged_files=len(ref['files'])-len(changed),approved_s7_finalization_changes=changed,
    figure_2_A_to_D_placement_unchanged=True,protected_main_text_blocks=len(ref["main_text_blocks"]),baseline_manifest_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(ROOT/'generated/figure2_confounding/scope_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(f"PASS: {len(ref['files'])-len(changed)} protected files byte-identical; {len(changed)} approved S7 finalization changes; Figure 2A–D placement unchanged.")
