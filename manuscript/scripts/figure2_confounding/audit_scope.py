#!/usr/bin/env python3
"""Verify immutable non-confounding results against the pre-revision baseline."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
path=ROOT/'generated/figure2_confounding/protected_baseline.json'
ref=json.loads(path.read_text())
for r in ref['files']:
    assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
blocks=(ROOT/'main.tex').read_text().split('\n\n')
for r in ref['main_text_blocks']:
    hits=[b for b in blocks if r['anchor'] in b]
    assert len(hits)==1 and hashlib.sha256(hits[0].encode()).hexdigest()==r['sha256'],r['anchor']
new=list(csv.DictReader((ROOT/'generated/main_figure_assembly_manifest.csv').open()))
for row in ref['assembly']:
    match=next(r for r in new if r['panel']==row['panel'])
    for key in row:
        assert row[key]==match[key] or (key in ('width_mm','height_mm','x_mm','y_mm') and abs(float(row[key])-float(match[key]))<1e-12),(key,row,match)
report=dict(baseline_commit=ref['baseline_commit'],protected_files=len(ref['files']),all_identical=True,
    figure_2_A_to_D_placement_unchanged=True,protected_main_text_blocks=len(ref["main_text_blocks"]),baseline_manifest_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(ROOT/'generated/figure2_confounding/scope_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(f"PASS: {len(ref['files'])} protected result/figure files byte-identical; Figure 2A–D placement unchanged.")
