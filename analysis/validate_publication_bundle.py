#!/usr/bin/env python3
"""Validate the committed publication data and its complete file manifest."""
from pathlib import Path
import csv
import hashlib
import gzip
import math
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / 'manuscript'
MANIFEST = ROOT / 'publication_manifest.csv'

def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))

def require(test, message):
    if not test:
        raise AssertionError(message)

def tracked_files():
    names = subprocess.check_output(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=ROOT).decode().split('\0')
    return sorted({name for name in names if name and name != MANIFEST.name and (ROOT / name).is_file()})

def validate():
    metrics = read(ROOT / 'results/final/tables/simulation_seed_method_metrics.csv')
    source = read(MANUSCRIPT / 'generated/revision_panel_inputs/simulation_seed_method_metrics.csv')
    observed = {(r['seed_id'], r['method_id'], r['metric']): r['value'] for r in metrics}
    expected = {(r['seed_id'], r['method_id'], r['metric']): r['value'] for r in source}
    require(observed == expected, 'Public and manuscript simulation seed metrics differ')
    summary = read(ROOT / 'results/final/tables/simulation_seed_summary.csv')
    for row in summary:
        values = [float(r['value']) for r in metrics if r['method_id'] == row['method_id'] and r['metric'] == row['metric'] and r['value'] not in ('', 'NA', 'NaN') and math.isfinite(float(r['value']))]
        require(len(values) == int(row['n']), 'Simulation summary count differs')
        require(abs(statistics.mean(values) - float(row['mean'])) < 1e-12, 'Simulation mean differs')
    for row in read(MANUSCRIPT / 'generated/review_validation/qcrfsc_completed_simulation_summary.csv'):
        match = next(r for r in summary if r['method_id'] == 'qc_rfsc' and r['metric'] == row['metric'])
        require(int(match['n']) == 10 and abs(float(match['mean']) - float(row['mean'])) < 1e-12, 'Corrected RFSC summary differs')
    checks = read(MANUSCRIPT / 'generated/review_validation/qcrfsc_completed_integration_checks.csv')
    require(len(checks) >= 225 and all(r['passed'].lower() == 'true' for r in checks), 'RFSC integration checks failed')
    endpoints = read(ROOT / 'results/final/tables/qcrfsc_completed_empirical_all_endpoints.csv')
    require(len(endpoints) == 84 and len({r['dataset_key'] for r in endpoints}) == 5, 'Incomplete empirical endpoint export')
    manifest = read(ROOT / 'analysis/config/endpoint_free_selection_manifest.csv')
    manuscript_manifest = read(MANUSCRIPT / 'generated/supplementary_table_data/endpoint_free_selection_manifest.csv')
    require(manifest == manuscript_manifest, 'Selected-parameter manifests differ')
    for row in read(MANUSCRIPT / 'generated/review_validation/qcrfsc_completed_empirical_summary.csv'):
        selected = next(r for r in manifest if r['dataset_key'] == row['dataset_key'] and r['method_id'] == 'qc_rfsc')
        require(selected['selected_parameters'] == row['selected_parameters'], 'Selected RFSC configuration differs')
    for f in (ROOT / 'analysis/config/revision_20261027').glob('*.csv'):
        if 'manifest' in f.name:
            require(all('Clinical' not in r.get('dataset', '') for r in read(f)), 'Private clinical fold identity present')
    # Clinical aggregate tables must contain no identity-level fields or matrices.
    identity_columns = {'sample_id', 'patient_id', 'participant_id', 'qc_id', 'identity', 'identity_id', 'feature_id'}
    roots = [MANUSCRIPT / 'generated', ROOT / 'results/final/revision_inputs']
    for f in [f for root in roots for pattern in ('*.csv', '*.csv.gz') for f in root.rglob(pattern)]:
        rows = read(f)
        for row in rows:
            clinical = 'clinical_fiams' in row.values() or any('Clinical FIA' in str(v) for v in row.values())
            if clinical:
                require(not any(row.get(k) for k in identity_columns), 'Clinical identifiers in ' + str(f.relative_to(ROOT)))
    # Input preparation is unchanged since the 27 October revision. Table regeneration was
    # repeated for the WiNN 0.1.6 release from the bundled supplement_20261111 inputs.
    for folder, name, count in [('revision_20261027', 'input_preparation_checks.csv', 5),
                                ('revision_20261111', 'table_reproduction_checks.csv', 56)]:
        rows = read(ROOT / 'analysis/config' / folder / name)
        require(len(rows) == count, 'Incomplete recorded validation: ' + name)
        if name == 'table_reproduction_checks.csv':
            require(all(row['passed'].lower() == 'true' for row in rows), 'Table regeneration failed')
        else:
            require(all(row['exact_matrix_match'] == 'TRUE' and row['roles_features_candidates_folds_match'] == 'TRUE' for row in rows), 'Input preparation differs')
    baseline_checks = read(MANUSCRIPT / 'generated/review_validation/baseline_biology_validation_checks.csv')
    # The WiNN 0.1.6 refresh reuses the unchanged baseline fits and bootstrap draws, so it
    # records the 91 checks that apply to recomputed WiNN rows and reused baseline rows.
    require(len(baseline_checks) == 91 and all(row['passed'] == 'TRUE' for row in baseline_checks), 'Baseline/identity validation failed')
    require(not list(MANUSCRIPT.rglob('*.rds')), 'Matrix file unexpectedly included in manuscript bundle')
    for name in ['F1_winn_overview.pdf', 'F2_simulation.pdf', 'F3_replicates.pdf', 'F4_biological_structure.pdf', 'Main_figures.pdf', 'WiNN_figures.pdf']:
        require((MANUSCRIPT / 'figures' / name).is_file(), 'Missing main figure ' + name)
    for name in ['references.bib', 'oup-authoring-template.cls']:
        require((MANUSCRIPT / name).is_file(), 'Missing manuscript source ' + name)
    # The immutable analysis revision precedes the documentation commit that
    # cites it. The later commit supplies both documents together.
    documents = [(MANUSCRIPT / name).is_file() for name in ['main.tex', 'supplementary.tex']]
    require(documents[0] == documents[1], 'Main and supplementary source must be released together')
    require(not (ROOT / 'results/final/figures').exists(), 'Stale previous-release figure tree remains')

validate()
if '--write-manifest' in sys.argv:
    with MANIFEST.open('w', newline='') as f:
        writer = csv.DictWriter(f, ['path', 'size_bytes', 'sha256'])
        writer.writeheader()
        for name in tracked_files():
            data = (ROOT / name).read_bytes()
            writer.writerow(dict(path=name, size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
rows = read(MANIFEST)
require({r['path'] for r in rows} == set(tracked_files()), 'Manifest does not cover the current public file set')
for row in rows:
    data = (ROOT / row['path']).read_bytes()
    require(len(data) == int(row['size_bytes']) and hashlib.sha256(data).hexdigest() == row['sha256'], 'Manifest mismatch: ' + row['path'])
print('PASS: WiNN 0.1.6 results, corrected RFSC results, selection settings, table regeneration, publication assets, privacy checks and', len(rows), 'file hashes')
print('Document source:', 'included' if (MANUSCRIPT / 'main.tex').exists() else 'separate citation/documentation commit; supporting assets included')
