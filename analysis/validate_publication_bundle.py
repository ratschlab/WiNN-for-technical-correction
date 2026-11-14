#!/usr/bin/env python3
"""Validate the committed publication data and its complete file manifest."""
from pathlib import Path
import csv
import hashlib
import gzip
import io
import json
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
    package = json.loads((ROOT / 'package/frozen_package.json').read_text())
    require(package['version'] == '0.1.7.9999', 'Current package version differs')
    require(package['drift_calibration'] == 'permutation' and package['n_permutations'] == 999 and package['permutation_seed'] == 104729, 'Calibration defaults differ')
    require(hashlib.sha256((ROOT / package['source_archive']).read_bytes()).hexdigest() == package['source_archive_sha256'], 'Frozen package archive differs')
    evidence = ROOT / 'analysis/config/calibration_revision'
    require(read(evidence / 'validation_helper_body_sha256.csv') == read(evidence / 'production_helper_body_sha256.csv'), 'Validation and production algorithms differ')
    selected = read(evidence / 'primary_selected_parameter_preflight.csv')
    require(len(selected) == 18 and all(r['test'] == 'Ljung-Box' and r['drift_calibration'] == 'permutation' and int(r['n_permutations']) == 999 and int(r['permutation_seed']) == 104729 for r in selected), 'Primary fits used unexpected gate settings')
    for name in ['supplement_20261114', 'derivative_20261114']:
        folder = ROOT / 'results/final/revision_inputs' / name
        manifest = read(folder / 'input_manifest.csv')
        require(len(manifest) >= 13, 'Incomplete aggregate archive: ' + name)
        for row in manifest:
            data = gzip.decompress((folder / (row['file'] + '.gz')).read_bytes())
            require(hashlib.sha256(data).hexdigest() == row['uncompressed_sha256'] and len(data) == int(row['uncompressed_bytes']), 'Archive input differs: ' + row['file'])
            require(len(list(csv.DictReader(io.StringIO(data.decode())))) == int(row['rows']), 'Archive row count differs: ' + row['file'])
    release_checks = read(ROOT / 'results/final/revision_inputs/supplement_20261114/validation/release_validation.csv.gz')
    require(len(release_checks) == 19 and all(r['passed'].lower() == 'true' for r in release_checks), 'Current release validation failed')
    current_tables = read(ROOT / 'analysis/config/calibration_revision/table_reproduction_checks.csv')
    require(len(current_tables) >= 56 and all(r['passed'].lower() == 'true' for r in current_tables), 'Current table reproduction failed')
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
    winn_selection = [r for r in manifest if r['method_id'].startswith('winn_')]
    require(len(winn_selection) == 18, 'Incomplete WiNN selection manifest')
    require(all(r['source_version'] == package['version'] and int(r['candidate_count']) == (1 if r['method_id'] == 'winn_fixed_no_qc' else 36) for r in winn_selection), 'WiNN selection metadata differs')
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
    # Keep the historical input-preparation and 0.1.6 table proofs alongside the
    # current release and table checks above.
    for folder, name, count in [('revision_20261027', 'input_preparation_checks.csv', 5),
                                ('revision_20261111', 'table_reproduction_checks.csv', 56)]:
        rows = read(ROOT / 'analysis/config' / folder / name)
        require(len(rows) == count, 'Incomplete recorded validation: ' + name)
        if name == 'table_reproduction_checks.csv':
            require(all(row['passed'].lower() == 'true' for row in rows), 'Table regeneration failed')
        else:
            require(all(row['exact_matrix_match'] == 'TRUE' and row['roles_features_candidates_folds_match'] == 'TRUE' for row in rows), 'Input preparation differs')
    baseline_checks = read(MANUSCRIPT / 'generated/review_validation/baseline_biology_validation_checks.csv')
    # The calibrated WiNN refresh reuses unchanged baseline fits and bootstrap draws, so it
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
print('PASS: permutation-calibrated WiNN results, corrected RFSC results, selection settings, table regeneration, publication assets, privacy checks and', len(rows), 'file hashes')
print('Document source:', 'included' if (MANUSCRIPT / 'main.tex').exists() else 'separate citation/documentation commit; supporting assets included')
