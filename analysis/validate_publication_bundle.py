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
    # Current RFSC execution preserves acquisition ranks. Historical adapter
    # checks above are retained as provenance, not substituted for these checks.
    review = MANUSCRIPT / 'generated/review_validation'
    current_checks = read(review / 'qcrfsc_chronology_simulation/validation_checks.csv')
    require(len(current_checks) == 220 and all(r['passed'].lower() == 'true' for r in current_checks),
            'Current chronology-preserving simulation validation failed')
    for seed in [f'SIM{n:02d}' for n in range(1, 11)]:
        fit = review / 'qcrfsc_chronology_simulation/fits' / seed
        assertions = read(fit / 'REGfit_training_index_assertions.csv')
        require(bool(assertions) and all(r['chronological'] == 'TRUE' and r['unauthorized'] == '0'
                and r['used'] == r['authorized'] for r in assertions), 'RFSC simulation training/order audit failed')
        audit = read(review / 'qcrfsc_chronology_simulation/evaluation' / seed / 'ljung_box_validity_audit.csv')
        require(len(audit) == 1 and audit[0]['profiles'] == audit[0]['finite_adjusted_tests'] == '5000'
                and audit[0]['unavailable_tests'] == '0', 'Simulation test denominator audit differs')
    chronological = read(review / 'qcrfsc_chronology_matrix_audit.csv')
    datasets = {'mtbls79', 'clinical_fiams', 'batchcorr_set1', 'sacurine', 'waveica'}
    require(len(chronological) == 5 and {r['dataset_key'] for r in chronological} == datasets,
            'Current empirical chronology audit incomplete')
    hashes = read(review / 'qcrfsc_bootstrap_input_manifest.csv')
    for row in chronological:
        require(row['chronological_ids_preserved'] == 'TRUE', 'Empirical acquisition order differs')
        matched = [r for r in hashes if r['dataset_key'] == row['dataset_key'] and r['object'] == 'corrected_matrix']
        require(len(matched) == 1 and matched[0]['sha256'] == row['new_matrix_file_sha256'],
                'RFSC uncertainty uses a different matrix from current corrected endpoints')
        a = read(review / 'chronology_empirical_20261122/results' / row['dataset_key'] / 'REGfit_training_index_assertions.csv')
        require(bool(a) and all(r['chronological'] == 'TRUE' and r['unauthorized'] == '0'
                and r['used'] == r['authorized'] for r in a), 'Empirical RFSC native training/order audit failed')
    preprocessing = [r for r in read(review / 'qcrfsc_preprocessing_input_audit.csv') if r['injection_subset'] == 'all']
    require(len(preprocessing) == 15, 'RFSC preprocessing input audit incomplete')
    for row in preprocessing:
        dataset = row['dataset']
        if dataset in datasets:
            matched = [r for r in hashes if r['dataset_key'] == dataset and r['object'] == 'prepared_input']
            require(len(matched) == 1 and matched[0]['sha256'] == row['input_file_sha256'],
                    'RFSC preprocessing audit uses a different empirical input')
        else:
            fit_hashes = read(review / 'qcrfsc_chronology_simulation/fits' / dataset / 'input_identity_checks.csv')
            raw = [r for r in fit_hashes if r['object'] == 'raw_intensity']
            require(len(raw) == 1 and raw[0]['observed_sha256'] == row['matrix_sha256'],
                    'RFSC preprocessing audit uses a different simulation input')
        require(int(row['n_na']) == int(row['n_nonfinite']) == int(row['n_negative']) == 0,
                'Unexpected unavailable or negative prepared RFSC input')
    wave = next(r for r in preprocessing if r['dataset'] == 'waveica')
    require(int(wave['n_zero']) == 48820 and int(wave['n_missing_after_FilterMV']) == 47912 and
            int(wave['n_features']) - int(wave['n_features_after_FilterMV']) == 9,
            'WaveICA native preprocessing audit differs from documented counts')
    boot = read(review / 'qcrfsc_bootstrap_estimates.csv')
    require(len(boot) == 18 and {r['dataset_key'] for r in boot} == datasets and
            all(int(r['valid_bootstrap']) == 2000 and math.isfinite(float(r['bootstrap_se'])) for r in boot),
            'QC-RFSC biological uncertainty incomplete')
    for name, minimum in [('qcrfsc_bootstrap_validation.csv', 82), ('qcrfsc_bootstrap_pairing_validation.csv', 41)]:
        rows = read(review / name)
        require(len(rows) >= minimum and all(r['passed'].lower() == 'true' for r in rows),
                'Current paired bootstrap checks failed: ' + name)
    endpoints = read(ROOT / 'results/final/tables/qcrfsc_completed_empirical_all_endpoints.csv')
    require(len(endpoints) == 84 and len({r['dataset_key'] for r in endpoints}) == 5, 'Incomplete empirical endpoint export')
    require(endpoints == read(review / 'qcrfsc_completed_empirical_all_endpoints.csv'),
            'Root and manuscript empirical QC-RFSC endpoints differ')
    representative = [r for r in source if r['method_id'] == 'qc_rfsc' and r['seed_id'] == 'SIM01']
    rf_expected = {(r['dataset_key'], r['metric']): r['value'] for r in endpoints + representative}
    combined = read(ROOT / 'results/final/tables/combined_primary_benchmark.csv')
    rf_combined = {(r['dataset_key'], r['metric']): r['value'] for r in combined if r['method_id'] == 'qc_rfsc'}
    require(rf_combined == rf_expected, 'Combined public RFSC points differ from current matrices')
    for row in read(ROOT / 'results/final/tables/primary_method_metrics_wide.csv'):
        if row['method_id'] == 'qc_rfsc':
            for (dataset, metric), value in rf_expected.items():
                if dataset == row['dataset_key'] and metric in row:
                    require(row[metric] == value, 'Wide public RFSC points differ')
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
            clinical = 'clinical_fiams' in f.parts or 'clinical_fiams' in row.values() or any('Clinical FIA' in str(v) for v in row.values())
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
    # S7 must use the same frozen 0.1.6 fold results as the manuscript text.
    fold_path = MANUSCRIPT / 'generated/supplementary_table_data/table_s12_lobo_fold_level.csv'
    folded = read(fold_path)
    archive = read(ROOT / 'results/final/revision_inputs/supplement_20261111/derivative/tables/leave_one_batch_out_bootstrap_summary.csv.gz')
    keys = ('dataset_key', 'endpoint', 'method_id', 'metric', 'heldout_batch')
    frozen = {tuple(r[k] for k in keys): r for r in archive}
    require(len(folded) == 408, 'Incomplete current S7 fold table')
    for row in folded:
        original = frozen[tuple(row[k] for k in keys)]
        for field in ('bootstrap_mean', 'ci_lower', 'ci_upper', 'valid_bootstrap'):
            require(math.isclose(float(row[field]), float(original[field]), rel_tol=0, abs_tol=1e-12),
                    'S7 differs from frozen 0.1.6 results: ' + str(tuple(row[k] for k in keys)))
    subprocess.run([sys.executable, str(MANUSCRIPT / 'scripts/audit_supplement_finalization.py')], check=True)
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
print('PASS: WiNN 0.1.6 results, chronology-preserving RFSC results and paired uncertainty, selection settings, table regeneration, publication assets, privacy checks and', len(rows), 'file hashes')
print('Document source:', 'included' if (MANUSCRIPT / 'main.tex').exists() else 'separate citation/documentation commit; supporting assets included')
