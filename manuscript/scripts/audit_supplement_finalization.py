#!/usr/bin/env python3
"""Check that Figure S7 and its text use the frozen WiNN 0.1.6 fold results."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def read(path):
    return list(csv.DictReader(path.open()))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
canonical = ROOT / 'generated/supplementary_table_data/table_s12_lobo_fold_level.csv'
rows = read(canonical)
keys = ('dataset_key', 'endpoint', 'method_id', 'metric', 'heldout_batch')
expected = {tuple(r[k] for k in keys): r for r in rows}
assert len(rows) == len(expected) == 408
assert canonical.read_bytes() == (canonical.parent / 'leave_one_batch_out_bootstrap_summary.csv').read_bytes()
metric_names = {'Pearson effect\ncorrelation': 'pearson', 'Spearman effect\ncorrelation': 'spearman', 'Sign\nconcordance': 'sign_concordance'}
seen = set()
for panel in sorted((ROOT / 'figures/source_data').glob('S7?_lobo_*.csv')):
    for row in read(panel):
        key = tuple(metric_names[row[k]] if k == 'metric' else row[k] for k in keys)
        assert key in expected and key not in seen, (panel, key)
        seen.add(key)
        for field in ('bootstrap_mean', 'ci_lower', 'ci_upper', 'valid_bootstrap'):
            assert math.isclose(float(row[field]), float(expected[key][field]), rel_tol=0, abs_tol=1e-12), (panel, key, field)
assert seen == set(expected)
provenance = json.loads((ROOT / 'generated/review_validation/supplement_finalization_provenance.json').read_text())
assert sha(canonical) == provenance['canonical_fold_source_sha256']
for name, digest in provenance['unchanged_main_figures'].items():
    assert sha(ROOT / name) == digest, name
for row in provenance['updated_s7_files']:
    assert sha(ROOT / row['path']) == row['sha256'], row['path']
print('PASS: all 408 S7 fold rows and intervals match WiNN 0.1.6; all main figures unchanged.')
