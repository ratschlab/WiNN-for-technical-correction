#!/usr/bin/env python3
"""Refresh compact WiNN result rows while preserving published comparators exactly."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--final-root', type=Path, required=True)
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
target_root = root / 'results/final/tables'
specs = {
    'combined_primary_benchmark.csv': ['dataset_key','method_id','metric'],
    'endpoint_free_selection_manifest.csv': ['dataset_key','method_id'],
    'primary_method_metrics_wide.csv': ['dataset_key','method_id'],
    'simulation_seed_method_metrics.csv': ['seed_id','method_id','metric'],
    'simulation_seed_summary.csv': ['method_id','metric'],
}
def read(path):
    with path.open(newline='') as f:
        r = csv.DictReader(f); return r.fieldnames, list(r)
def winn(row):
    return row['method_id'].startswith('winn_')
def digest(rows):
    return hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest()
pending = []
proof = []
for name, fields in specs.items():
    columns, old = read(target_root / name)
    fresh_columns, fresh = read(a.final_root / 'tables' / name)
    assert set(columns) == set(fresh_columns), name + ': schema differs'
    key = lambda r: tuple(r[k] for k in fields)
    old_winn = [r for r in old if winn(r)]
    new_winn = [r for r in fresh if winn(r)]
    assert len({key(r) for r in old_winn}) == len(old_winn), name
    assert len({key(r) for r in new_winn}) == len(new_winn), name
    lookup = {key(r): r for r in new_winn}
    assert set(lookup) == {key(r) for r in old_winn}, name + ': WiNN keys differ'
    updated = [lookup[key(r)] if winn(r) else r for r in old]
    before = [r for r in old if not winn(r)]
    after = [r for r in updated if not winn(r)]
    assert before == after, name + ': comparator rows changed'
    pending.append((target_root / name,columns,updated))
    proof.append(dict(file=name,refreshed_winn_rows=len(new_winn),unchanged_other_rows=len(before),unchanged_rows_sha256=digest(before)))
for path, columns, rows in pending:
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(rows)
(root / 'analysis/config/calibration_revision/public_result_table_staging.json').write_text(json.dumps(proof,indent=2)+'\n')
print('Updated',len(pending),'compact result tables; all comparator dictionaries retained exactly.')
