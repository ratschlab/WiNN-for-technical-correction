#!/usr/bin/env python3
"""Check row-order independence and refusal of changed comparator draws."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

script = Path(__file__).with_name('restore_derivative_comparator_rows.py')
fields = ['method_id', 'metric', 'bootstrap', 'value']
previous = [dict(method_id=m, metric='endpoint', bootstrap=str(b), value=str(v))
            for b in (1, 2) for m, v in [('raw', b / 10), ('combat', b / 20), ('winn_fixed_no_qc', b / 30)]]

def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fields)
        writer.writeheader()
        writer.writerows(rows)

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    new, old = root / 'new', root / 'old'
    new.mkdir(); old.mkdir()
    write(old / 'bootstrap_replicates.csv', previous)
    changed = [dict(row) for row in reversed(previous)]
    for row in changed:
        if row['method_id'].startswith('winn_'):
            row['value'] = '0.777'
    write(new / 'bootstrap_replicates.csv', changed)
    call = [sys.executable, str(script), '--new', str(new), '--previous', str(old)]
    subprocess.run(call, check=True, capture_output=True, text=True)
    with (new / 'bootstrap_replicates.csv').open(newline='') as f:
        restored = list(csv.DictReader(f))
    assert restored == changed, 'Draw identifiers or new WiNN values changed'
    restored[0]['value'] = '0.888'
    raw = next(row for row in restored if row['method_id'] == 'raw')
    raw['value'] = '0.999'
    write(new / 'bootstrap_replicates.csv', restored)
    rejected = subprocess.run(call, capture_output=True, text=True)
    assert rejected.returncode != 0 and 'comparator mismatch' in rejected.stderr
print('PASS: shuffled draws preserved; altered comparator draws rejected.')
