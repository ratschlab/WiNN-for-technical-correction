#!/usr/bin/env python3
"""Archive aggregate/public derivative inputs; never copy clinical identities or matrices."""
import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--final-root', required=True, type=Path)
p.add_argument('--destination-root', type=Path, default=Path(__file__).resolve().parents[1] / 'results/final/revision_inputs')
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
for name, source in [('supplement', a.final_root), ('derivative', a.final_root / 'derivative')]:
    previous = root / 'results/final/revision_inputs' / (name + '_20261111')
    destination = a.destination_root / (name + '_20261114')
    files = [f.relative_to(previous).with_suffix('') for f in previous.rglob('*.csv.gz')]
    if name == 'derivative':
        files += [Path(ds) / 'biological_point_metrics.csv' for ds in ['batchcorr_set1', 'clinical_fiams', 'mtbls79', 'sacurine', 'waveica']]
    records = []
    for rel in sorted(set(files)):
        raw = (source / rel).read_bytes()
        rows = list(csv.DictReader(io.StringIO(raw.decode())))
        for row in rows:
            clinical = 'clinical_fiams' in row.values() or 'clinical_fiams' in rel.parts or any('Clinical FIA' in str(v) for v in row.values())
            if clinical:
                forbidden = {'sample_id','patient_id','participant_id','qc_id','identity','identity_id','feature_id'}
                assert not any(row.get(k) for k in forbidden), str(rel)
        target = destination / (str(rel) + '.gz')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(gzip.compress(raw, mtime=0))
        records.append(dict(file=str(rel), rows=len(rows), uncompressed_sha256=hashlib.sha256(raw).hexdigest(), uncompressed_bytes=len(raw), stored_bytes=target.stat().st_size))
    with (destination / 'input_manifest.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    print(name, len(records), 'files archived')
