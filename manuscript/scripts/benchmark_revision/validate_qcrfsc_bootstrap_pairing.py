#!/usr/bin/env python3
"""Compare raw replay to archived per-draw results, and draw hashes to QC-RFSC."""
import argparse
import csv
import gzip
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
V=ROOT/'generated/review_validation'
parser=argparse.ArgumentParser()
parser.add_argument('--archive-root',required=True,type=Path)
args=parser.parse_args()
def read(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:
        return list(csv.DictReader(f))
checks=[]
def check(name,passed,detail=''):
    checks.append(dict(check=name,passed=passed,detail=detail))
    if not passed: raise ValueError(name+': '+str(detail))
query_archive=read(V/'retrieval_query_weighted_bootstrap_replicates.csv.gz')
for ds in ['mtbls79','clinical_fiams','batchcorr_set1','sacurine','waveica']:
    root=V/'qcrfsc_raw_validation_20261122'/ds
    raw=read(root/'bootstrap_replicates.csv.gz')
    archive=read(args.archive_root/ds/'bootstrap_replicates.csv')
    archive=[r for r in archive if r['method_id']=='raw']
    for metric in sorted(set(r['metric'] for r in raw)):
        a={int(r['bootstrap']):float(r['value']) for r in raw if r['metric']==metric}
        source=query_archive if metric.startswith('cross_batch_') else archive
        b={int(r['bootstrap']):float(r['value']) for r in source if r['method_id']=='raw' and r['metric']==metric and (not metric.startswith('cross_batch_') or r['dataset']==ds)}
        check(ds+' '+metric+' complete archived pairing',set(a)==set(b)==set(range(1,2001)))
        error=max(abs(a[i]-b[i]) for i in a)
        check(ds+' '+metric+' all raw draws reproduce',error<1e-8,error)
    rm=read(root/'input_manifest.csv')
    qm=read(V/'qcrfsc_bootstrap_20261122'/ds/'input_manifest.csv')
    rh=next(r['sha256'] for r in rm if r['object']=='bootstrap_draws')
    qh=next(r['sha256'] for r in qm if r['object']=='bootstrap_draws')
    check(ds+' raw and QC-RFSC draw hashes identical',rh==qh,rh)
with open(V/'qcrfsc_bootstrap_pairing_validation.csv','w') as f:
    out=csv.DictWriter(f,fieldnames=list(checks[0]));out.writeheader();out.writerows(checks)
print(f'PASS: {len(checks)} pairing checks; archived raw draws and identical RFSC draw hashes verified.')
