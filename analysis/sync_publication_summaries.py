#!/usr/bin/env python3
"""Synchronize public summary copies from the current manuscript aggregates.

Only QC-RFSC point cells/summary rows are recalculated; all other method rows
are checked against the pre-existing public copies. No matrix fitting occurs.
"""
import csv, math, shutil, statistics
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
M=ROOT/'manuscript'; V=M/'generated/review_validation'; T=ROOT/'results/final/tables'
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def write(p,rows,fields=None):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
def quantile(values,p):
    v=sorted(values);h=(len(v)-1)*p;j=int(h);return v[j]+(h-j)*(v[min(j+1,len(v)-1)]-v[j])
source=M/'generated/revision_panel_inputs/simulation_seed_method_metrics.csv'
old=read(T/'simulation_seed_method_metrics.csv');current=read(source)
def numeric_map(rows):return {(r['seed_id'],r['method_id'],r['metric']):r['value'] for r in rows if r['method_id']!='qc_rfsc'}
assert numeric_map(old)==numeric_map(current),'Non-RFSC seed metrics changed unexpectedly'
shutil.copy2(source,T/'simulation_seed_method_metrics.csv')
summary=read(T/'simulation_seed_summary.csv');groups=defaultdict(list)
for r in current:
    if r['method_id']=='qc_rfsc' and r['value'] not in ('','NA','NaN') and math.isfinite(float(r['value'])):groups[r['metric']].append(r)
nonrf_summary=[r for r in summary if r['method_id']!='qc_rfsc'];new=[]
for metric,rows in groups.items():
    values=[float(r['value']) for r in rows];assert len(values)==10
    new.append(dict(method_id='qc_rfsc',method='QC-RFSC',metric=metric,direction=rows[0]['direction'],units=rows[0]['units'],n=10,mean=statistics.mean(values),sd=statistics.stdev(values),median=statistics.median(values),q1=quantile(values,.25),q3=quantile(values,.75),minimum=min(values),maximum=max(values)))
write(T/'simulation_seed_summary.csv',nonrf_summary+new,list(summary[0]))
for name in ['qcrfsc_completed_empirical_all_endpoints.csv','qcrfsc_completed_empirical_summary.csv','qcrfsc_completed_simulation_summary.csv']:
    shutil.copy2(V/name,T/name)
settings=M/'generated/supplementary_table_data/endpoint_free_selection_manifest.csv'
for dst in [T/settings.name,ROOT/'analysis/config'/settings.name]:shutil.copy2(settings,dst)
# Root benchmark-wide tables use representative SIM01, not the ten-seed mean.
combined_path=T/'combined_primary_benchmark.csv';combined=read(combined_path)
dataset_labels={r['dataset_key']:r['dataset'] for r in combined}
replacements=read(V/'qcrfsc_completed_empirical_all_endpoints.csv')+[r for r in current if r['method_id']=='qc_rfsc' and r['seed_id']=='SIM01']
keys=[(r['dataset_key'],r['metric']) for r in replacements];assert len(keys)==len(set(keys))
new_combined=[]
for r in replacements:
    row={k:r.get(k,'') for k in combined[0]};row.update(dataset=dataset_labels[r['dataset_key']],method_id='qc_rfsc',method='QC-RFSC');new_combined.append(row)
assert {r['dataset_key'] for r in new_combined}=={'simulation','mtbls79','clinical_fiams','batchcorr_set1','sacurine','waveica'}
nonrf=[r for r in combined if r['method_id']!='qc_rfsc'];write(combined_path,nonrf+new_combined,list(combined[0]))
assert [r for r in read(combined_path) if r['method_id']!='qc_rfsc']==nonrf
wide_path=T/'primary_method_metrics_wide.csv';wide=read(wide_path);prior_nonrf=[r.copy() for r in wide if r['method_id']!='qc_rfsc']
for row in wide:
    if row['method_id']!='qc_rfsc':continue
    values={r['metric']:r['value'] for r in new_combined if r['dataset_key']==row['dataset_key']}
    for k in row:
        if k not in ('dataset_key','dataset','method_id','method'):row[k]=values.get(k,'')
write(wide_path,wide)
assert [r for r in read(wide_path) if r['method_id']!='qc_rfsc']==prior_nonrf
print('PASS: current simulation, empirical RFSC and selection copies synchronized; other methods unchanged.')
