#!/usr/bin/env python3
"""Aggregate the complete prespecified confounding grid; no scenario selection."""
import csv, hashlib, json, math, statistics, sys
from collections import defaultdict
from pathlib import Path

root, out = map(Path, sys.argv[1:3]); out.mkdir(parents=True, exist_ok=True)
metrics=[]; designs=[]; manifests=[]; sources=[]
for task in range(1,191):
    d=root/f'task_{task:03d}'
    meta=json.loads((d/'complete.json').read_text())
    assert meta['task']==task and meta['package']['version']=='0.1.6'
    assert meta['phenotype_given_to_correction'] is False
    manifests.append(meta)
    for name in ('metrics.csv','design_diagnostics.csv','allocation.csv','feature_effects.csv.gz','complete.json'):
        p=d/name; sources.append(dict(file=f'{d.name}/{name}',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    with (d/'metrics.csv').open() as f: rows=list(csv.DictReader(f))
    assert len(rows)==48
    metrics.extend(rows)
    with (d/'design_diagnostics.csv').open() as f: ds=list(csv.DictReader(f))
    assert len(ds)==1
    ds[0].update(seed_id=meta['seed_id'],axis=meta['axis'],strength=meta['strength']);designs.extend(ds)

by_seed=defaultdict(list)
for m in manifests: by_seed[m['seed_id']].append(m)
assert len(by_seed)==10
for seed,rows in by_seed.items():
    assert len(rows)==19
    for key in ('technical_sha256','baseline_sha256','feature_truth_sha256','protocol_sha256'):
        assert len({r[key] for r in rows})==1,(seed,key)
assert len({m['protocol_sha256'] for m in manifests})==1
for d in designs:
    assert (int(d['n_case']),int(d['n_control']))==(180,270)
    boundary=float(d['strength'])==1 and d['axis'] in ('batch','joint')
    assert (d['design_full_rank']=='FALSE')==boundary
    if boundary: assert int(d['overlap_samples'])==0 and abs(float(d['cramers_v_phenotype_batch'])-1)<1e-12

groups=defaultdict(list)
for r in metrics:
    axes=('batch','run_order','joint') if r['axis']=='none' else (r['axis'],)
    for axis in axes:
        key=(axis,float(r['strength']),r['method'],r['endpoint'],r['metric'])
        groups[key].append(r)
summary=[]
for key,rows in sorted(groups.items()):
    assert len(rows)==10 and len({r['seed_id'] for r in rows})==10
    vals=[float(r['value']) for r in rows if r['value'] not in ('NA','')]
    if key[3]=='marginal': assert len(vals)==10 and all(map(math.isfinite,vals))
    mean=statistics.mean(vals) if vals else None
    sd=statistics.stdev(vals) if len(vals)>1 else None
    summary.append(dict(zip(('axis','strength','method','endpoint','metric'),key))|
                   dict(mean=mean,sd=sd,se=sd/math.sqrt(len(vals)) if sd is not None else None,n_seeds=len(vals)))
for r in metrics:
    if r['axis']=='batch' and float(r['strength'])==1:
        match=[s for s in metrics if s['axis']=='joint' and float(s['strength'])==1 and all(s[k]==r[k] for k in ('seed_id','method','endpoint','metric'))]
        assert len(match)==1 and match[0]['value']==r['value'],r
def write(name,rows):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
write('seed_metrics.csv',metrics);write('summary.csv',summary);write('design_diagnostics.csv',designs);write('input_manifest.csv',sources)
(out/'run_manifests.json').write_text(json.dumps(manifests,indent=2)+'\n')
(out/'validation.json').write_text(json.dumps(dict(tasks=190,seeds=10,conditions_per_seed=19,
    technical_artifacts_fixed_within_seed=True,package_version='0.1.6',complete_batch_separation_verified=True,
    joint_and_batch_full_boundary_identical=True,primary_metric_rows=len([r for r in metrics if r['endpoint']=='marginal']),
    protocol_sha256=manifests[0]['protocol_sha256']),indent=2)+'\n')
print('PASS: 190 tasks; fixed artifacts and implementation; every prespecified condition; complete-confounding boundary')
