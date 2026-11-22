#!/usr/bin/env python3
"""Replace only QC-RFSC simulation cells from the acquisition-order-preserving rerun."""
import csv, math, statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
V=ROOT/'generated/review_validation'
E=V/'qcrfsc_chronology_simulation'
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def write(p,rows,fields=None):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
checks=[]
def check(name,condition,detail=''):
    checks.append(dict(check=name,passed=bool(condition),detail=detail))
    if not condition: raise ValueError(name+': '+str(detail))
source=ROOT/'generated/revision_panel_inputs/simulation_seed_method_metrics.csv'
old=read(source)
archive=E/'superseded_reordered_adapter_metrics.csv'
if not archive.exists():write(archive,[r for r in old if r['method_id']=='qc_rfsc'])
nonrf=[r for r in old if r['method_id']!='qc_rfsc']; corrected=[]
truth_metrics=['truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean']
for n in range(1,11):
    seed=f'SIM{n:02d}'; fit=E/'fits'/seed; ev=E/'evaluation'/seed
    check(seed+' completed','status=complete' in (fit/'COMPLETE.txt').read_text() and 'status=complete' in (ev/'COMPLETE.txt').read_text())
    check(seed+' chronology asserted','acquisition_chronology_preserved=TRUE' in (fit/'COMPLETE.txt').read_text())
    ids=read(fit/'input_identity_checks.csv')
    check(seed+' exact frozen raw/truth',len(ids)==2 and all(r['passed']=='TRUE' and r['observed_sha256']==r['expected_sha256'] for r in ids))
    a=read(fit/'REGfit_training_index_assertions.csv')
    check(seed+' training labels isolated',bool(a) and all(int(r['unauthorized'])==0 and r['used']==r['authorized'] for r in a))
    validity=read(ev/'ljung_box_validity_audit.csv')
    check(seed+' all 5000 Ljung-Box tests finite',len(validity)==1 and int(validity[0]['profiles'])==5000 and int(validity[0]['unavailable_tests'])==0 and int(validity[0]['finite_adjusted_tests'])==5000)
    metrics=read(ev/'metrics.csv'); check(seed+' metric names unique',len({r['metric'] for r in metrics})==len(metrics))
    required=truth_metrics+['heldout_qc_cv_mean','residual_run_order_gam_mean_deviance','batch_weighted_pc_r2_categorical']
    check(seed+' required endpoints finite',all(len([r for r in metrics if r['metric']==m])==1 and math.isfinite(float(next(r['value'] for r in metrics if r['metric']==m))) for m in required))
    check(seed+' native RF input chronology',all(r.get('chronological')=='TRUE' for r in a))
    simple=read(fit/'simulation_truth_metrics.csv')
    for method in ['raw','corrected']:
        sr=[r for r in simple if r['method']==method];check(seed+' '+method+' single row',len(sr)==1)
        for metric in truth_metrics+['heldout_qc_cv_mean']:
            expected=[float(r['value']) for r in (old if method=='raw' else metrics) if r['metric']==metric and (method!='raw' or (r['seed_id']==seed and r['method_id']=='raw'))]
            check(seed+' '+method+' independent '+metric,len(expected)==1 and abs(float(sr[0][metric])-expected[0])<1e-11)
    for r in metrics:
        row={k:r.get(k,'') for k in old[0]};row.update(method_id='qc_rfsc',method='QC-RFSC',dataset_key='simulation',seed_id=seed);corrected.append(row)
check('exact ten seed sets',set(r['seed_id'] for r in corrected)=={f'SIM{n:02d}' for n in range(1,11)})
write(source,nonrf+corrected)
check('all non-RFSC simulation results unchanged',[r for r in read(source) if r['method_id']!='qc_rfsc']==nonrf)
labels={'heldout_qc_cv_mean':'QC CV','residual_run_order_gam_mean_deviance':'GAM','batch_weighted_pc_r2_categorical':'Batch R2',truth_metrics[0]:'RMSE',truth_metrics[1]:'ICC',truth_metrics[2]:'Pearson'}
tablepath=ROOT/'generated/supplementary_table_data/table_s4_simulation.csv';table=read(tablepath)
rf=[r for r in table if r['method_id']=='qc_rfsc'];check('one table RFSC row',len(rf)==1);summary=[]
for metric,label in labels.items():
    values=[float(r['value']) for r in corrected if r['metric']==metric];check(metric+' ten values',len(values)==10)
    mean=statistics.mean(values);se=statistics.stdev(values)/math.sqrt(10);digits=2 if label=='QC CV' else 4
    rf[0][label]=f'{mean:.{digits}f} $\\pm$ {se:.{digits}f}';summary.append(dict(metric=metric,n=10,mean=mean,se=se))
write(tablepath,table);write(V/'qcrfsc_completed_simulation_summary.csv',summary);write(E/'summary.csv',summary)
p=ROOT/'generated/supplementary_tables/table_s4_simulation.tex';lines=p.read_text().splitlines()
check('one TeX RFSC row',sum(l.startswith('QC-RFSC &') for l in lines)==1)
line=' & '.join(['QC-RFSC']+[rf[0][x] for x in ['QC CV','GAM','Batch R2','RMSE','ICC','Pearson']])+r' \\'
p.write_text('\n'.join(line if l.startswith('QC-RFSC &') else l for l in lines)+'\n')
ranking=[]
for r in corrected:
    if r['metric'] not in truth_metrics:continue
    w=[a for a in old if a['method_id']=='winn_fixed_no_qc' and a['seed_id']==r['seed_id'] and a['metric']==r['metric']];check('WiNN reference unique '+r['seed_id']+r['metric'],len(w)==1)
    wr,qr=float(w[0]['value']),float(r['value']);ranking.append(dict(seed_id=r['seed_id'],metric=r['metric'],winn=wr,qcrfsc=qr,winn_better=(wr<qr if 'rmse' in r['metric'] else wr>qr)))
write(E/'winn_comparison.csv',ranking);write(E/'validation_checks.csv',checks)
print('Validated',len(checks),'checks; updated only QC-RFSC simulation outputs.')
for s in summary: print(s)
print('WiNN better on',sum(r['winn_better'] for r in ranking),'of',len(ranking),'paired truth comparisons; this is reported, not enforced.')
