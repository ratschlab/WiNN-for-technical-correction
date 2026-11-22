#!/usr/bin/env python3
"""Integrate aggregate-only corrected-QC-RFSC uncertainty; preserve other methods."""
import csv
import gzip
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V = ROOT / 'generated/review_validation'
PANELS = ROOT / 'generated/revision_panel_inputs'
SOURCE = V / 'qcrfsc_bootstrap_20261122'
DATASETS = ['mtbls79', 'clinical_fiams', 'batchcorr_set1', 'sacurine', 'waveica']
LABEL = dict(zip(DATASETS, ['MTBLS79', 'Clinical FIA-MS', 'BatchCorr Set 1', 'Sacurine', 'WaveICA']))

def read(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as stream:
        return list(csv.DictReader(stream))

def write(path, rows, fields=None):
    fields = fields or list(rows[0])
    with (gzip.open(path, 'wt') if str(path).endswith('.gz') else open(path, 'w')) as stream:
        out = csv.DictWriter(stream, fieldnames=fields)
        out.writeheader()
        out.writerows([{k: row.get(k, '') for k in fields} for row in rows])

def quantile(values, p):
    values = sorted(values); index = (len(values)-1)*p; lo = int(index)
    return values[lo] + (index-lo)*(values[min(lo+1,len(values)-1)]-values[lo])

def append_current(path, rows, method_field='method_id'):
    old = read(path); kept = [r for r in old if r[method_field] != 'qc_rfsc']
    write(path, kept+rows, list(old[0]))

checks = []
def check(name, condition, detail=''):
    checks.append(dict(check=name, passed=bool(condition), detail=detail))
    if not condition:
        raise ValueError(f'{name}: {detail}')

estimates, replicates, manifests, protocols = [], [], [], []
current_matrices = {r['dataset_key']:r['new_matrix_file_sha256']
                    for r in read(V/'qcrfsc_chronology_matrix_audit.csv')}
for ds in DATASETS:
    check(ds+' complete', (SOURCE/ds/'COMPLETE.txt').exists())
    est = read(SOURCE/ds/'bootstrap_estimates.csv')
    rep = read(SOURCE/ds/'bootstrap_replicates.csv.gz')
    for row in est:
        values = [float(r['value']) for r in rep if r['metric']==row['metric']]
        check(ds+' '+row['metric']+' complete draws', len(values)==2000 and all(math.isfinite(x) for x in values))
        check(ds+' '+row['metric']+' SE', abs(statistics.stdev(values)-float(row['bootstrap_se']))<1e-12)
        check(ds+' '+row['metric']+' interval', abs(quantile(values,.025)-float(row['ci_lower']))<1e-12 and abs(quantile(values,.975)-float(row['ci_upper']))<1e-12)
    manifest=read(SOURCE/ds/'input_manifest.csv')
    matrix_hash=next(r['sha256'] for r in manifest if r['object']=='corrected_matrix')
    check(ds+' bootstrap uses chronological matrix',matrix_hash==current_matrices[ds],matrix_hash)
    estimates += est; replicates += rep; manifests += manifest
    if ds in ['batchcorr_set1','clinical_fiams']:
        protocols += read(SOURCE/ds/'matching_protocol.csv')
write(V/'qcrfsc_bootstrap_estimates.csv',estimates)
write(V/'qcrfsc_bootstrap_replicates.csv.gz',replicates)
write(V/'qcrfsc_bootstrap_input_manifest.csv',manifests)
write(V/'qcrfsc_matching_protocol.csv',protocols)
lookup = {(r['dataset_key'],r['metric']):r for r in estimates}
oldpoints = read(V/'qcrfsc_completed_empirical_all_endpoints.csv')
technical = []
for ds in DATASETS:
    for mean, sd in [('heldout_qc_cv_mean','heldout_qc_cv_sd'),
                     ('residual_run_order_gam_mean_deviance','residual_run_order_gam_sd_deviance')]:
        a=next(r for r in oldpoints if r['dataset_key']==ds and r['metric']==mean)
        b=next(r for r in oldpoints if r['dataset_key']==ds and r['metric']==sd)
        check(ds+' '+mean+' descriptive denominator',a['denominator']==b['denominator'])
        technical.append(dict(dataset_key=ds,method_id='qc_rfsc',metric=mean,point_estimate=a['value'],
            descriptive_sem=float(b['value'])/math.sqrt(int(b['denominator'])),denominator=a['denominator'],units=a['units'],
            uncertainty_scope='descriptive feature/profile SEM; features are not independent experimental replicates'))
write(V/'qcrfsc_descriptive_sem.csv',technical)
for row in oldpoints:
    key = row['dataset_key'],row['metric']
    if key not in lookup or key==('mtbls79','biological_weighted_pc_r2'):
        continue
    # Historic generic Pearson uses all pairs; manuscript Pearson uses cross-batch pairs.
    if row['metric']=='sample_replicate_pearson_median':
        continue
    error = abs(float(row['value'])-float(lookup[key]['point_estimate']))
    check('Existing point '+str(key),error<1e-8,error)
mt_rows = read(V/'qcrfsc_mtbls79_sample_means.csv')
for r in mt_rows:
    r['uncertainty']='paired biological-unit intervals supplied in qcrfsc_bootstrap_estimates.csv'
write(V/'qcrfsc_mtbls79_sample_means.csv',mt_rows)
mt = next(r for r in mt_rows if r['method_id']=='qc_rfsc')
check('MTBLS79 sample-mean biology point',abs(float(mt['biological_pc_variance_percent'])/100-float(lookup['mtbls79','biological_weighted_pc_r2']['point_estimate']))<1e-8)

# Add the newly evaluated comparator to panels with identical conditional endpoints.
for path in PANELS.glob('S8*.csv'):
    old=read(path); ds=old[0]['dataset_key']; metric=old[0]['metric']
    row=dict(lookup[ds,metric],dataset_key=ds,dataset=LABEL[ds],method_label='QC-RFSC')
    append_current(path,[row])
for letter,metric in [('A','age_weighted_pc_r2'),('B','bmi_weighted_pc_r2'),('C','gender_weighted_pc_r2')]:
    bx=lookup['sacurine','batch_weighted_pc_r2_categorical']; by=lookup['sacurine',metric]
    row=dict(method_id='qc_rfsc',method='QC-RFSC',method_label='QC-RFSC',
             x=100*float(bx['point_estimate']),x_se=100*float(bx['bootstrap_se']),
             y=100*float(by['point_estimate']),y_se=100*float(by['bootstrap_se']))
    append_current(PANELS/f'S10{letter}_sacurine_{metric}.csv',[row])

# Match query-weighted retrieval intervals to the current WiNN paired draws.
old_rep=read(V/'retrieval_query_weighted_bootstrap_replicates.csv.gz')
old_est=read(V/'frozen_identity_retrieval_estimates.csv')
new_est=[]; new_rep=[]; new_delta=[]
for ds in ['batchcorr_set1','clinical_fiams']:
    for metric in ['cross_batch_top1_accuracy','cross_batch_mean_reciprocal_rank']:
        estimate=lookup[ds,metric]
        new_est.append(dict(estimate,dataset=ds,uncertainty_scope='paired 2000-draw identity-cluster bootstrap retaining query-count weights; fitted correction and candidate set fixed'))
        rfsc={int(r['bootstrap']):float(r['value']) for r in replicates if r['dataset_key']==ds and r['metric']==metric}
        winn={int(r['bootstrap']):float(r['value']) for r in old_rep if r['dataset']==ds and r['metric']==metric and r['method_id']=='winn_fixed_no_qc'}
        check(ds+' '+metric+' paired draw index',set(rfsc)==set(winn)==set(range(1,2001)))
        delta=[winn[i]-rfsc[i] for i in range(1,2001)]
        wp=next(float(r['point_estimate']) for r in old_est if r['dataset']==ds and r['method_id']=='winn_fixed_no_qc' and r['metric']==metric)
        new_delta.append(dict(metric=metric,direction='higher',winn_method_id='winn_fixed_no_qc',winn_method='WiNN',comparator_method_id='qc_rfsc',comparator_method='QC-RFSC',point_direction_normalized_difference=wp-float(estimate['point_estimate']),bootstrap_mean_difference=statistics.mean(delta),ci_lower=quantile(delta,.025),ci_upper=quantile(delta,.975),probability_winn_advantage=statistics.mean(x>1e-12 for x in delta),valid_bootstrap=2000,difference_definition='WiNN minus comparator; positive favors WiNN',dataset=ds,uncertainty_scope='paired identity resampling with query-count weights; fitted correction and candidate set fixed'))
        new_rep += [dict(dataset=ds,bootstrap=i,method_id='qc_rfsc',metric=metric,value=rfsc[i]) for i in range(1,2001)]
append_current(V/'frozen_identity_retrieval_estimates.csv',new_est)
append_current(V/'frozen_identity_retrieval_paired_differences.csv',new_delta,'comparator_method_id')
append_current(V/'retrieval_query_weighted_bootstrap_replicates.csv.gz',new_rep)

# Preserve the table's wording and every other row; insert the two new rows.
path=ROOT/'generated/supplementary_tables/table_s23_identity_retrieval.tex'
lines=[line for line in path.read_text().splitlines() if ' & QC-RFSC & ' not in line]
def fmt(r,key='point_estimate'):
    return f"{float(r[key]):.3f} [{float(r['ci_lower']):.3f}, {float(r['ci_upper']):.3f}]"
for ds in ['batchcorr_set1','clinical_fiams']:
    a=lookup[ds,'cross_batch_top1_accuracy']; b=lookup[ds,'cross_batch_mean_reciprocal_rank']
    d=next(r for r in new_delta if r['dataset']==ds and r['metric']=='cross_batch_top1_accuracy')
    where=next(i for i,line in enumerate(lines) if line.startswith(LABEL[ds]+' & QC-RLSC & '))+1
    lines.insert(where,' & '.join([LABEL[ds],'QC-RFSC',fmt(a),fmt(b),fmt(d,'point_direction_normalized_difference')])+r' \\')
path.write_text('\n'.join(lines)+'\n')
write(V/'qcrfsc_bootstrap_validation.csv',checks)
print(f'PASS: {len(checks)} checks; {len(estimates)} endpoints from five corrected matrices; 2000 valid draws each.')
