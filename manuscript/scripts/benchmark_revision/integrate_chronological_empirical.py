#!/usr/bin/env python3
"""Update only current empirical QC-RFSC points/selection from chronology-safe fits."""
import csv
import hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
V=ROOT/'generated/review_validation'
E=V/'chronology_empirical_20261122'
labels={'mtbls79':'MTBLS79','batchcorr_set1':'BatchCorr Set 1','clinical_fiams':'Clinical FIA-MS','sacurine':'Sacurine','waveica':'WaveICA'}
def read(p):
    with open(p) as f:return list(csv.DictReader(f))
def write(p,rows,fields=None):
    with open(p,'w') as f:
        out=csv.DictWriter(f,fieldnames=fields or list(rows[0]));out.writeheader();out.writerows(rows)
points=[]; summary=[]; audit=[]; selections=[]
for ds,label in labels.items():
    fit=E/'results'/ds; evaluation=E/'empirical_evaluation'/ds
    assert 'status=complete' in (fit/'COMPLETE.txt').read_text()
    assert 'status=complete' in (evaluation/'COMPLETE.txt').read_text()
    assertions=read(fit/'REGfit_training_index_assertions.csv')
    assert assertions and all(r['chronological']=='TRUE' and int(r['unauthorized'])==0 and r['used']==r['authorized'] for r in assertions)
    selected=read(fit/'selected_parameters.csv')[0]
    scores=read(fit/'selection_summary.csv')
    best=[r for r in scores if r['selected']=='TRUE']
    assert len(best)==1 and best[0]['candidate_id']==selected['candidate_id']
    diag=read(fit/'revision_diagnostic_metrics.csv')
    raw=next(r for r in diag if r['method']=='raw');corrected=next(r for r in diag if r['method']=='corrected')
    assert int(corrected['nonfinite'])==0
    metrics=read(evaluation/'corrected/metrics.csv')
    values={r['metric']:float(r['value']) for r in metrics if r['value'] not in ('NA','','NaN')}
    assert abs(values['heldout_qc_cv_mean']-float(corrected['heldout_cv_percent']))<1e-10
    points+=metrics
    summary.append(dict(dataset_key=ds,dataset=label,raw_cv_percent=float(raw['heldout_cv_percent']),qcrfsc_cv_percent=values['heldout_qc_cv_mean'],fallback_features=int(corrected['uncorrected_features']),retained_features=int(values['retained_features']),selected_parameters=selected['params'],residual_gam_percent=100*values['residual_run_order_gam_mean_deviance'],batch_weighted_pc_r2=values['batch_weighted_pc_r2_categorical'],replicate_profile_icc=values.get('genuine_replicate_icc_median',''),replicate_profile_pearson=values.get('sample_replicate_pearson_median','')))
    selections.append(dict(dataset_key=ds,candidate_id=selected['candidate_id'],parameters=selected['params'],asserted_calls=len(assertions),unauthorized_training=0,all_chronological=True))
    audit+=read(fit/'chronology_matrix_audit.csv')
write(V/'qcrfsc_completed_empirical_all_endpoints.csv',points)
write(V/'qcrfsc_completed_empirical_summary.csv',summary)
write(ROOT/'generated/supplementary_table_data/table_s26_qcrfsc_empirical.csv',summary)
write(V/'qcrfsc_chronological_selection_summary.csv',selections)
write(V/'qcrfsc_chronology_matrix_audit.csv',audit)
# New table cells must match the ordinary model used by the other table rows,
# not the generic limma count in the separate full point export.
extra=read(V/'qcrfsc_waveica_biology_points.csv')
raw=next(r for r in extra if r['method_id']=='raw')
reference=next(r for r in read(ROOT/'generated/revision_panel_inputs/biological_association_summary.csv')
               if r['dataset_key']=='waveica' and r['method_id']=='raw')
assert int(raw['biological_associated_features'])==int(reference['significant_features'])
assert abs(float(raw['standardized_effect_magnitude'])-float(reference['median_absolute_standardized_effect']))<1e-12
assert abs(float(raw['within_group_variance'])-float(reference['median_within_group_or_residual_variance']))<1e-12
feature_icc=read(V/'qcrfsc_mtbls79_feature_icc.csv')
raw_icc=next(float(r['metabolite_icc_a1_median']) for r in feature_icc if r['method_id']=='raw')
assert abs(raw_icc-0.548997409303439)<1e-12
rfsc_icc=next(float(r['metabolite_icc_a1_median']) for r in feature_icc if r['method_id']=='qc_rfsc')
assert abs(rfsc_icc-next(float(r['value']) for r in points if r['dataset_key']=='mtbls79' and r['metric']=='metabolite_icc_a1_median'))<1e-12
manifest_path=ROOT/'generated/supplementary_table_data/endpoint_free_selection_manifest.csv'
manifest=read(manifest_path)
for row in manifest:
    if row['method_id']!='qc_rfsc':continue
    if row['dataset_key']=='simulation':
        row['selection_criterion']='prespecified before SIM01; chronological injection ranks; no parameter selection'
        row['departure_from_native_default']='TRUE'
        row['departure_justification']='Native RF algorithm with only first/last-column QC guards removed in isolated function copies; no QC relocation. Prespecified ntree=500, coCV=30, Frule=0.8; no grid search.'
        continue
    ds=row['dataset_key'];s=next(x for x in selections if x['dataset_key']==ds)
    row['selected_parameters']=s['parameters']
    key='::'.join((ds,'QC-RFSC',s['candidate_id'],'full_training','post_selection'))
    row['selection_seed']=str(int(hashlib.sha256(key.encode()).hexdigest()[:7],16)+1)
    row['selection_criterion']='training-QC score with chronological injection ranks, masked reference IDs and abundance restoration'
    row['departure_justification']='Native RF algorithm with only first/last-column QC guards removed in isolated function copies; no QC relocation. Four-candidate training-QC selection.'
write(manifest_path,manifest)
path=ROOT/'generated/supplementary_tables/table_s26_qcrfsc_empirical.tex'
lines=path.read_text().splitlines()
for i,line in enumerate(lines):
    for r in summary:
        if line.startswith(r['dataset']+' & ') and 'multicolumn' not in line:
            lines[i]=f"{r['dataset']} & {r['raw_cv_percent']:.3f} & {r['qcrfsc_cv_percent']:.3f} & {r['fallback_features']} & {r['retained_features']} & {r['residual_gam_percent']:.3f} & {r['batch_weighted_pc_r2']:.5f}"+r' \\'
        elif line.startswith(r['dataset']+' & ') and r['dataset_key'] in ('batchcorr_set1','clinical_fiams'):
            lines[i]=r['dataset']+r' & \multicolumn{3}{l}{ICC(A,1): '+f"{r['replicate_profile_icc']:.4f}"+r'} & \multicolumn{3}{l}{Pearson correlation: '+f"{r['replicate_profile_pearson']:.4f}"+r'} \\'
path.write_text('\n'.join(lines)+'\n')
path=ROOT/'generated/supplementary_tables/table_s21_selected_parameters.tex'
lines=path.read_text().splitlines()
for i,line in enumerate(lines):
    for r in summary:
        if line.startswith(r['dataset']+' & QC-RFSC &'):
            lines[i]=r['dataset']+' & QC-RFSC & '+r['selected_parameters']+r' & 4 & Training-QC score \\'
path.write_text('\n'.join(lines)+'\n')
print('PASS: five chronological empirical fits, exact authorized references and current selection integrated; other methods unchanged.')
