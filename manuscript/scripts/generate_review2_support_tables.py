#!/usr/bin/env python3
"""Create transparent stage and acquisition-design tables from aggregate CSVs."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
V=ROOT/'generated/review_validation';T=ROOT/'generated/supplementary_tables'
def read(n):
 with (V/n).open() as f:return list(csv.DictReader(f))
def table(name,caption,label,cols,head,rows):
 s='\\footnotesize\n\\begin{longtable}{'+cols+'}\n\\caption{'+caption+'}\\label{'+label+'}\\\\\n\\toprule\n'+head+' \\\\\n\\midrule\n\\endfirsthead\n\\toprule\n'+head+' \\\\\n\\midrule\n\\endhead\n\\bottomrule\n\\endlastfoot\n'
 s+='\n'.join(' & '.join(r)+r' \\' for r in rows)+'\n\\end{longtable}\n'
 s=s.replace('\\endfirsthead\n\\toprule', '\\endfirsthead\n\\multicolumn{'+str(head.count('&')+1)+'}{l}{Table~\\thetable\\ (continued)} \\\\\n\\toprule')
 (T/name).write_text(s)
def esc(x):return str(x).replace('_',r'\_').replace('%',r'\%')
d=read('preparation_biology_summary.csv');rows=[]
for sc,lab in [('balanced_sparse','Balanced sparse'),('rare_subgroup','Rare subgroup'),('broad_coordinated','Broad coordinated')]:
 vals=[next(x for x in d if x['scenario']==sc and x['stage']==st) for st in ['mad_only','pqn_only','mad_then_pqn']]
 rows.append([lab]+[f"{float(x['median_recovery_ratio']):.4f}" for x in vals])
table('table_s29_preparation_response_patterns.tex','Stage-specific biological effect recovery in artifact-free simulations. Each cell is the mean over 100 seeds of the median signed corrected-to-clean effect ratio across responsive features. A ratio of one denotes unchanged effects. These tests isolate preparation and thresholded PQN from drift and batch adjustment.','tab:supp_stage_response_patterns','lrrr','Response design & MAD alone & PQN alone & MAD + PQN',rows)
names={'mtbls79':'MTBLS79','sacurine':'Sacurine','waveica':'WaveICA','batchcorr_set1':'BatchCorr','clinical_fiams':'Clinical FIA-MS'}
g=read('prepared_design_group_by_batch.csv');a=read('prepared_design_run_order_balance.csv');rows=[]
for ds in ['mtbls79','sacurine','waveica']:
 batches=list(dict.fromkeys(x['batch'] for x in g if x['dataset']==ds))
 for b in batches:
  counts=[x for x in g if x['dataset']==ds and x['batch']==b];counts.sort(key=lambda x:x['group'])
  countstr='; '.join(esc(x['group'])+': '+x['count'] for x in counts)
  r=next(x for x in a if x['dataset']==ds and x['batch']==b and x['statistic']=='run_order_R2_by_categorical_group')
  rows.append([names[ds],esc(b),countstr,f"{float(r['value']):.4f}"])
table('table_s30_acquisition_group_balance.tex','Biological group counts by analytical batch and linear association with acquisition order. The final column is the ordinary $R^2$ from regressing within-batch run order on the group indicator. C/S denote cow/sheep; Sacurine uses deposited F/M gender labels; WaveICA retains deposited groups 0/1. Quartile-level group counts and continuous-covariate medians and interquartile ranges are supplied as aggregate source CSVs.','tab:supp_acquisition_balance','llp{60mm}r','Dataset & Batch & Study-group counts & Order $R^2$',rows)
r=read('prepared_design_reference_placement.csv');rows=[]
for ds in names:
 for b in dict.fromkeys(x['batch'] for x in r if x['dataset']==ds):
  v=[x for x in r if x['dataset']==ds and x['batch']==b and x['assignment'] in ['training','hidden','external']]
  for x in v:
   role={'training':'Training','hidden':'Evaluation','external':'NIST evaluation'}[x['assignment']]
   gap=x['median_spacing'];gap='--' if gap=='NA' else f'{float(gap):g}'
   rows.append([names[ds],esc(b),role,x['count'],x['first_position']+'--'+x['last_position'],gap])
table('table_s31_reference_acquisition_placement.tex','Reference placement within each analytical batch. Positions count all correction units in acquisition order; clinical units are technical-injection pair means. The spacing column is the median distance between consecutive references of that role, in correction units. A dash denotes a single reference.','tab:supp_reference_placement','lllrrr','Dataset & Batch & Role & Count & First--last position & Median spacing',rows)
print('Generated stage-response, group-balance and reference-placement tables.')
