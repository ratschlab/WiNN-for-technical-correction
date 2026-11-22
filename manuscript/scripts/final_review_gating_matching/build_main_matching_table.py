#!/usr/bin/env python3
"""Compact main-text table; all values are percentages from paired matching outputs."""
import csv
from pathlib import Path
root=Path(__file__).resolve().parents[2]
rows=list(csv.DictReader((root/'generated/final_review_gating_matching/main_matching_table.csv').open()))
rf=root/'generated/review_validation/qcrfsc_bootstrap_estimates.csv'
assert rf.exists(), "Completed QC-RFSC bootstrap is required"
estimates=[r for r in csv.DictReader(rf.open()) if r['metric']=='cross_batch_top1_accuracy' and r['method_id']=='qc_rfsc']
assert len(estimates)==2 and {r['dataset_key'] for r in estimates}=={'batchcorr_set1','clinical_fiams'}
v={r['dataset_key']:float(r['point_estimate'])*100 for r in estimates}
rows.insert(3,dict(method='QC-RFSC',BatchCorr_accuracy_percent=v['batchcorr_set1'],clinical_accuracy_percent=v['clinical_fiams']))
assert len(rows)==7 and len({r['method'] for r in rows})==7
baseline_estimates=list(csv.DictReader((root/'generated/review_validation/frozen_identity_retrieval_estimates.csv').open()))
for method_id,label in [('mean_batch_all','Mean-only batch'),('wehrens_ancova_all','All-injection ANCOVA')]:
    b=[r for r in baseline_estimates if r['method_id']==method_id and r['metric']=='cross_batch_top1_accuracy']
    assert len(b)==2 and {r['dataset'] for r in b}=={'batchcorr_set1','clinical_fiams'}
    v={r['dataset']:100*float(r['point_estimate']) for r in b}
    rows.append(dict(method_id=method_id,method=label,BatchCorr_accuracy_percent=v['batchcorr_set1'],clinical_accuracy_percent=v['clinical_fiams']))
assert len(rows)==9 and len({r['method'] for r in rows})==9
lines=[r'\begin{table}[t]',r'\centering\small',r'\caption{Cross-batch sample matching (top-1 accuracy, \%).}',r'\label{tab:matching}',r'\begin{tabular}{lrr}',r'\toprule',r'Method & BatchCorr & Clinical \\',r'\midrule']
for r in rows:lines.append(f"{r['method']} & {float(r['BatchCorr_accuracy_percent']):.1f} & {float(r['clinical_accuracy_percent']):.1f} "+r'\\')
lines += [r'\midrule',r'Chance & 0.16 & 0.22 \\',r'\bottomrule',r'\end{tabular}',r'\par\smallskip{\raggedright\footnotesize Queries: 693 from 346 plant lines (BatchCorr); 536 pair means from 267 samples (clinical). WiNN minus raw: 24.8 percentage points (95\% paired interval 20.5--29.2) and 11.4 (7.3--15.7), respectively. Candidate pools are fixed across methods. QC-RFSC includes native filtering and raw-value fallback (Supplementary Methods S3). Full paired intervals: Supplementary~Table~S10.\par}',r'\end{table}']
(root/'generated/main_matching_table.tex').write_text('\n'.join(lines)+'\n')
