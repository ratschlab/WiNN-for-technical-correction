#!/usr/bin/env python3
"""Check whether unavailable profile tests affect archived RFSC scalar summaries."""
import argparse
import csv
import math
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--evaluation-root',required=True,type=Path)
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args()
rows=[]
for dataset in ['mtbls79','clinical_fiams','batchcorr_set1','sacurine','waveica']:
    root=args.evaluation_root/dataset/'corrected'
    profiles=list(csv.DictReader((root/'ljung_box.csv').open()))
    def finite(value):
        try: return math.isfinite(float(value))
        except ValueError: return False
    valid=[r for r in profiles if finite(r['p_adj'])]
    count=sum(r['is_autocorrelated']=='TRUE' for r in valid)
    reported=next(r for r in csv.DictReader((root/'metrics.csv').open()) if r['metric']=='residual_ljung_box_proportion_significant')
    value=count/len(valid)
    assert int(reported['denominator'])==len(valid)
    assert abs(value-float(reported['value']))<1e-12
    rows.append(dict(dataset_key=dataset,total_profiles=len(profiles),finite_tests=len(valid),nonfinite_tests=len(profiles)-len(valid),significant_tests=count,reported_fraction=reported['value'],finite_only_fraction=value,changed=False))
with args.output.open('w') as f:
    out=csv.DictWriter(f,fieldnames=list(rows[0]));out.writeheader();out.writerows(rows)
print('PASS: all current RFSC datasets have identical reported and finite-only summaries.')
