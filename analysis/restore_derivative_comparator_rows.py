#!/usr/bin/env python3
"""Verify recomputed derivative comparator values, then retain their frozen rows."""
import csv, math, argparse
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--new',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
a=p.parse_args();checks=[]
def read(path):
    with path.open(newline='') as f:
        reader=csv.DictReader(f);return reader.fieldnames,list(reader)
def close(x,y):
    if x==y:return True,0.
    try:
        u,v=float(x),float(y)
        if math.isnan(u) and math.isnan(v):return True,0.
        return math.isclose(u,v,rel_tol=1e-9,abs_tol=1e-9),abs(u-v)
    except (ValueError,TypeError):return False,0.
for target in sorted(a.new.glob('*.csv')):
    # LOBO runs independently and already splices its frozen comparator rows.
    if target.name.startswith('leave_one_batch_out_bootstrap'):continue
    old=a.previous/target.name
    if not old.exists():continue
    fields,rows=read(target);old_fields,prior=read(old)
    if 'method_id' not in fields:continue
    keep=lambda r:not r['method_id'].startswith('winn_')
    actual=list(filter(keep,rows));expected=list(filter(keep,prior))
    if not expected:continue
    assert set(fields)==set(old_fields),(target.name,'column mismatch')
    assert len(actual)==len(expected),(target.name,'comparator row count')
    keys=[k for k in ('method_id','metric','endpoint','feature_id','identity','bootstrap',
                     'permutation','heldout_batch','feature_panel','seed_id') if k in fields]
    key=lambda row:tuple(row[k] for k in keys)
    actual_by_key={key(row):row for row in actual}
    expected_by_key={key(row):row for row in expected}
    assert len(actual_by_key)==len(actual),(target.name,'duplicate actual comparator key',keys)
    assert len(expected_by_key)==len(expected),(target.name,'duplicate historical comparator key',keys)
    assert actual_by_key.keys()==expected_by_key.keys(),(target.name,'comparator key set differs',keys)
    worst=0.
    for row_key,x in actual_by_key.items():
        y=expected_by_key[row_key]
        for k in fields:
            ok,delta=close(x[k],y[k]);worst=max(worst,delta)
            assert ok,(target.name,row_key,k,'comparator mismatch')
    # Partial reruns can group rows by method; match exact scientific identifiers,
    # including the bootstrap/permutation index, rather than relying on row order.
    output=[expected_by_key[key(row)] if keep(row) else row for row in rows]
    with target.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(output)
    checks.append({'file':target.name,'rows':len(expected),'max_numeric_difference_before_restore':worst,'passed':True})
with (a.new/'comparator_row_restoration.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['file','rows','max_numeric_difference_before_restore','passed']);w.writeheader();w.writerows(checks)
print('PASS:',len(checks),'derivative files retain verified frozen comparator rows')
