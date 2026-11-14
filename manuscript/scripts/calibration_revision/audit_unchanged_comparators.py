#!/usr/bin/env python3
"""Record/check comparator rows so a WiNN-only revision cannot silently change them."""
import argparse, csv, gzip, hashlib, json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('mode',choices=['record','check']);p.add_argument('snapshot')
a=p.parse_args();root=Path(__file__).resolve().parents[2]
ids={'raw','combat','qc_rlsc','qc_rfsc','tiger','serrf'}
manifest={}
for path in sorted((root/'generated/revision_panel_inputs').rglob('*')):
    if not(path.name.endswith('.csv') or path.name.endswith('.csv.gz')):continue
    opener=gzip.open if path.name.endswith('.gz') else open
    with opener(path,'rt',newline='') as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or 'method_id' not in reader.fieldnames:continue
        rows=[r for r in reader if r['method_id'] in ids]
    if not rows:continue
    encoded=sorted(json.dumps(r,sort_keys=True,separators=(',',':')) for r in rows)
    manifest[str(path.relative_to(root))]={'rows':len(rows),'sha256':hashlib.sha256('\n'.join(encoded).encode()).hexdigest()}
snapshot=Path(a.snapshot)
if a.mode=='record':snapshot.write_text(json.dumps(manifest,indent=2)+'\n')
else:
    old=json.loads(snapshot.read_text())
    changed=[key for key,value in old.items() if manifest.get(key)!=value]
    if changed:raise SystemExit('Comparator rows changed:\n'+'\n'.join(changed))
print(f'{a.mode}: {len(manifest)} files; {sum(x["rows"] for x in manifest.values())} comparator rows')
