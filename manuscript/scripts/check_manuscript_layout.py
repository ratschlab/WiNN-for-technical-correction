#!/usr/bin/env python3
"""Audit compiled document boundaries and LaTeX diagnostics; visual review is also required."""
from pathlib import Path
import hashlib,json,re
import pdfplumber
ROOT=Path(__file__).resolve().parents[1]
report={}
for name in ('main','supplementary'):
    path=ROOT/f'{name}.pdf'
    log=(ROOT/'build/publication'/f'{name}.log').read_text(errors='replace')
    outside=[]
    with pdfplumber.open(path) as pdf:
        pages=len(pdf.pages)
        for i,page in enumerate(pdf.pages,1):
            for c in page.chars:
                if c.get('text','').strip() and (c['x0'] < -1 or c['x1'] > page.width+1 or c['top'] < -1 or c['bottom'] > page.height+1):
                    outside.append(dict(page=i,text=c['text'],bbox=[c['x0'],c['top'],c['x1'],c['bottom']]))
    errors=re.findall(r'^!.*|^.*(?:undefined references|Citation .*undefined|Reference .*undefined|Float too large).*$',log,re.M)
    overflows=re.findall(r'Overfull [^\n]+',log)
    report[name]=dict(pages=pages,pdf_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256((ROOT/f'{name}.tex').read_bytes()).hexdigest(),
        out_of_page_characters=outside,unresolved_errors=errors,overfull_boxes=overflows)
    assert not outside,(name,outside[:3])
    assert not errors,(name,errors)
assert report['main']['pages']<=7,report['main']['pages']
main=(ROOT/'main.tex').read_text()
assert '[AUTHOR NOTE:' not in main
assert '[AUTHOR NOTE:' not in (ROOT/'supplementary.tex').read_text()
report['inline_author_notes_remaining']=0
report['main_figure_count']=len(re.findall(r'\\begin\{figure\*?\}',main))
report['main_alt_text_count']=main.count(r'\textbf{Alt text:}')
report['guidance_url']='https://academic.oup.com/bioinformatics/pages/author-guidelines'
report['guidance_checked']='2026-11-14'
(ROOT/'revision_notes/final_package_qa.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:dict(pages=v['pages'],outside=len(v['out_of_page_characters']),errors=len(v['unresolved_errors']),overfull=v['overfull_boxes']) for k,v in report.items() if isinstance(v,dict)},indent=2))
