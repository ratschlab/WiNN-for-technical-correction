#!/usr/bin/env python3
"""Validate figure references, displayed labels and vector/PNG provenance.

Uses pypdf and Pillow. Tesseract OCR is used when available. This complements
visual review; it does not infer aesthetic quality from text extraction.
"""
from pathlib import Path
import csv
import hashlib
import re
import shutil
import subprocess
from pypdf import PdfReader
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OLD_TERMINOLOGY = re.compile(r'clean[ -](?:truth|reference)|to clean|held-out reference CV|deposited groups|weighted[ -]PC R|run-order deviance', re.I)

OLD = re.compile(r'\bwinn\s*(?:fixed|auto)|\bfixed[\s-]+winn\b', re.I)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def normalize(text):
    return text.replace('\u2011','-').replace('\u2013','-').replace('\u2014','-')

def main():
    checks=[]
    def add(name, path, passed, detail=''):
        checks.append(dict(check=name, path=str(path), passed=bool(passed), detail=detail))
    for texname in ('main.tex','supplementary.tex'):
        tex=(ROOT/texname).read_text()
        refs=re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^{}]+\.pdf)\}',tex)
        refs+=['figures/supplementary_panels/'+name for name in
               re.findall(r'\\supppanel\{[^{}]+\}\{([^{}]+\.pdf)\}',tex)]
        for ref in refs:
            candidates=[ROOT/ref, ROOT/'figures'/ref, ROOT/'figures/main_panels'/ref, ROOT/'figures/supplementary_panels'/ref]
            resolved=next((p for p in candidates if p.is_file()),None)
            add('referenced PDF exists',ref,resolved is not None,texname + (' -> '+str(resolved.relative_to(ROOT)) if resolved else ''))
    pdfs=sorted((ROOT/'figures').rglob('*.pdf'))
    for path in pdfs:
        reader=PdfReader(path)
        txt=normalize('\n'.join(p.extract_text() or '' for p in reader.pages))
        matches=OLD.findall(txt)
        old_terms=OLD_TERMINOLOGY.findall(txt)
        add('current displayed endpoint terminology',path.relative_to(ROOT),not old_terms,
            '; '.join(old_terms))
        add('no obsolete displayed method label',path.relative_to(ROOT),not matches,
            f'{len(reader.pages)} pages; {len(matches)} obsolete labels')
    for path in sorted((ROOT/'figures/source_data').glob('*.csv')):
        rows=list(csv.DictReader(path.open()))
        hits=[(i+2,key,value) for i,row in enumerate(rows) for key,value in row.items()
              if value and OLD.search(normalize(value))]
        add('no obsolete source-data display label',path.relative_to(ROOT),not hits,
            f'{len(rows)} rows; {len(hits)} obsolete labels')
    manifest={r['png']:r for r in csv.DictReader((ROOT/'generated/figure_preview_manifest.csv').open())}
    ocr=shutil.which('tesseract')
    for path in sorted((ROOT/'figures').rglob('*.png')):
        rel=str(path.relative_to(ROOT)); record=manifest.get(rel)
        add('PNG has provenance entry',rel,record is not None)
        if record:
            pdf=ROOT/record['source_pdf']
            add('PNG source PDF unchanged',rel,digest(pdf)==record['source_pdf_sha256'])
            add('PNG bytes match rendered preview',rel,digest(path)==record['png_sha256'])
            with Image.open(path) as im:
                add('PNG readable',rel,im.width>0 and im.height>0,f'{im.width}x{im.height}')
        if ocr:
            result=subprocess.run([ocr,str(path),'stdout','--psm','11'],capture_output=True,text=True,check=True)
            hits=OLD.findall(normalize(result.stdout))
            add('no obsolete raster method label',rel,not hits,f'{len(hits)} obsolete labels')
    atlas=PdfReader(ROOT/'figures/WiNN_figures.pdf')
    titles=[str(item['/Title']) for item in atlas.outline if isinstance(item,dict)]
    add('figure atlas bookmarks','figures/WiNN_figures.pdf',len(titles)==14,'; '.join(titles))
    out=ROOT/'revision_notes/figure_package_audit.csv'
    with out.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(checks[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(checks)
    failed=[r for r in checks if not r['passed']]
    for row in failed:print('FAILED:',row)
    print(f'Audited {len(pdfs)} vector PDFs and {len(manifest)} raster previews; '
          f'{len(checks)} checks, {len(failed)} failed.')
    raise SystemExit(bool(failed))

if __name__=='__main__':
    main()
