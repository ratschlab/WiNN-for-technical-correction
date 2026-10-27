#!/usr/bin/env python3
"""Render every packaged PNG preview from its current vector PDF.

Run after regenerating panels and assembling figures. Poppler's pdftoppm is
required. The manifest binds each PNG to the exact source PDF bytes.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import hashlib
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DPI = 300

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def render(path):
    pdf = path.with_suffix('.pdf')
    assert pdf.is_file(), f'PNG has no corresponding vector source: {path}'
    subprocess.run(['pdftoppm', '-png', '-singlefile', '-r', str(DPI),
                    str(pdf), str(path.with_suffix(''))], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return dict(png=str(path.relative_to(ROOT)), source_pdf=str(pdf.relative_to(ROOT)),
                dpi=DPI, source_pdf_sha256=digest(pdf), png_sha256=digest(path))

def main():
    paths = sorted((ROOT / 'figures').rglob('*.png'))
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(render, paths))
    out = ROOT / 'generated/figure_preview_manifest.csv'
    with out.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print(f'Rendered {len(rows)} PNG previews from current vector PDFs at {DPI} dpi.')

if __name__ == '__main__':
    main()
