#!/usr/bin/env python3
"""Reassemble Figure 2 only, then refresh bundles containing that figure."""
from pathlib import Path
import csv, hashlib
from pypdf import PageObject, PdfReader, PdfWriter, Transformation
ROOT=Path(__file__).resolve().parents[2]; F=ROOT/'figures'; MM=72/25.4
names=['F2A_selective_forced_truth_rmse.pdf','F2B_truth_log1p_rmse.pdf',
       'F2C_truth_metabolite_profile_icc_mean.pdf','F2D_truth_sample_profile_pearson_mean.pdf',
       'F2E_median_attenuation_ratio_clean_responsive.pdf','F2F_empirical_fdr.pdf']
page=PageObject.create_blank_page(width=175*MM,height=156*MM)
rows=[]
for i,name in enumerate(names):
    p=F/'main_panels'/name; q=PdfReader(p).pages[0]; scale=85*MM/float(q.mediabox.width)
    x=(i%2)*90*MM; y=(156-(i//2+1)*50-(i//2)*3)*MM+50*MM-float(q.mediabox.height)*scale
    page.merge_transformed_page(q,Transformation().scale(scale).translate(x,y))
    rows.append(dict(figure='F2_simulation.pdf',panel=name,panel_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                     width_mm=85,height_mm=50,x_mm=x/MM,y_mm=y/MM))
w=PdfWriter();w.add_page(page);w.add_metadata({'/Title':'F2 simulation','/Creator':'WiNN Figure 2 confounding-only assembly'});w.write(F/'F2_simulation.pdf')
main=['F1_winn_overview.pdf','F2_simulation.pdf','F3_replicates.pdf','F4_biological_structure.pdf']
w=PdfWriter()
for name in main:w.append(F/name)
w.add_metadata({'/Title':'WiNN main figures'});w.write(F/'Main_figures.pdf')
# Supplement-only comparison retained at the user's request; existing S1–S10 stay intact.
calibration=['S11A_null_calibration.pdf','S11B_correlated_calibration.pdf','S11C_drift_sensitivity.pdf']
if all((F/'supplementary_panels'/name).is_file() for name in calibration):
    cp=PageObject.create_blank_page(width=174*MM,height=176*MM)
    for i,name in enumerate(calibration):
        q=PdfReader(F/'supplementary_panels'/name).pages[0]
        sc=174*MM/float(q.mediabox.width)
        cp.merge_transformed_page(q,Transformation().scale(sc).translate(0,(176-i*59)*MM-float(q.mediabox.height)*sc))
    cw=PdfWriter();cw.add_page(cp);cw.write(F/'compiled_supplementary/Figure_S11.pdf')
sw=PdfWriter()
for i in range(1,12):sw.append(F/'compiled_supplementary'/f'Figure_S{i}.pdf')
sw.add_metadata({'/Title':'WiNN supplementary figures S1-S11'});sw.write(F/'compiled_supplementary/Supplementary_figures.pdf')
w=PdfWriter()
for i,name in enumerate(main,1):w.append(F/name,outline_item=f'Figure {i}')
for i in range(1,12):w.append(F/'compiled_supplementary'/f'Figure_S{i}.pdf',outline_item=f'Figure S{i}')
w.add_metadata({'/Title':'WiNN: main and supplementary figures'});w.write(F/'WiNN_figures.pdf')
p=ROOT/'generated/main_figure_assembly_manifest.csv'
with p.open() as f: old=list(csv.DictReader(f))
rows+= [r for r in old if r['figure']!='F2_simulation.pdf']
with p.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
print('Assembled Figure 2 and containing bundles; all other individual figures untouched.')

p=ROOT/'generated/revision_main_panel_manifest.csv'
s=p.read_text().replace('F2,E,WiNN: effect recovery,','F2,E,WiNN: retained group contrast,').replace('F2,F,WiNN: empirical FDR,','F2,F,WiNN: detection and false discoveries,')
p.write_text(s)
