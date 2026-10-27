#!/usr/bin/env python3
"""Assemble publication figures directly from the vector panel PDFs.

Panels retain their native 85 mm width; each two-column figure is 175 mm wide.
No data are reread or recalculated, and no rasterization is used.
Run with Python and pypdf from any working directory.
"""
from pathlib import Path
import csv
import hashlib
from pypdf import PageObject, PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / 'figures'
PANELS = FIGURES / 'main_panels'
MM = 72 / 25.4
SPECS = [
    ('F2_simulation.pdf', [
        'F2A_selective_forced_truth_rmse.pdf', 'F2B_truth_log1p_rmse.pdf',
        'F2C_truth_metabolite_profile_icc_mean.pdf', 'F2D_truth_sample_profile_pearson_mean.pdf',
        'F2E_median_attenuation_ratio_clean_responsive.pdf', 'F2F_empirical_fdr.pdf']),
    ('F3_replicates.pdf', [
        'F3A_batchcorr_qc_cv_vs_run_order_deviance.pdf', 'F3B_clinical_qc_cv_vs_run_order_deviance.pdf',
        'F3C_batchcorr_batch_variance_vs_replicate_icc.pdf', 'F3D_clinical_plate_variance_vs_replicate_icc.pdf']),
    ('F4_biological_structure.pdf', [
        'F4A_waveica_batch_vs_group.pdf', 'F4B_mtbls79_batch_vs_group.pdf']),
]
manifest = []
for filename, panel_names in SPECS:
    rows = len(panel_names) // 2
    height_mm = rows * 55 + (rows - 1) * 3
    page = PageObject.create_blank_page(width=175 * MM, height=height_mm * MM)
    for index, name in enumerate(panel_names):
        source = PANELS / name
        reader = PdfReader(source)
        assert len(reader.pages) == 1
        panel = reader.pages[0]
        assert abs(float(panel.mediabox.width) - 85 * MM) < 1.1
        assert abs(float(panel.mediabox.height) - 55 * MM) < 1.1
        scale = 85 * MM / float(panel.mediabox.width)
        row, column = divmod(index, 2)
        x = column * 90 * MM
        y = (height_mm - (row + 1) * 55 - row * 3) * MM
        y += 55 * MM - float(panel.mediabox.height) * scale
        page.merge_transformed_page(panel, Transformation().scale(scale).translate(x, y))
        manifest.append(dict(figure=filename, panel=name,
            panel_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            width_mm=85, height_mm=55, x_mm=x/MM, y_mm=y/MM))
    writer = PdfWriter()
    writer.add_page(page)
    writer.add_metadata({'/Title':filename.removesuffix('.pdf').replace('_',' '),
                         '/Creator':'WiNN reproducible vector figure assembly'})
    writer.write(FIGURES / filename)

bundle = PdfWriter()
for filename in ['F1_winn_overview.pdf'] + [spec[0] for spec in SPECS]:
    bundle.append(FIGURES / filename)
bundle.add_metadata({'/Title':'WiNN main figures', '/Creator':'WiNN reproducible vector figure assembly'})
bundle.write(FIGURES / 'Main_figures.pdf')
with (ROOT / 'generated/main_figure_assembly_manifest.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(manifest[0]))
    writer.writeheader()
    writer.writerows(manifest)
print('Assembled F2-F4 and the four-page main-figure bundle without rasterization.')
