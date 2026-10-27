#!/usr/bin/env python3
"""Assemble vector supplementary figures from the current manuscript panels.

The figure blocks and row breaks in supplementary.tex define the panel order
and page grouping. Every panel is placed at its native PDF size; no rasterizing
or downscaling is performed. Continued figures remain multipage PDFs.

Requires pypdf. Run from any directory with:
    python3 scripts/assemble_supplementary_figures.py
"""
from collections import OrderedDict
from pathlib import Path
import re
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parent.parent
PANEL_ROOT = ROOT / "figures/supplementary_panels"
OUTPUT = ROOT / "figures/compiled_supplementary"
MM = 72.0 / 25.4
MARGIN = 3.0 * MM
GAP = 3.0 * MM


def figure_layouts():
    """Read all S1-S10 figure blocks, retaining explicit TeX row breaks."""
    layouts = OrderedDict((f"S{n}", []) for n in range(1, 11))
    tex = (ROOT / "supplementary.tex").read_text()
    blocks = re.findall(r"\\begin\{figure\*\}.*?\\end\{figure\*\}", tex, re.S)
    for block in blocks:
        graphics = block.split(r"\caption", 1)[0]
        rows = []
        for row in re.split(r"\\\\(?:\[[^\]]*\])?", graphics):
            names = re.findall(r"\{(S\d+[A-Z]_[^{}]+\.pdf)\}", row)
            if names:
                rows.append(names)
        if not rows:
            continue
        figure_id = re.match(r"S\d+", rows[0][0])[0]
        assert all(re.match(r"S\d+", name)[0] == figure_id for row in rows for name in row)
        layouts[figure_id].append(rows)
    assert all(layouts.values()), "Every supplementary figure must have a layout."
    assert sum(len(row) for pages in layouts.values() for rows in pages for row in rows) == 57
    return layouts


def native_panel(path):
    reader = PdfReader(path)
    assert len(reader.pages) == 1, path
    page = reader.pages[0]
    assert page.rotation == 0, path
    box = page.mediabox
    return page, float(box.width), float(box.height), float(box.left), float(box.bottom)


def assemble_page(writer, rows):
    panels = [[native_panel(PANEL_ROOT / name) for name in row] for row in rows]
    ncols = max(map(len, rows))
    column_widths = [max((row[j][1] for row in panels if j < len(row)), default=0) for j in range(ncols)]
    row_heights = [max(p[2] for p in row) for row in panels]
    width = 2 * MARGIN + sum(column_widths) + GAP * (ncols - 1)
    height = 2 * MARGIN + sum(row_heights) + GAP * (len(rows) - 1)
    output = writer.add_blank_page(width, height)
    top = height - MARGIN
    for row, row_height in zip(panels, row_heights):
        short_row = len(row) < ncols
        row_width = sum(p[1] for p in row) + GAP * (len(row) - 1)
        left = (width - row_width) / 2 if short_row else MARGIN
        for j, (panel, pw, ph, px, py) in enumerate(row):
            slot_width = pw if short_row else column_widths[j]
            x = left + (slot_width - pw) / 2
            y = top - ph
            assert x >= MARGIN - 1e-8 and x + pw <= width - MARGIN + 1e-8
            assert y >= MARGIN - 1e-8 and y + ph <= height - MARGIN + 1e-8
            output.merge_transformed_page(panel, Transformation().translate(x - px, y - py))
            left += slot_width + GAP
        top -= row_height + GAP
    return width / MM, height / MM


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    layouts = figure_layouts()
    paths = []
    for figure_id, pages in layouts.items():
        writer = PdfWriter()
        dimensions = [assemble_page(writer, rows) for rows in pages]
        writer.add_metadata({"/Title": f"Supplementary Figure {figure_id}",
                             "/Author": "WiNN manuscript authors",
                             "/Subject": "Vector panels; captions and methods in supplementary.tex",
                             "/Creator": "assemble_supplementary_figures.py"})
        path = OUTPUT / f"Figure_{figure_id}.pdf"
        writer.write(path)
        check = PdfReader(path)
        assert len(check.pages) == len(pages)
        expected_labels = [re.match(r"S\d+([A-Z])", name)[1] for rows in pages for row in rows for name in row]
        assert len(expected_labels) == len(set(expected_labels)), (figure_id, expected_labels)
        paths.append(path)
        dims = "; ".join(f"{w:.1f} x {h:.1f} mm" for w, h in dimensions)
        print(f"{path.relative_to(ROOT)}: {len(pages)} page(s), {dims}")
    bundle = PdfWriter()
    for path in paths:
        bundle.append(path)
    bundle.add_metadata({"/Title": "WiNN supplementary figures S1-S10", "/Creator": "assemble_supplementary_figures.py"})
    bundle_path = OUTPUT / "Supplementary_figures.pdf"
    bundle.write(bundle_path)
    assert len(PdfReader(bundle_path).pages) == sum(len(pages) for pages in layouts.values())
    print(f"{bundle_path.relative_to(ROOT)}: {len(bundle.pages)} pages")
    main_names = ["F1_winn_overview.pdf", "F2_simulation.pdf", "F3_replicates.pdf", "F4_biological_structure.pdf"]
    main_paths = [ROOT / "figures" / name for name in main_names]
    if all(path.exists() for path in main_paths):
        atlas = PdfWriter()
        for n, path in enumerate(main_paths, 1):
            atlas.append(path, outline_item=f"Figure {n}")
        for n, path in enumerate(paths, 1):
            atlas.append(path, outline_item=f"Figure S{n}")
        atlas.add_metadata({"/Title": "WiNN: main and supplementary figures", "/Creator": "assemble_supplementary_figures.py"})
        atlas_path = ROOT / "figures/WiNN_figures.pdf"
        atlas.write(atlas_path)
        check = PdfReader(atlas_path)
        assert len(check.outline) == 14
        assert len(check.pages) == sum(len(PdfReader(path).pages) for path in main_paths + paths)
        print(f"{atlas_path.relative_to(ROOT)}: {len(check.pages)} pages; 14 figure bookmarks")


if __name__ == "__main__":
    main()
