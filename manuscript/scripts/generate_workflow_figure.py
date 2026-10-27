#!/usr/bin/env python3
"""Draw the WiNN workflow at its final 175 mm print width.

The four deterministic illustrative profiles depict the selection logic only;
no empirical observations or benchmark estimates enter this schematic.
"""
from pathlib import Path
import math
import random
import reportlab
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor

ROOT = Path(__file__).resolve().parents[1]
MM = 72 / 25.4
W, H = 175 * MM, 68 * MM
FONTS = Path(reportlab.__file__).parent / 'fonts'
pdfmetrics.registerFont(TTFont('WorkflowSans', str(FONTS / 'Vera.ttf')))
pdfmetrics.registerFont(TTFont('WorkflowSansBold', str(FONTS / 'VeraBd.ttf')))
OUT = ROOT / 'figures/F1_winn_overview.pdf'
c = canvas.Canvas(str(OUT), pagesize=(W, H), pageCompression=1)
c.setTitle('WiNN workflow: randomized acquisition and selective correction')
c.setAuthor('WiNN manuscript authors')
INK = '#263247'; GRID = '#CBD3DC'; BLUE = '#0072B2'; RED = '#B44C46'


def label(x, y, value, bold=False, align='left', size=9, color=INK):
    c.setFillColor(HexColor(color))
    c.setFont('WorkflowSansBold' if bold else 'WorkflowSans', size)
    draw = {'left': c.drawString, 'center': c.drawCentredString, 'right': c.drawRightString}[align]
    draw(x, y, value)


def line(x1, y1, x2, y2, color=INK, width=.5, dash=None):
    c.setStrokeColor(HexColor(color)); c.setLineWidth(width)
    c.setDash(dash or [])
    c.line(x1, y1, x2, y2)
    c.setDash([])


def arrow(x1, x2, y):
    line(x1, y, x2, y, '#8894A5', .65)
    line(x2 - 2.5, y + 1.8, x2, y, '#8894A5', .65)
    line(x2 - 2.5, y - 1.8, x2, y, '#8894A5', .65)


# Equal biological-group counts in each batch, with randomized within-batch order.
label(5, 181, 'Randomized acquisition', bold=True)
label(5, 167, 'Balance groups across batches;')
label(5, 155, 'mix groups along run order.')
group_colors = ['#0072B2', '#D55E00', '#009E73']
for batch, x0 in enumerate([208, 280, 352, 424], 1):
    label(x0 + 26, 181, f'Batch {batch}', align='center')
    order = [0, 0, 1, 1, 2, 2]
    random.Random(20 + batch).shuffle(order)
    for j, group in enumerate(order):
        c.setFillColor(HexColor(group_colors[group])); c.setStrokeColor(HexColor(group_colors[group]))
        c.circle(x0 + j * 10, 167, 3.2, fill=1, stroke=0)
    arrow(x0 - 3, x0 + 54, 157)
line(5, 149, W - 5, 149, GRID)

# Four rows distinguish whether drift and batch shifts are selected.
xs = [54, 165, 276, 387]
pw = 102
heads = [('Input feature', 'profiles'), ('Serial dependence?', 'GAM detrending'),
         ('Residual batch shift?', 'Batch correction'), ('Corrected feature', 'profiles')]
for x, (a, b) in zip(xs, heads):
    label(x + pw/2, 135, a, align='center', bold=False)
    label(x + pw/2, 123, b, align='center', bold=True)
ys = [107, 81, 55, 29]
case_names = [('Neither', 'effect'), ('Drift', ''), ('Batch', ''), ('Both', '')]
for y, (a, b) in zip(ys, case_names):
    label(5, y + (2 if b else -3), a)
    if b: label(5, y - 9, b)

rng = random.Random(107)
n = 96
noise = [rng.gauss(0, .055) for _ in range(n)]
shift = [0, .24, -.18, .31]
drift = [.22 * math.sin(2*math.pi*(j%24)/23) + .13*(j%24)/23 - .065 for j in range(n)]


def profile(x, y, values, fitted=None, fit_color=BLUE):
    height = 18
    c.setStrokeColor(HexColor('#A6AEB9')); c.setLineWidth(.35)
    c.rect(x, y-height/2, pw, height, fill=0, stroke=1)
    for k in (1, 2, 3):
        line(x + pw*k/4, y-height/2, x + pw*k/4, y+height/2, GRID, .35, [2,2])
    def px(j): return x + 2 + (pw - 4) * j/(n-1)
    def py(v): return y + 15*v - 1
    c.setFillColor(HexColor('#4C5564'))
    for j, value in enumerate(values): c.circle(px(j), py(value), .45, fill=1, stroke=0)
    if fitted is not None:
        for k in range(4):
            for j in range(k*24, (k+1)*24-1):
                line(px(j), py(fitted[j]), px(j+1), py(fitted[j+1]), fit_color, .75)


for row, y in enumerate(ys):
    has_drift = row in (1, 3)
    has_batch = row in (2, 3)
    offsets = [shift[j//24] if has_batch else 0 for j in range(n)]
    trend = [drift[j] if has_drift else 0 for j in range(n)]
    raw = [noise[j]+offsets[j]+trend[j] for j in range(n)]
    profile(xs[0], y, raw)
    if has_drift:
        profile(xs[1], y, raw, [offsets[j]+trend[j] for j in range(n)], BLUE)
    else: label(xs[1]+pw/2, y-3, 'Unchanged', align='center', color='#667487')
    if has_batch:
        profile(xs[2], y, [noise[j]+offsets[j] for j in range(n)], offsets, RED)
    else: label(xs[2]+pw/2, y-3, 'Unchanged', align='center', color='#667487')
    profile(xs[3], y, [value + sum(offsets)/n for value in noise])
    for j in range(3): arrow(xs[j]+pw+1, xs[j+1]-2, y)
label(5, 5, 'Illustrative profiles: injection order runs left to right; dashed lines mark batches.', size=9)
c.showPage(); c.save()
print(OUT.relative_to(ROOT), '175 x 68 mm; all text 9 pt')
