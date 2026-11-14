WiNN manuscript package

Article: main.tex
Supplementary Methods, Results and tables: supplementary.tex
References: references.bib

In Overleaf, compile main.tex for the article. Select supplementary.tex as the
main document to compile the supplement, then restore main.tex for the article.
Both sources use the included journal class and bibliography style.

Compiled main figures: figures/F1_winn_overview.pdf, F2_simulation.pdf,
F3_replicates.pdf and F4_biological_structure.pdf. These are the same vector
figures included by main.tex. Individual panels and plotted source CSVs are
under figures/main_panels, figures/supplementary_panels and figures/source_data.
Compiled supplementary figures are under figures/compiled_supplementary.
For direct review of every figure, open figures/WiNN_figures.pdf, with bookmarks
for Figures 1-4 and Supplementary Figures S1-S10.
Captions and endpoint definitions are in the manuscript and supplement.

WiNN is the primary configuration. Supplementary comparisons distinguish
WiNN (QC-tuned) and WiNN (inferred segments); Methods define both modes.
Analysis identifiers and API values retain their original names.

Reproducible figure and table generators are in scripts/. Supporting tables,
validation results and provenance are in generated/. Private clinical matrices
and identity maps are excluded. All benchmark results use WiNN 0.1.6 (package
commit 50b72a4, frozen 11 November 2026), which corrects the mgcv smoother calls
of earlier versions; only WiNN outputs were recomputed, and raw and comparator
outputs are reused from the validated 0.1.4 release.

The supplement prints 13 condensed tables assembled by
scripts/condense_supplementary_tables.py; generated/supplementary_table_data/README.md
maps them to the detailed generator outputs, which remain available as source data.
Every numerical claim in the text was audited against these source files
before release.

Code and results snapshot (WiNN 0.1.6 analysis release):
https://github.com/ratschlab/WiNN-for-technical-correction/tree/491a6f9027f4a7ad7d6d17fc11d3b986bd5b0b13
The later documentation commit records this citation in the manuscript sources
without changing the numerical results.

Author-note revision (14 November 2026)
--------------------------------------
All 33 inline author notes were assessed; decisions are recorded in
revision_notes/author_note_decisions_20261114.csv. Reader-facing endpoint
terminology is shared by the text, figures and supplementary tables; original
CSV metric keys and filenames remain stable.

Build both PDFs with: bash scripts/build_manuscript.sh
Run numerical audits with:
  python3 revision_notes/audit_main_text_claims.py
  python3 revision_notes/audit_new_claims.py
The main-text audit now uses bundled source CSVs, with no private release paths.
The full null/preparation reruns reproduced 14,877 archived rows exactly.

Formatting was checked against Bioinformatics author guidance accessed on
14 November 2026: https://academic.oup.com/bioinformatics/pages/author-guidelines
The article uses the OUP two-column template, a structured abstract, four main
figures, figure alt text and an acknowledgement of AI assistance. The disclosure
must also be included in the journal cover letter at submission.

Updated analysis code (portable confounding validation and comparable QC-RFSC evaluation):
https://github.com/ratschlab/WiNN-for-technical-correction/tree/853541434d507e9cfe0a6a94c205b0e7d1b39f29
The numerical-results provenance remains revision 491a6f9; this code fix does
not change archived results. Package documentation is updated at d9e81b5,
while the analyzed implementation remains package commit 50b72a4.
