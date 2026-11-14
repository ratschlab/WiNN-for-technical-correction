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
for Figures 1-4 and Supplementary Figures S1-S11.
Captions and endpoint definitions are in the manuscript and supplement.

WiNN is the primary configuration. Supplementary comparisons distinguish
WiNN (QC-tuned) and WiNN (inferred segments); Methods define both modes.
Analysis identifiers and API values retain their original names.

Reproducible figure and table generators are in scripts/. Supporting tables,
validation results and provenance are in generated/. Private clinical matrices
and identity maps are excluded. WiNN corrections use development version 0.1.7.9999
(package commit ab4969792c98bd86208a28a4e83b7a69e204b9fc), which replaces the
asymptotic drift gate with 999 shared acquisition-order permutations and maxT
familywise adjustment. Smoothing, MAD preparation, batch adjustment and PQN
are unchanged. Only WiNN outputs were recomputed; raw and comparator outputs
are preserved, including the corrected QC-RFSC evaluation. CRAN version 0.1.6
does not include the new calibration.

The supplement prints 13 condensed tables assembled by
scripts/condense_supplementary_tables.py; generated/supplementary_table_data/README.md
maps them to the detailed generator outputs, which remain available as source data.
Every numerical claim in the text was audited against these source files
before release.

Code and results snapshot (permutation-calibrated WiNN development release):
https://github.com/ratschlab/WiNN-for-technical-correction/tree/3794b8408070950cd59c0eed9fe207c6c8fbb5f6
The frozen supplement_20261114 and derivative_20261114 archives contain the
aggregate inputs used to regenerate the publication tables and figures.
The final manuscript, table/figure generators and publication checks are at:
https://github.com/ratschlab/WiNN-for-technical-correction/tree/main/manuscript

Author-note revision (14 November 2026)
--------------------------------------
All 33 inline author notes were assessed. Reader-facing endpoint
terminology is shared by the text, figures and supplementary tables; original
CSV metric keys and filenames remain stable.

Build both PDFs with: bash scripts/build_manuscript.sh
Run numerical audits with:
  python3 revision_notes/audit_main_text_claims.py
  python3 revision_notes/audit_new_claims.py
The main-text audit now uses bundled source CSVs, with no private release paths.
The legacy null/preparation verification reproduced 14,877 archived rows exactly.

Formatting was checked against Bioinformatics author guidance accessed on
14 November 2026: https://academic.oup.com/bioinformatics/pages/author-guidelines
The article uses the OUP two-column template, a structured abstract, four main
figures, figure alt text and an acknowledgement of AI assistance. The disclosure
must also be included in the journal cover letter at submission.
