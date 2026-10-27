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
and identity maps are excluded. The numerical benchmark uses WiNN 0.1.4; the
current package distribution is available on CRAN.

Current code/results snapshot:
https://github.com/ratschlab/WiNN-for-technical-correction/tree/97e5ea3cadb44fe74ee914528ff19a2940fe95d7
The later documentation commit links these sources to that immutable analysis
revision; the benchmark source remains WiNN 0.1.4.
