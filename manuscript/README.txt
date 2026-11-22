WiNN manuscript and reproducibility package

main.tex contains the article; supplementary.tex contains supplementary methods,
results, figures and tables. Both use the included Bioinformatics/OUP template.
Build both documents with: bash scripts/build_manuscript.sh
In Overleaf, choose main.tex or supplementary.tex as the main document.

Figures are supplied as vector PDFs with source data under figures/. The complete
bookmarked atlas is figures/WiNN_figures.pdf. figure_alt_text.txt contains four
figure descriptions for the journal's accessibility/submission fields; these are
not printed in the article.

The correction benchmarks use WiNN 0.1.6, package commit 50b72a4. The experimental
permutation gate in Supplementary Methods S10 is a separate calibration analysis
and is not used in Figure 2 or the empirical correction benchmarks.

Reproduction of the final focused analyses:
- scripts/final_review_gating_matching/README.md: paired simulation policy effects,
  null-data perturbation and empirical sample matching.
- scripts/benchmark_revision/README.md: chronological, identifier-safe QC-RFSC
  implementation, training-QC selection, abundance scaling and paired uncertainty.
- scripts/figure2_confounding/: phenotype-allocation experiments underlying
  Figure 2E-F; raw technical artifacts and WiNN configuration are fixed.

Numerical and layout checks:
  python3 scripts/audit_supplement_finalization.py
  python3 scripts/check_manuscript_layout.py
  python3 scripts/audit_figure_package.py
Input hashes, method settings, finite-sample checks, aggregate bootstrap results
and figure source CSVs are supplied under generated/ and figures/source_data/.
Clinical raw measurements and identity maps are not included.

Formatting guidance was checked on 22 November 2026:
https://academic.oup.com/bioinformatics/pages/author-guidelines
