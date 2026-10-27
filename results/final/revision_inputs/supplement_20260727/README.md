# Historical aggregate inputs for supplementary table assembly

These 15 files are the exact archived numerical inputs used by the original S4–S20 assembler. Primary benchmark inputs originated from analysis base `21de7c4`; derivative aggregate files were added in the local 27 July 2026 analysis export and were not included in that public base. This release publishes the exact required inputs with uncompressed SHA-256 checksums. Clinical records contain dataset/method/metric summaries only, never patient, injection, feature or identity-level data.

These historical inputs include the previous RFSC adapter and original uncertainty. They are not current publication results. From the repository root, run:

```sh
Rscript manuscript/scripts/generate_supplementary_tables.R --final-root=results/final/revision_inputs/supplement_20260727
Rscript manuscript/scripts/nearest_baseline/generate_validation_tables.R
python3 manuscript/scripts/generate_review2_support_tables.py
```

The first command automatically runs the revised RFSC integration, which replaces simulation points, removes obsolete empirical RFSC rows and updates selected settings/S26. The next commands render the current retrieval/baseline and S29–S31 support tables from the included aggregate outputs. All clinical correction matrices and biological identifiers remain excluded.
