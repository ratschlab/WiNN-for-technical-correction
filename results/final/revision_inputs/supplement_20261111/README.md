# Aggregate inputs for supplementary table assembly (WiNN 0.1.6, 11 November 2026)

These 15 files have the same layout as `supplement_20260727` and replace it for the current tables. They come from the WiNN 0.1.6 re-run (frozen source tarball SHA-256 `e40a2c3fd709f11a456ac3721379da3549e89c619e5a164743e203fe1e9be7ab`). In 0.1.6 the mgcv smoothers are called correctly; earlier versions fell back to LOESS.

Only WiNN was recomputed. The Raw and competitor rows (ComBat, QC-RLSC, QC-RFSC, TIGER, SERRF) are copied unchanged from the validated July release by `analysis/splice_winn_rerun.R`. As in the July set, the QC-RFSC rows are the previous-adapter values; the supplement assembler replaces them with the corrected results automatically.

The derivative summaries come from the July derivative code. In it, only the three WiNN modes were recomputed; all other methods' rows reuse July's draws, which come from the same fixed seeds and common panel, so the paired comparisons remain exactly paired. `validation/release_validation.csv` combines two validations:

- the WiNN-only re-run, from `analysis/validate_release.R` with `WINN_RERUN_METHODS` set;
- the validation of the reused competitor outputs, from the July release.

`analysis/combine_winn_rerun_validation.R` merges them.

From the manuscript directory:

```sh
Rscript scripts/generate_supplementary_tables.R --final-root=../public_analysis_release/results/final/revision_inputs/supplement_20261111
Rscript scripts/nearest_baseline/generate_validation_tables.R
python3 scripts/generate_review2_support_tables.py
Rscript scripts/generate_s24_paired_effects_table.R
```

Clinical records are dataset/method/metric summaries only. No patient, injection, feature or identity-level data are included. SHA-256 checksums in `input_manifest.csv` refer to the exact uncompressed bytes.
