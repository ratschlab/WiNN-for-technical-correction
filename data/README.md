# Data layout

`simulated/canonical/SIM01/` contains the truth-known realization shared by the primary comparison and main WiNN ablation.

The public-data download and preprocessing scripts create the following directories:

- `public/raw/` and `public/processed/` for MTBLS79;
- `public/batchcorr_set1/`;
- `public/sacurine/`;
- `public/waveica_adenocarcinoma/`.

Processed files are not committed. Each downloader records source provenance and checksums. Raw spectra are not downloaded.

The internal cohort is not distributed. Authorized users may place its standardized prepared bundle at `private/clinical_fiams/prepared_injection_level.rds`; the entire `private/` tree is ignored by Git.
