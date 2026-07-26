# Suitability audit: BatchCorrMetabolomics Set 1

**Verdict: SUITABLE WITH LIMITATIONS — PROCEED TO THE CORRECTION BENCHMARK.**

The package supplies an uncorrected LC-MS metabolite-intensity matrix, exact batch and acquisition sequence, pooled references, and accession identities with genuine biological replicates. The matrix is suitable for the requested benchmark after documented missingness filtering and imputation.

## Provenance

- Repository: [`rwehrens/BatchCorrMetabolomics`](https://github.com/rwehrens/BatchCorrMetabolomics), pinned commit `e0c7668140e206dcdae9afa602dd2e1b337ac4f6`.
- Source object: `data/BC.RData`, SHA-256 `a5d3918c69a902af3886c0141292a614cc20c01526a68930853157e5c9aec113`.
- Objects: `set.1` and `set.1.Y`.
- Paper: Wehrens et al., “Improved batch correction in untargeted MS-based metabolomics,” *Metabolomics* 12:88 (2016), DOI [10.1007/s11306-016-1015-8](https://doi.org/10.1007/s11306-016-1015-8).
- Assay: negative-ion untargeted LC-MS of Arabidopsis natural accessions.

Only the 3.3 MB package data object, documentation, demo, and article XML were downloaded. No raw spectra were needed or downloaded.

## Matrix and processing state

`set.1` is a numeric matrix with 761 injections as rows and 567 reconstructed metabolites as columns. The package documentation states that Metalign and MSClust were used for mass-feature extraction, alignment, and clustering. The paper states that metabolites occurring in fewer than 20 genotypes were removed. These are acceptable upstream feature-construction steps rather than correction.

The paper and figures explicitly call Set I “uncorrected,” and the package demo loads `set.1` before applying its correction models. There is no evidence of prior run-order correction, batch correction, or QC normalization.

The deposited matrix has 206,376 NAs (47.83%) and no zeros. The paper documents the missing values as non-detects. Applying the benchmark’s 80%-observation rule retains 199 features with 3.11% missingness before imputation; 368 features are removed. Remaining values are imputed with missForest (`set.seed(42)`, 100 trees) on the deposited log scale.

## Design reconstruction

`set.1.Y` has one row per matrix row and three fields: `SeqNr`, `Batch`, and `SCode`.

- `SeqNr` is globally unique and is used as the acquisition order.
- `Batch` defines 10 analytical batches.
- `SCode == "ref"` identifies the pooled reference QC; all other values are anonymized accession IDs.
- Matrix and metadata row names are identical and unique.

The chronological batch sequence, recovered from minimum `SeqNr` rather than from the batch names, is:

`B1, B2, B9, B10, B7, B8, B3, B4, B5, B6`.

| Batch | Total | Study | QC | SeqNr range | Maximum study injections between QCs |
|---|---:|---:|---:|---:|---:|
| B1 | 79 | 74 | 5 | 10–88 | 19 |
| B2 | 78 | 72 | 6 | 89–171 | 19 |
| B9 | 79 | 74 | 5 | 183–261 | 19 |
| B10 | 48 | 43 | 5 | 262–309 | 11 |
| B7 | 79 | 74 | 5 | 322–401 | 19 |
| B8 | 80 | 75 | 5 | 402–494 | 19 |
| B3 | 79 | 74 | 5 | 507–585 | 19 |
| B4 | 79 | 74 | 5 | 586–665 | 19 |
| B5 | 80 | 75 | 5 | 678–757 | 19 |
| B6 | 80 | 75 | 5 | 758–837 | 19 |

Every batch begins and ends with a QC in the retained injection table. A fixed holdout of one interior QC per batch is therefore feasible while retaining four or five training QCs and both QC endpoints in every batch.

## Biological replication and balance

There are 710 study injections from 357 accession IDs:

- 351 accessions have two injections;
- one accession has three injections;
- five accessions have one injection;
- 346 accessions have measurements in more than one batch.

This is strong genuine biological replication for accession-level agreement and feature repeatability. Metrics that require at least two measurements exclude the five singletons. The 11 accessions not spanning batches are reported rather than treated as cross-batch replicates.

## Limitations

1. The package demo explicitly says Set 1 is already log-scaled, but neither the documentation nor the paper states the log base. Imputation is therefore done directly on the deposited scale. A reversible `expm1` linearization is used only to feed methods defined for positive multiplicative intensities; all evaluations use `log1p`, exactly returning to the deposited/imputed analysis scale.
2. The anonymized matrix has no metabolite column names, m/z values, or retention times. Stable IDs `Feature_0001` through `Feature_0567` preserve original column positions.
3. The original matrix is sparse. The 80% rule reduces the benchmark to 199 features, so conclusions apply to the well-observed portion of Set 1.
4. `SeqNr` contains gaps between or within several batches. These are retained as documented acquisition numbers, not renumbered globally; `within_batch_order` is a separate rank field.

## QC-RLSC implementation safeguard

Set 1 retains only four or five training QCs per batch after the common holdout, and its deposited values are already log-scaled. The benchmark therefore uses a dataset-specific, sample-ID-preserving QC-RLSC implementation: subtractive per-batch LOESS on the exact `log1p` analysis scale with fixed `span = 1` and `degree = 1`, followed by batch shifting on that log scale. Parameters are fixed from the QC design and never selected using hidden-QC values.

The wrapper preserves input sample order explicitly and forwards the fixed span to each within-batch LOESS fit. These safeguards are checked before the corrected matrix is accepted.

## Audit artifacts

- Download script: `scripts/download_batchcorr_set1.py`
- Preprocessing script: `scripts/preprocess_batchcorr_set1.R`
- Download manifest: `data/public/batchcorr_set1/source/download_manifest.csv`
- Standardized metadata: `data/public/batchcorr_set1/processed/BatchCorr_set1_metadata.csv`
- Feature audit: `data/public/batchcorr_set1/processed/BatchCorr_set1_feature_metadata.csv`
- Batch/QC counts: `data/public/batchcorr_set1/audit/batch_and_qc_counts.csv`
- QC positions: `data/public/batchcorr_set1/audit/qc_positions.csv`
- Accession/batch balance: `data/public/batchcorr_set1/audit/accession_batch_balance.csv`
- Machine-readable criterion report: `reports/batchcorr_set1_suitability.csv`
