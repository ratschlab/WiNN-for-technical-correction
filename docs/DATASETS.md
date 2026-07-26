# Benchmark datasets

## Truth-known simulation

The canonical simulation contains 1,000 features and 500 injections across five analytical batches, including 50 reference injections. Clean intensities and feature-level technical-effect annotations are known. `SIM01` is used for both the primary method comparison and the main ablation. Nine additional prespecified realizations provide a ten-seed stability analysis in total.

## MTBLS79

The benchmark uses `Dataset07__SFPM.xlsx`, the pre-normalization processed workbook deposited with MetaboLights study MTBLS79. It contains the feature-intensity matrix, eight batches, pooled QCs, acquisition order, and repeated cow/sheep sample identities. Zeros are documented as non-detections for this table. Features must be present in at least 80% of injections, and remaining missing values are imputed with missForest using seed 42.

The downloader retrieves the workbook only; raw mass-spectrometry archives are not required.

## Pair-aware clinical FIA-MS cohort

Each sample in this restricted cohort was measured as two adjacent technical injections. Correction and drift detection use the log-scale mean of each pair as one independent unit. The learned pair-level correction is applied identically to both original injections when replicate preservation is evaluated. This design prevents adjacency of technical duplicates from generating artificial lag-1 drift.

The public repository includes the generic pair-aware transform but contains no clinical data, identifiers, feature names, or row-level output. Only aggregate result tables are distributed.

## BatchCorrMetabolomics Set 1

Set 1 is loaded from `set.1` and `set.1.Y` in the pinned `rwehrens/BatchCorrMetabolomics` repository. It supplies a minimally processed negative-ion LC-MS matrix, exact sequence numbers, ten batch labels, pooled references, and accession identities with biological replication. The deposited values are already log-scaled, although the base is not documented. Imputation is performed on that scale before conversion to a nonnegative correction representation.

The principal limitation is missingness: 199 features remain after the 80% presence rule.

## Sacurine / W4M00001 / MTBLS404

The benchmark uses the negative-ion tables mirrored by the `phenomis` project. They contain a pre-correction intensity matrix, injection metadata, feature annotation, two analytical batches, and pooled QCs. Biological analyses use age, BMI, and gender. The same osmolality handling is applied after every technical-correction method.

Only 113 features pass the analysis filters, so feature coverage is reported alongside all performance metrics.

## WaveICA adenocarcinoma plasma data

The benchmark uses the XCMS-derived `Amide_data.rda` object from the WaveICA 2.0 repository. It contains a nonnegative pre-correction feature matrix, three batches, acquisition order, pooled QCs, and anonymized case/control labels. Deposited feature identifiers and group labels are retained.

This dataset is suitable with a limitation: the input is a processed XCMS feature table rather than raw peak data, and its acquisition documentation is less detailed than a full ISA-Tab deposition.

## Audited but excluded datasets

MTBLS146 and the oral-cavity dataset were excluded because the small deposited files did not provide a usable uncorrected injection-level feature-intensity matrix with the required technical metadata. Spectra were not downloaded, and no analysis was forced from already corrected or incomplete material.

## Source terms

External matrices are downloaded from their original repositories and are not committed here. BatchCorrMetabolomics declares GPL (>= 2), and the pinned phenomis source carrying the Sacurine tables declares CeCILL. The WaveICA repository includes an MIT license file, although its package DESCRIPTION contains the inconsistent value `License: No`; for that reason the downloader records the discrepancy and the matrix remains at its source. MetaboLights files remain subject to the repository's data-use terms. The download manifests retain the exact source commit or URL and checksums needed to identify every file used.
