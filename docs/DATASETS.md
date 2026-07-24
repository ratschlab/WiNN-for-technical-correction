# Benchmark datasets

## Simulation

The display simulation contains 1,000 features and 500 injections across five plates, including 50 controls. It has explicit clean ground truth and feature-level drift and batch-effect annotations. The repeated-seed robustness analysis regenerates this design with 30 independent seed ledgers.

## MTBLS79

The benchmark uses `Dataset07__SFPM.xlsx`, the pre-normalization processed workbook deposited with MetaboLights study MTBLS79. The file contains the feature-intensity matrix, eight batches, pooled QCs, acquisition order, and repeated cow/sheep sample identities. Zeros are treated as non-detections, features must be present in at least 80% of injections, and remaining values are imputed with missForest using seed 42.

The downloader pins the workbook's SHA-256 checksum and does not retrieve the much larger raw-data archives.

## BatchCorrMetabolomics Set 1

Set 1 is loaded from `set.1` and `set.1.Y` in the pinned `rwehrens/BatchCorrMetabolomics` repository. It supplies an uncorrected negative-ion LC-MS matrix, exact sequence numbers, ten batch labels, pooled references, and accession identities with genuine biological replication. The deposited values are already log-scaled, although the base is undocumented. Imputation is performed on that scale, followed by a reversible `expm1` representation for correction.

The high deposited missingness is the main limitation: the 80% presence rule leaves 199 features.

## Sacurine

The Sacurine benchmark uses the small W4M00001/MTBLS404 negative-ion tables pinned from the `phenomis` repository. The tables contain the pre-correction intensity matrix, injection metadata, feature annotation, two analytical batches, and pooled QCs. Biological analyses use age, BMI, and gender, with osmolality normalization applied consistently after each technical-correction method.

Only 113 features pass the analysis filters, so all method-specific feature loss and common-panel results are reported explicitly.

## WaveICA adenocarcinoma plasma data

The benchmark uses the XCMS-derived `Amide_data.rda` object pinned from the WaveICA 2.0 repository. It contains a nonnegative pre-correction feature matrix, three batches, acquisition order, pooled QCs, and anonymized case/control group labels. The deposited feature identifiers and group labels are retained exactly.

The matrix is suitable with limitations because it is a processed XCMS table rather than raw peak data and its package documentation is less detailed than a full ISA-Tab deposition.

## Internal clinical cohort

The sixth benchmark is a non-public FIA-MS cohort. Each biological sample was injected as two adjacent technical dips. Correction and run-order diagnostics therefore operate on one independent log-mean unit per adjacent pair, and the fitted correction is then mapped back to both injections for replicate-agreement evaluation. This prevents the paired injection design from appearing as artificial lag-1 drift.

No clinical inputs, identifiers, feature names, row-level outputs, or cohort-specific loader code are distributed here. Only aggregate, manuscript-facing metrics are included in the separate analysis-results release.

## Audited but excluded candidates

MTBLS146 and the oral-cavity study were audited and rejected because a usable uncorrected injection-level feature-intensity matrix could not be recovered from the small deposited files. Their spectra were not downloaded, and no correction results were manufactured from already corrected values or incomplete metadata.
