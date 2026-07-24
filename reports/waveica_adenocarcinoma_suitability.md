# WaveICA adenocarcinoma suitability audit

**Verdict: SUITABLE WITH LIMITATIONS.** `Amide_data` provides an aligned XCMS-derived, pre-correction intensity matrix with global injection order, three batches, pooled QCs, and two balanced-within-batch biological groups. The benchmark may proceed.

## Source and provenance

Pinned source: `https://github.com/dengkuistat/WaveICA_2.0` at commit `56fff2e6b0410b5957c6ea83bc658241df222f41`; object `Amide_data`.
Only the 33 MB R data object and small package source/documentation were downloaded. No mzXML, mzML, vendor, or raw spectra were downloaded.

## Design and scale

- Matrix: 6461 peaks × 642 injections; 568 study samples and 74 QCs.
- Batches: 192 study + 25 QC, 192 + 25, and 184 + 24.
- Study groups by batch: group 1 = 168/168/161; group 0 = 24/24/23. All batches contain both groups in essentially identical proportions.
- Injection order is globally unique 1–642 and batches are contiguous (1–217, 218–434, 435–642).
- Intensities range from 0 to 2.34e+08 (median 3.48e+04); 1.177% of entries are zero and no values are missing or nonfinite.
- Zeros are retained because the deposit does not establish that zero is an NA placeholder. No missingness filtering or imputation was needed.

## Pre-correction evidence and label mapping

The package source labels `Amide_data` as original data, removes only its first three metadata columns, and passes the remaining peak table directly to WaveICA 2.0. The associated paper describes XCMS preprocessing followed by batch-effect correction. Raw-intensity magnitude, zeros, and visible batch/run-order structure are consistent with this provenance; there is no evidence of prior drift or batch correction.
The object stores groups only as `0` and `1`. It contains exactly 497 group-1 and 71 group-0 study samples, while the publication reports cohort totals of 497 colorectal-cancer and 71 chronic-enteritis patients. This count agreement does not prove the numeric label direction, so disease names are not assigned and all analyses use the neutral deposited groups.

## Limitations

The source lacks m/z/retention-time annotations, a complete XCMS parameter record, and a source-backed mapping from numeric groups to diagnoses. Group-batch Cramer's V is 0.0000 (p=1); sequence diagnostics and logistic models are saved under `data/public/waveica_adenocarcinoma/audit/`.
Patients are not repeated across batches, so no replicate correlation or patient ICC will be calculated.

## Criterion ledger

- 1. **PASS** — Numerical feature-by-injection matrix: 6461 numeric peak rows by 642 aligned injections after orientation.
- 2. **PASS WITH LIMITATION** — Uncorrected or minimally processed intensities: The paper/repository describe an XCMS-derived peak table; m/z and retention-time annotations and the complete XCMS script are not included.
- 3. **PASS WITH LIMITATION** — No prior batch/drift/QC correction: The package example treats Amide_data as original input to WaveICA; values are raw-intensity-like with zeros and strong batch/drift structure, with no evidence of prior correction.
- 4. **PASS** — Study and pooled QC in same matrix: 568 study samples and 74 pooled QCs coexist in Amide_data.
- 5. **PASS** — Unambiguous matrix-metadata matching: Metadata and features are columns of the same deposited data frame; row names are unique sample IDs.
- 6. **PASS** — Reliable supplied batch labels: Three supplied batches occupy contiguous acquisition intervals.
- 7. **PASS** — Reliable acquisition order: Injection_order is globally unique 1..642 and equals deposited row order.
- 8. **PASS** — Identifiable QCs: group=QC identifies 25, 25, and 24 QCs by batch.
- 9. **PASS WITH LIMITATION** — Usable biological labels: Neutral groups are direct. Disease labels are not assigned: group counts match the two published cohort totals but do not prove which numeric value denotes which diagnosis.
- 10. **PASS** — Sufficient study samples per batch: 192, 192, and 184 study samples in batches 1, 2, and 3.
- 11. **PASS** — Sufficient training QCs per batch: Five hidden QCs per batch leave 20, 20, and 19 training QCs.
- 12. **PASS** — Boundary QCs can remain in training: The fixed holdout excludes the first and last QC of each batch.
- 13. **PASS** — No severe biological-batch/order confounding: Both groups occur in every batch in virtually identical proportions; Cramer's V=0.0000 and chi-squared p=1.
- 14. **PASS** — Safe input scale for WiNN: All values are finite and nonnegative; no inverse transformation is required.
