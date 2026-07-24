# Sacurine negative-ion suitability audit

**Verdict: SUITABLE WITH LIMITATIONS.** The W4M00001 files provide an aligned, pre-correction positive-intensity matrix with reliable QCs, two analytical batches, within-batch acquisition order, and complete study covariates. The benchmark may proceed.

## Source and provenance

Pinned source: `https://github.com/SciDoPhenIA/phenomis` at commit `1e83ce6997a8d16b89ce5f0f899a1570004ebc0e`.
The three deposited W4M tables are used directly. No raw spectra, mzML, vendor files, ropls corrected object, or post-correction W4M output was downloaded.

## Design and scale

- Matrix: 113 retained identified metabolites × 210 injections; 184 study samples and 26 pooled QCs.
- Batch ne1: 90 study + 14 QC; batch ne2: 94 study + 12 QC.
- `injectionOrder` resets across batches but is unique and strictly increasing within each batch. The original values are preserved; within-batch ranks and a ne1→ne2 concatenated global order are derived explicitly.
- Intensities range from 506.5 to 678587872.0 (median 673549.5), with no NA, nonfinite, negative, or zero values.
- Repeated positive lower-bound values are retained; the source does not designate them as missing. No filtering or imputation was required.

## Evidence that the table is pre-correction

The pinned phenomis vignette first reads these exact three files and subsequently calls `correcting()` for QC-loess drift and batch adjustment, then performs pool-CV filtering, removes pools, divides by osmolality, and applies log10. This ordering is direct evidence that the deposited table precedes those operations.

## Balance and limitations

Gender versus batch: chi-squared p=0.1094, Cramer's V=0.118. Age and BMI batch standardized mean differences are -0.017 and 0.016. Within-batch correlations and the complete sequence are saved under `data/public/sacurine/audit/`.
The panel contains only 113 identified metabolites, so conclusions do not generalize to all detected peaks. Each volunteer occurs once, so no biological-replicate correlation or participant ICC will be calculated. Technical metrics will use pre-osmolality values; a common post-correction osmolality normalization will be used only for biological analyses.

## Criterion ledger

- 1. **PASS** — Numerical feature-by-injection matrix: 113 numeric metabolite rows by 210 aligned injections.
- 2. **PASS WITH LIMITATION** — Uncorrected or minimally processed intensities: Positive identified-metabolite peak intensities; the table is a selected 113-metabolite panel rather than the full untargeted peak table.
- 3. **PASS** — No prior batch/drift/QC correction: The phenomis vignette reads these files and only then invokes correcting(), QC-CV filtering, osmolality normalization, and log10 transformation.
- 4. **PASS** — Study and pooled QC in same matrix: 184 study samples and 26 pool QCs coexist in the same table.
- 5. **PASS** — Unambiguous matrix-metadata matching: Matrix column names exactly equal metadata row names in deposited order.
- 6. **PASS** — Reliable supplied batch labels: Two deposited batches: ne1 (104 injections) and ne2 (106 injections).
- 7. **PASS** — Reliable acquisition order: Deposited injectionOrder is unique and strictly increasing within each batch; it resets across batches, so global order is a documented concatenation.
- 8. **PASS** — Identifiable QCs: sampleType=pool identifies 14 ne1 and 12 ne2 QCs.
- 9. **PASS** — Usable biological labels: Age, BMI, gender, sampling, and osmolality are present for every study sample.
- 10. **PASS** — Sufficient study samples per batch: 90 and 94 study injections in ne1 and ne2.
- 11. **PASS** — Sufficient training QCs per batch: Four hidden QCs per batch leave 10 and 8 training QCs.
- 12. **PASS** — Boundary QCs can remain in training: First/last deposited QCs in both batches are retained by the holdout constraints.
- 13. **PASS** — No severe biological-batch/order confounding: Gender-batch Cramer's V=0.118; age SMD=-0.017; BMI SMD=0.016; both genders occur throughout both batches.
- 14. **PASS** — Safe input scale for WiNN: All values are finite and nonnegative; no undocumented inverse transformation is required.
