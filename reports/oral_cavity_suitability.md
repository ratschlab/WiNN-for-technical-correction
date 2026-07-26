# Suitability audit: spatially and temporally resolved human oral-cavity metabolome

**Verdict: UNSUITABLE — STOP BEFORE CORRECTION BENCHMARKING.**

The blocking reason is the absence of a deposited feature-by-injection intensity matrix. No correction method was run, no raw spectrum was downloaded, and no metadata was fabricated.

## Study and audited sources

- Article: Ciurli et al., *iScience* 27 (2024), 108884, DOI [10.1016/j.isci.2024.108884](https://doi.org/10.1016/j.isci.2024.108884), PMC10839270.
- Positive ion, batch 1: Mendeley Data [10.17632/fjkwnhmdjp.1](https://doi.org/10.17632/fjkwnhmdjp.1).
- Positive ion, batch 2: Mendeley Data [10.17632/3ybszwwfww.1](https://doi.org/10.17632/3ybszwwfww.1).
- Negative ion, batch 1: Mendeley Data [10.17632/tnbksjfv36.1](https://doi.org/10.17632/tnbksjfv36.1).
- Negative ion, batch 2: Mendeley Data [10.17632/9v5bw5zr9j.1](https://doi.org/10.17632/9v5bw5zr9j.1).
- Europe PMC full-text XML and supplementary package.

The four Mendeley file inventories contain exactly the same 18,096-byte `MSDIAL setting.docx` plus one raw ZIP per record:

| Assay | Raw archive | Deposited size | Downloaded? |
|---|---:|---:|---:|
| positive, batch 1 | `pos1.zip` | 4,253,676,252 bytes | no |
| positive, batch 2 | `pos2.zip` | 8,609,811,300 bytes | no |
| negative, batch 1 | `neg1.zip` | 8,178,311,256 bytes | no |
| negative, batch 2 | `neg2.zip` | 6,651,379,595 bytes | no |

Together these raw archives are 27.69 GB. They were inventoried through the public Mendeley API and were not downloaded.

## What the small files contain

The MS-DIAL document contains peak detection, alignment, gap-filling, and identification settings only. It contains no intensities or injection records.

The article supplements are downstream summaries:

| File | Dimensions | Contents |
|---|---:|---|
| `mmc2.xlsx` | 179 × 17 | metabolite names, identification level, group-comparison p-values/fold changes, PLS-DA VIP |
| `mmc3.xlsx` | 179 × 10 | AT-stratified p-values and fold changes |
| `mmc4.xlsx` | 179 × 11 | BT-stratified p-values and fold changes |
| `mmc5.xlsx` | 179 × 10 | CK-stratified p-values and fold changes |
| `mmc6.csv` | 405 × 26 | metabolite-library plate positions, structures, and classifications |

None has one row or column per analytical injection. In particular, no supplement contains the 159 biological samples, pooled-QC injections, or blanks as intensity profiles.

## Metadata that is documented but not independently tabulated

Mendeley documents the raw filename structure as:

`[polarity]_[batch]_[sample type]_[subject id]_[time]_[oral location]_[injection number]`

This is good evidence that batch, sample type, biological identity, and acquisition order existed in the raw filenames. The article also reports 159 saliva samples from 20 subjects, three oral locations (AT, BT, CK), three collection times (M, A, E), randomized injections, and a pooled QC after every five samples. However, the deposit has no small standalone injection manifest. The actual filename list is inside the raw ZIPs, so per-batch counts, precise QC identities, and exact sequences were not reconstructed under the prohibition on downloading raw spectra.

## Criterion results

| Requirement | Result | Reason |
|---|---|---|
| Uncorrected/minimally processed intensity matrix | **FAIL** | no processed matrix is deposited |
| Reliable batch and acquisition order | **FAIL** | schema is documented, but no row-level non-spectral manifest is available |
| Identifiable QCs/reference samples | **PASS WITH LIMITATION** | QC type and intended spacing are documented; individual identities are not available outside raw archives |
| Biological labels/replicate identities | **PASS WITH LIMITATION** | subject/time/location are encoded by design, but no row-level non-spectral manifest is available |
| Sufficient samples/QCs per batch | **FAIL** | intended QC cadence is reported, but usable counts cannot be verified |
| No prior batch/drift correction | **FAIL** | there is no matrix whose processing state can be assessed |
| Study samples and QCs in one matrix | **FAIL** | no matrix is deposited |
| Unambiguous matrix-to-metadata match | **FAIL** | neither matrix columns nor a standalone injection table is deposited |

## Smallest additional resource that would make the study auditable

The minimum resource is an original-scale MS-DIAL alignment export (`.txt` or `.csv`) containing peak areas or peak heights for every biological sample, pooled QC, and blank, with columns named by the deposited raw filenames. If those column names are not preserved, a separate injection manifest is also required with polarity, batch, injection number, sample type, subject ID, collection time, and oral location. The export must be prior to run-order correction, batch correction, QC normalization, and global scaling.

With that export, the existing raw archives would not need to be downloaded or reprocessed. Without it, the requested benchmark cannot be run responsibly.
