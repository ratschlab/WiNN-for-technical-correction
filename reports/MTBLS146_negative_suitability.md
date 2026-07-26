# MTBLS146 negative-ion benchmark suitability audit

**Verdict: FAIL — MTBLS146 is not suitable for the correction benchmark.**

Audit date: 2026-07-22

The correction benchmark was not run. Complete remote inventories at MetaboLights and GigaDB contain converted LC-MS spectra and ISA-Tab metadata, but no numerical feature-by-injection intensity matrix. Per the benchmark's critical stop rule, the spectral collection was not downloaded and no peak picking, preprocessing, imputation, correction, tuning, or evaluation was attempted.

## Sources inspected

- [MetaboLights study API](https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS146)
- [MetaboLights public study directory](https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS146/)
- [MetaboLights complete spectral-file directory](https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS146/FILES/)
- [MetaboLights SHA-256 records](https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS146/HASHES/)
- [GigaDB DOI record](https://doi.org/10.5524/100108)
- GigaDB public object-store listing for prefix `live/pub/10.5524/100001_101000/100108/`
- [Published graphical run-order supplement](https://static-content.springer.com/esm/art%3A10.1186%2Fs13742-015-0054-9/MediaObjects/13742_2015_54_MOESM1_ESM.pdf)

Only transient copies of small metadata/listing files were inspected. No MTBLS146 files remain under `data/` in this repository.

## Repository inventory and matrix search

The current MetaboLights top-level directory contains the investigation file, study file, three assay files, archive metadata, and checksums. Its `FILES/` directory contains 634 spectra:

- 227 negative-ion untargeted LC-MS mzXML files;
- 216 positive-ion untargeted LC-MS mzXML files;
- 191 targeted-positive lipidomics/carnitine mzML files.

The complete GigaDB public object listing reports `KeyCount=636` and `IsTruncated=false`. It contains:

- 634 mzML spectral files;
- `ISA_pregnancy.zip`;
- `readme_100108.txt`.

No CSV, TSV, Excel, R object, PeakML, featureXML, or other feature-abundance table occurs in either complete inventory. The GigaDB API itself returns at most 50 file records, so the conclusion is based on the untruncated object-store listing and repository readme, not the truncated API response.

The article and ISA protocol state that XCMS/CAMERA was used to generate an RT-m/z-by-sample ion-intensity matrix. That matrix and its exact processing parameters were not deposited. In the negative assay, the `Derived Spectral Data File` field points to mzXML spectra; these are converted spectral files, not feature-intensity tables.

## Negative-ion assay identification and recoverable design

The ISA investigation identifies `a_MTBLS146_pregnancy_negative.txt` as metabolite profiling by mass spectrometry on an LTQ Orbitrap Velos. All 227 assay rows have scan polarity `negative`. Identification therefore comes from ISA-Tab, not filename guessing.

The assay metadata itself is strong:

- 227 total negative-ion injections;
- 161 study samples and 66 pooled-QC injections;
- five analytical batches with observed sizes 42, 43, 35, 44, and 63;
- per-batch QC counts 13, 14, 15, 13, and 11;
- gestational-class counts A=28, B=26, C=27, D=29, E=26, and F=25;
- unique `(batch, within-batch run order)` pairs for all 227 injections;
- first and last observed positions in every batch are QCs.

Run order comes directly from `Factor Value[Run order]` and batch from `Factor Value[Batch number]`; the published supplement corroborates the five sequential batches and QC layout. No alphabetical order is needed.

Exact gestational week is not recoverable. The study table contains interval labels. A-E correspond to 9-12, 13-16, 17-20, 21-24, and 25-28 weeks. Current deposited F rows say 29-32 weeks, whereas the publication and supplement define F as 29-40 weeks. A future analysis should retain A-F as the primary endpoint and set exact week to missing unless a more precise source is supplied.

## Randomization and balance

Class-by-batch counts among the 161 study samples are:

| Class | Batch 1 | Batch 2 | Batch 3 | Batch 4 | Batch 5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 4 | 5 | 5 | 5 | 9 |
| B | 7 | 7 | 3 | 3 | 6 |
| C | 4 | 2 | 4 | 8 | 9 |
| D | 3 | 5 | 4 | 7 | 10 |
| E | 4 | 6 | 2 | 6 | 8 |
| F | 7 | 4 | 2 | 2 | 10 |

Pearson's chi-squared statistic is 14.684 (df=20, p=0.794); Cramer's V is 0.151. Every class occurs in every batch, and the graphical acquisition sequence is visibly interspersed rather than class-sorted. There is no evidence of severe structural class-batch confounding. This does not rescue the dataset because the critical matrix criterion fails first.

## Suitability criteria

| ID | Criterion | Status | Evidence |
| ---: | --- | --- | --- |
| 1 | Numerical feature-by-injection abundance matrix is available | **FAIL** | Both complete inventories contain spectra plus ISA/readme files only. |
| 2 | Matrix represents uncorrected or minimally processed intensities | **FAIL** | No matrix exists to assess. |
| 3 | Matrix has not already undergone correction or destructive normalization | **FAIL** | No matrix exists to certify. |
| 4 | Study samples and pooled QCs are in the same matrix | **FAIL** | Study and QC spectra coexist, but no shared matrix was deposited. |
| 5 | Matrix columns match injection metadata unambiguously | **FAIL** | Spectra map to ISA injections, but there are no matrix columns. |
| 6 | Five batches and within-batch acquisition order are recoverable | **PASS** | ISA-Tab supplies complete, unique batch/order pairs and the supplement corroborates them. |
| 7 | QC identities and positions are recoverable | **PASS** | Sixty-six pooled QCs are explicitly identified and positioned. |
| 8 | Gestational classes are assignable | **PASS WITH LIMITATION** | A-F is reliable; exact week is unavailable and the deposited F interval differs from the publication. |
| 9 | Every batch has sufficient study samples for drift testing | **PASS** | Study counts are 29, 29, 20, 31, and 52. |
| 10 | Every batch has enough QCs for QC-based methods | **PASS** | QC counts are 13, 14, 15, 13, and 11. |
| 11 | First and last usable QCs can be retained for QC-RFSC | **PASS** | Every batch begins and ends with an observed QC. |
| 12 | Gestational class is not severely confounded with batch/order | **PASS** | All classes occur in all batches; Cramer's V=0.151; the published sequence is interspersed. |
| 13 | Data are not merely final normalized study values | **PASS WITH LIMITATION** | Deposited spectra are upstream of normalization, but using them requires a new raw-data-processing workflow. |

## Smallest additional resource needed

The smallest resource that could make MTBLS146 suitable is an **uncorrected XCMS/CAMERA feature-intensity matrix covering all 227 negative-ion injections**, with:

1. stable feature IDs and the original m/z and retention time;
2. one column per injection, named by `MS Assay Name` or deposited spectral filename;
3. peak areas before run-order correction, batch correction, QC normalization, or global scaling;
4. processing provenance, including XCMS/CAMERA parameters and zero/missing-value semantics.

The alternative would be a newly approved raw-data workflow over the 227 negative spectra, including peak picking, retention-time alignment, correspondence/grouping, gap filling, isotope/adduct annotation, and matrix-level validation. That would create a new processed dataset rather than recover the authors' matrix and was not authorized.

## Final validation

- Spectral dataset downloaded: **no**.
- Correction variants attempted: **0 of 9**, by design after the critical failure.
- Held-out QC labels exposed to methods: **not applicable; no method ran**.
- Biological-repeat metrics fabricated: **no**.
- MTBLS79 and simulation files overwritten: **no**; their pre-audit SHA-256 hashes were unchanged after inspection.
- Preprocessing outputs/notebook/render created: **no**, because doing so would violate the stop rule.
