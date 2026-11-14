# Nearest study-sample regression baseline

`doBC.R` is copied unchanged from the authors' BatchCorrMetabolomics 0.1.14 implementation, commit `e0c7668140e206dcdae9afa602dd2e1b337ac4f6`:
https://github.com/rwehrens/BatchCorrMetabolomics/blob/e0c7668140e206dcdae9afa602dd2e1b337ac4f6/R/doBC.R

Upstream authors: Ron Wehrens, Roland Mumm, Joost Keurentjes and Ric de Vos. Upstream license GPL (>=2); GPLv3 included in COPYING. The project describes its implementation as illustrative. The revision evaluates these baselines on hidden-reference precision, exact simulation truth, repeated-identity discrimination, and empirical feature effects.

Run:

```sh
Rscript run_nearest_baselines.R --input-root=/path/to/prepared_inputs --output-root=/path/to/new_outputs --datasets=mtbls79,batchcorr_set1,sacurine,waveica
```

Each authorized input file is named `input_<dataset>.rds` and is a list with numeric, nonnegative feature-by-sample matrix `x`; sample-aligned `meta` containing `sample_id`, `batch`, `run_order`, `is_study`; and character vectors `training_ids`, `hidden_ids`, `external_ids`. Metadata may additionally contain `class` (MTBLS79) or `biological_group` (WaveICA). Clinical inputs are private and intentionally excluded.

All-injection ANCOVA fits `X ~ S * B` using `method="lm"`, supplied run order, batches and every intensity; no QC or biological labels enter fitting. Mean-only correction uses the same implementation with `seq.idx=NULL`. The study-only sensitivity excludes reference injections and thus receives additional reference/study identity information. All fits operate on log1p intensities; expm1 outputs are floored at zero, with affected cell counts recorded. There is no parameter selection or use of held-out endpoints during fitting.

Outputs include full corrected matrices and descriptive held-out reference CV, categorical batch weighted-PC R2 among study samples, and group weighted-PC R2 where available. The latter is a relative structural endpoint. The additional evaluation below supplies paired uncertainty for simulation truth and repeated-identity discrimination; effect comparisons against raw input remain descriptive. Four public matrices used in the revision match frozen-release SHA256 after preparation under R 4.5; the clinical reconstruction differs at floating-point precision and is explicitly not byte-identical to the frozen matrix.


## Complete validation for the second manuscript review

Run from the manuscript root after preparing the original benchmark inputs and baseline matrices:

```sh
Rscript scripts/nearest_baseline/run_simulation_baselines.R --input-root=/path/to/prepared_inputs --cores=4
Rscript scripts/nearest_baseline/evaluate_baseline_biology.R --input-root=/path/to/prepared_inputs --derivative-root=/path/to/analysis_repository/results/final/derivative
Rscript scripts/nearest_baseline/generate_validation_tables.R
```

`run_simulation_baselines.R` fits the two all-injection baselines on the exact ten canonical raw matrices, verifies all 20 raw/clean object SHA256 hashes against `scripts/benchmark_revision/simulation_bundle_hashes.csv`, and reproduces the archived raw truth metrics before adding frozen WiNN per-seed endpoints. Pass `--reference=/path/to/simulation_seed_method_metrics.csv` to compare against a newly validated WiNN release; the permutation-calibration revision uses this option with the complete new aggregate. Paired differences use the same ten seeds and t intervals with 9 degrees of freedom. Outputs are `generated/review_validation/baseline_simulation_*.csv`; Table S27 is generated from these aggregates.

To reproduce only public datasets, add `--datasets=batchcorr_set1,mtbls79,waveica,sacurine` to the empirical command. `--output-root=/path/to/new_aggregate_outputs` can keep a partial rerun separate from the complete manuscript aggregates. The evaluator accepts uncompressed `.csv` or corresponding `.csv.gz` frozen derivative files.

`evaluate_baseline_biology.R` evaluates all-injection baseline matrices already produced by `run_nearest_baselines.R`. It does not refit WiNN or other published methods. Biological-effect evaluation exactly reproduces the archived raw per-feature effects before comparison. WiNN effect coefficients and leave-one-batch-out (LOBO) results are read directly from the frozen derivative export: local caches are not substituted. MTBLS79 biological models first average study profiles within each of 20 independent extracts; the grouping is repeated within each training/held-out fold. WaveICA models use study-profile groups. Sacurine divides intensity by osmolality before log1p and fits age, BMI, sex and sampling jointly, reporting the three biological coefficients separately. LOBO uses jointly corrected matrices and measures coefficient consistency across batches, not prospective correction transfer. Table S28 contains these descriptive matched-feature comparisons; matching raw effects more closely does not establish matching biological truth.

### Retrieval estimand and corrected uncertainty

The retrieval query is a study injection (BatchCorr) or pair-collapsed correction unit (clinical FIA-MS). Similarity is Pearson correlation across the common features on log1p intensity. Candidates are all finite-similarity study profiles in other batches. A query is eligible only if a same-identity candidate exists. Rankings use stable decreasing correlation order; exact ties retain the original input-column order. The top-1 result is whether the first candidate has the same identity. Reciprocal rank uses the first same-identity candidate. Both point estimates average eligible queries, not identities.

There are 693 eligible queries from 346 identities among 710 BatchCorr study injections (357 total identities); there are 536 eligible clinical correction-unit queries from all 267 identities. New baseline and raw matrices have identical eligible-query sets and no exact top-score ties. BatchCorr candidate counts range 635--667 (median 636); clinical counts range 456--482 (median 458). Historical-method exact tie counts cannot be recovered from archived summary tables and are not fabricated; their stable sorting rule is unchanged.

The original public derivative script averaged per-identity retrieval means in its bootstrap, while its point estimates averaged queries. The corrected S23 resamples the same biological identities in the same 2,000 paired draws, retaining each resampled identity's number of eligible queries in the numerator and denominator. This aligns all intervals with the published query-weighted estimand. Correction fits, candidate sets and feature panels remain fixed. Original bootstrap seeds are `2026072701 + 199` (BatchCorr) and `2026072701 + 1594` (clinical). Original aggregates are retained with `_before_query_weighting` suffixes for audit; they are not current estimates.

The separate BatchCorr identity-balanced sensitivity export, `batchcorr_identity_retrieval_paired.csv`, retains its original seed `20261027` and 5,000 paired draws. Reproduce it from the complete derivative export without refitting:

```sh
Rscript scripts/nearest_baseline/refresh_identity_sensitivity.R <derivative-root> <output-root>
```

BatchCorr per-identity means are read from the archived derivative table. The clinical per-identity table was omitted from the original public release and the original private matrix export is no longer present. Historical comparator per-identity means are recovered from the complete archived bootstrap replicate series using the exact known resampling design. The 2,000-by-267 system has full column rank, giving a unique solution. Each method/endpoint must reproduce all 2,000 archived replicates to absolute error below 1e-11, lie within metric bounds, and reproduce its original query-weighted point. Direct raw identity summaries provide an additional independent check (agreement below 1e-10). Top-1 means are constrained to their known integer-success/query-count grid after floating-point inversion. Current WiNN per-identity summaries are read directly from the new derivative run; reconstruction is not used for these refreshed outputs. Reconstructed identity rows are never exported. Only aggregate checks, intervals, paired differences and bootstrap estimates are written.

`derivative_biology_helpers.R` is copied unchanged from the frozen local derivative analysis export (SHA256 `d687a7335bdfbdc6cc9339c82c5261cbbc7344ef857c98c28b649a2fd44a0f09`). The export records analysis base commit `21de7c448bfa5545d27c5a930161fb0bd8e40d85`; the helper itself was added after that base and is included explicitly in this revision release. The public inputs remain reproducible from the registered preparation pipeline; clinical reproduction requires authorized access to the forthcoming cohort data. The clinical aggregate outputs contain no sample identifiers or individual intensity profiles.

The original evaluation passed 132 empirical validation checks. The current permutation-calibration refresh passes 91 checks against the freshly recomputed WiNN derivative points and the unchanged baseline rows and resampling draws. Results support WiNN's improvement over both baselines on all three simulation truth metrics in every seed. Empirical retrieval differences against the baselines have intervals spanning zero, and several baseline effect slopes are closer to the raw effects. These results do not support universal empirical superiority over the simpler baselines.
