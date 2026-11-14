# Reproducibility and validation

## Frozen software and inputs

All final WiNN runs use the bundled 0.1.6 source archive. Each production task checks the installed package version, installation path, and source-archive SHA-256 before loading a dataset. Run manifests retain the package record, input and output hashes, seed, runtime, warnings, and session information.

WiNN 0.1.6 replaced 0.1.4 because the earlier versions silently fell back from the documented spline smoothers to LOESS or linear fits. Only WiNN outputs depend on that change, so the 0.1.6 release re-ran the three WiNN configurations, with Raw as a reproduction check, and reused every comparator output from the validated 0.1.4 release (`analysis/euler/submit_winn_rerun.sh`, `analysis/splice_winn_rerun.R`). The re-run used the same input matrices, seeds, reference partitions and resampling indices, so comparator and WiNN results remain paired. The Raw rows reproduced the previous release to within 5e-14, and `analysis/combine_winn_rerun_validation.R` combines the validation of the re-run with that of the reused comparator outputs.

Public downloads are pinned to a source commit or file checksum. Download scripts reuse a file only after validation and record its URL, size, checksum, and download date.

## Training and held-out references

Reference assignments are fixed before method fitting. Held-out references are relabelled as ordinary observations in method-facing metadata and are excluded from fitting, automatic selection, and candidate ranking. Competitor candidates are selected only from published/native defaults, method-native procedures, or training-QC behavior. The frozen result is evaluated once on held-out QCs and preservation endpoints.

The primary public-dataset holdouts are stored under `config/holdouts/`; the simulation holdout is stored with the canonical sample metadata. Both the holdout files and the larger task/split manifests are covered by `analysis/config/file_manifest.csv`.

The selection manifest records the candidate count, selected setting, selection seed, source package and version, information available during selection, and the justification for any departure from a native default. It also records that every configuration was frozen before final evaluation. Automatic-batch WiNN receives training-QC identities and run order but no supplied batch labels.

Reference-split stability repeats this separation across 10 shared splits per dataset. Split results describe within-dataset sensitivity; they are not treated as independent studies.

## Matrix invariants

Every task validates numeric type, finite nonnegative input, unique feature and sample identifiers, exact metadata alignment, batch and run order, output dimensions, and feature/sample retention. A method cannot silently use a smaller evaluation panel. Failures are written explicitly and successful tasks are marked only after all required outputs have been installed atomically.

All methods are evaluated on the same nonnegative intensity domain. Some regression-based competitors can extrapolate below zero even when their inputs are nonnegative. As in the original benchmark workflow, those negative corrected values are floored at zero before training-QC scoring or final evaluation. The number of affected values, features, and samples and the pre-floor minimum are retained with the method details. Non-finite canonical inputs and non-finite outputs outside the documented TIGER safeguard remain hard failures.

TIGER has one narrower method-specific safeguard. If a training-QC split makes TIGER return a non-finite prediction for a feature, the complete feature is restored from the aligned uncorrected input. This retains the prespecified feature panel without combining raw and corrected cells within a feature. The native non-finite count and the numbers of affected features and samples are saved for every TIGER run; all other finite TIGER output remains unchanged.

The released Euler analysis reused a correction artifact only when its input, output, dimensions, seed, package, and configuration were recoverable and checksum-validated. Because those large matrices are not committed, a clean public reproduction refits the corresponding method. Both decisions are recorded in the reuse audit; matrix loading time is never reported as fitting runtime.

The internal cohort is processed on independent adjacent-pair units. Its technical diagnostics never use the original alternating duplicate-injection sequence.

## Metrics

The common technical endpoints are held-out QC CV, residual run-order GAM deviance, residual Ljung-Box structure, categorical batch weighted-PC R-squared, runtime, and output coverage. Truth recovery is used only for the simulation. Biological and replicate endpoints follow each dataset's actual design; no pseudo-replicates are created. Metabolite ICC(A,1) is reported only where matched cross-batch feature ratings support that definition. Clinical FIA-MS and BatchCorr Set 1 instead report their explicitly named between/within feature repeatability ratio, alongside genuine replicate-profile agreement.

Batch is always passed to weighted-PC R-squared as a categorical factor. Automated validation fails if the categorical metric note is absent.

## WiNN ablation

The cumulative sequence is Raw, outlier shrinkage, selective drift correction, selective residual-batch correction, and PQN. The gate experiment compares selective/selective correction with all/selective, selective/all, and all/all correction while keeping other settings fixed. Step-impact tables use adjacent cumulative stages. Positive impact values always favor adding the stated stage after accounting for the metric's direction.

## Robustness analyses

The repeated simulation contains 10 prespecified seeds. `SIM01` is distributed because it is the common input for the primary comparison and main ablation; `SIM02` through `SIM10` are reconstructed from the component-level seed ledger. Reconstructed matrices are checked with the frozen identities in `analysis/config/simulation_bundle_hashes.csv`. The cross-platform gate uses a nine-decimal quantized matrix hash, while exact reference-platform object hashes are retained as provenance diagnostics. This accommodates harmless last-bit differences between R platforms without allowing changes in dimensions, names, or numerical content at the declared precision.

The reference analysis contains 10 splits for each dataset and method. The original partial-confounding grid contains 10 seeds and 16 conditions per seed, evaluating Raw and fixed/default WiNN; it remains in Supplementary Figure S5 and Table S7. The revised Figure 2E–F design uses the same technical artifacts and WiNN 0.1.6, 19 conditions per seed, fixed 180:270 group counts, and a complete batch-confounding boundary. Its primary endpoint is the retained marginal group contrast; Welch/BH discovery sensitivity accompanies FDR. The former adjusted coefficient is retained as a sensitivity analysis and marked non-estimable at complete batch separation. See `analysis/figure2_confounding/PROTOCOL.md` and `README.md` in that directory for the fixed design and reproduction commands.

## Figures and aggregate tables

Final tables are assembled only from completed dataset-level outputs. Aggregation checks declared task counts before writing. Every rendered figure has a source CSV and shares the method and dataset palettes used throughout the repository. `analysis/validate_release.R` audits task completeness, package provenance, unresolved failures, categorical batch coding, output coverage, and figure-source links.

## Privacy boundary

The restricted clinical input belongs under `data/private/clinical_fiams/`, which is ignored by Git. Only aggregate clinical metrics enter `results/final/`. Repository validation scans tracked paths and content before release.
