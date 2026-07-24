# Reproducibility and validation

## Frozen inputs

Public downloads are pinned to source commits or a file checksum. Download scripts reuse a local file only after checksum validation and write a machine-readable manifest. The synthetic analyses use explicit seed ledgers and content hashes.

## Hidden-QC protocol

Held-out controls are selected once and saved in `config/holdouts/`. Before correction, their QC labels are removed from method-facing metadata. They are excluded from parameter selection, drift fitting, method training, and WiNN automatic selection. The same held-out controls are used for every method within a dataset.

## Matrix invariants

Each runner checks numeric type, finite nonnegative values, feature and sample identifiers, metadata order, output dimensions, and feature/sample retention. A method-specific loss is recorded rather than silently imposing a smaller panel on other methods. Biological comparisons also include a common-feature panel where required.

## Metrics

The common technical endpoints are held-out QC CV, residual run-order GAM deviance, residual Ljung-Box structure, categorical batch weighted-PC R-squared, runtime, and output coverage. Biological endpoints follow each dataset's actual design. Genuine replicate metrics are calculated only when real repeated measurements exist.

The clinical benchmark uses adjacent-pair independent units for drift and batch modelling. No pseudo-replicates are introduced.

## Caching and reruns

Command-line runners are idempotent and reuse completed results unless `--force` is supplied. Robustness runs store source hashes, parameter hashes, matrix checksums, method status, errors, warnings, runtime, session information, and completion manifests. A forced rerun archives an incompatible completed result before replacement.

## Plotting

Method, dataset, gate, and stage colours are centralized in `scripts/plot_theme.R`. Figure source tables are saved beside the rendered outputs so plots can be regenerated without extracting values from graphics.

## Privacy boundary

This repository is safe to share publicly: it contains no internal clinical matrices, metadata, identifiers, paths, or cohort-specific adapter. The separate analysis-results package includes only aggregate clinical tables and publication figures after a privacy scan.
