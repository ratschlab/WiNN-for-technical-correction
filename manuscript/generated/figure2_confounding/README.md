# Figure 2E–F source data

The protocol was fixed before fitting. Only phenotype allocation and evaluation
changed; WiNN 0.1.6 and the original technical-artifact matrices stayed fixed.
There are ten paired seeds and 19 conditions per seed. See PROTOCOL.md for the
complete allocation rules, effect sizes, endpoints and boundary interpretation.

`seed_metrics.csv` contains every Raw, WiNN and Truth estimate, for marginal
(primary) and adjusted (sensitivity) endpoints. `summary.csv` contains means and
Monte Carlo standard errors. The shared zero-alignment condition is displayed
once on each curve; full batch and joint alignment are identical by design.
Adjusted coefficients at complete batch confounding are non-estimable, recorded
as missing. FDR is the mean seed-level false-discovery proportion, defined as
zero for zero discoveries; sensitivity and discovery counts must be read with it.

`design_diagnostics.csv`, `run_manifests.json` and `input_manifest.csv` document
realized allocation, package/input/protocol hashes and the compact task outputs.
`validation.json` checks completeness, fixed technical artifacts and boundaries.
`protected_baseline.json` and `scope_audit.json` verify that every pre-existing
result/figure outside the permitted Figure 2 outputs and the bundle containing the separately
requested supplementary S11 comparison remains byte-identical,
including Figure 2A–D and every selective/forced comparison.

Reproduction code and the compact synthetic task archive are in the analysis
repository under `analysis/figure2_confounding/` and
`results/final/revision_inputs/figure2_confounding_20261114/`. No clinical data
are involved. The original adjusted, balanced confounding experiment is retained
in Supplementary Figure S5 and Table S7.
