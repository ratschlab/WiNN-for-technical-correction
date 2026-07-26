# Results

`final/` contains the compact, validated result release used to generate the benchmark tables and figures. It includes aggregate clinical statistics but no row-level clinical data. Large corrected matrices, candidate fits, per-feature diagnostics, scheduler logs, and caches are generated locally beneath `results/` and are not committed.

## Reading the comparison

The primary WiNN result is the supplied-batch fixed/default mode without QC identities. The two automatic modes are secondary analyses and are not selected independently for each endpoint. Competitor settings are published or native defaults, or are selected using training QCs only. Held-out QCs and biological or replicate information are reserved for evaluation.

No method is expected to lead every metric. Held-out QC CV, residual run-order structure, and batch weighted-PC R-squared measure different technical effects; truth, biological, and replicate metrics guard against improving a technical statistic by removing useful variation. Feature and sample retention should therefore be read alongside every apparent gain.

## Main interpretation

Across the truth-known simulation, fixed/default WiNN gives the strongest overall balance of technical-effect removal, truth recovery, and complete feature coverage. This pattern is stable across the ten prespecified simulation seeds. The QC-assisted supplied-batch mode is very similar, while automatic batch discovery is less reliable when the supplied batch design is already known.

The empirical datasets are deliberately heterogeneous. WiNN gives a particularly strong balance in the pair-aware clinical FIA-MS cohort and BatchCorr Set 1. In MTBLS79, several methods are competitive on individual endpoints; QC-RFSC has the lowest held-out QC CV but leaves more residual batch structure, whereas the selected QC-RLSC fit is numerically unstable. In Sacurine, QC-RFSC leads several technical metrics, while WiNN retains a favorable balance with biological variation on a small feature panel. WaveICA contains large batch and run-order effects. QC-RFSC and TIGER reduce held-out QC CV more strongly there, but leave about half of the weighted-PC batch variance; fixed/default WiNN gives the lowest residual GAM and Ljung-Box structure, near-zero batch variance, strong biological variance, and complete matrix coverage. These results support balanced, design-aware performance rather than universal superiority on every endpoint.

## What the ablations show

Selective drift correction is important. Forcing drift correction across all eligible profiles can produce a lower residual drift statistic, but it generally worsens held-out QC CV and truth or replicate preservation. This explains part of WiNN's advantage in settings where many metabolites do not require the same drift adjustment.

The residual batch stage produces the largest technical improvement in most datasets. Selective rather than forced batch correction has a smaller numerical effect here, although the gate remains a useful protection when technical and biological structure overlap. PQN contributes relatively little after the drift and residual-batch stages in these benchmarks.

The partial-confounding experiment identifies the expected boundary: WiNN remains effective under partial technical-biological overlap, but severe or complete run-order and joint confounding attenuate recoverable biological effects because the two sources of variation are no longer identifiable from the observed design.

## Files

The combined and per-dataset tables are under `final/tables/`. Publication figures are under `final/figures/`, and every figure has a CSV source table in `final/figures/source_data/`. Validation summaries, the logical task ledger, scheduler attempts, package provenance, and output-domain audits are under `final/validation/`.
