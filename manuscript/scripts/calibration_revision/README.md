# Standalone selection-test comparison (Supplementary Methods S10 / Figure S11)

These files preserve the completed paired calibration experiment. The permutation
rule is an experimental alternative, **not the WiNN 0.1.6 implementation used in
Figure 2 or any other correction benchmark**. Do not install this source snapshot
or substitute its outputs into the primary benchmarks. The snapshot is sourced
only by the diagnostic runner, which compares the two rules on identical data.

The original null experiment is unchanged; all 8,000 original selection counts
must reproduce exactly. The paired experiment and raw results were recovered from
manuscript revision `6e6bff2`, without refitting or selecting conditions. Its
prespecified design has 11,600 datasets and 31,200 stage-by-algorithm rows.
The main manuscript's former Section 3.5 is now in Supplementary Methods S10.

From the manuscript root, reproduce the experiment (or split it into independent
shards using the optional final two arguments):

```sh
Rscript scripts/calibration_revision/run_gate_validation.R \
  scripts/calibration_revision/source_snapshot/winn /path/to/results 1 40
# Repeat shards 2 through 40, then aggregate all shards:
Rscript scripts/calibration_revision/summarize_gate_validation.R /path/to/results \
  generated/review_validation/null_gate_verification_replicates.csv
Rscript scripts/calibration_revision/plot_gate_validation.R \
  generated/calibration_revision/gate_validation_summary.csv figures/supplementary_panels
```

`validate_saved.R` independently regenerates all summary rows from the archived
replicates and checks parity with the original null experiment. It requires no
new correction fits. `source_snapshot/` retains the experimental version label
for provenance; production WiNN remains 0.1.6. The comparison changes the error
criterion from nominal BH FDR to per-segment maxT familywise error. Under the
global null these criteria coincide; under mixed signals they do not. The power
comparison and failures under scheduled sample-type differences must accompany
the improved null-selection rates.
