# Reproduce Figure 2E–F

This directory changes only the confounding allocation and evaluation. The fixed
scientific protocol is in `PROTOCOL.md` (committed as `7d67937` before the first
fit). It does not modify WiNN or the canonical technical-artifact generator.
The original adjusted experiment remains available in Supplementary Figure S5
and Table S7.

From the analysis repository root, restore the recorded R environment, run
`make install-winn`, and set `R_LIBS_USER="$PWD/Rlib"`. Generate the original
simulation bundles with `make simulations`. No empirical data are needed.

```sh
Rscript analysis/figure2_confounding/test_design.R
# Run each task 1 through 190; outputs are independent. Example:
Rscript analysis/figure2_confounding/run.R "$PWD" results/figure2_confounding 1 \
  "$PWD/analysis/figure2_confounding"
# After all tasks finish:
python3 analysis/figure2_confounding/summarize.py \
  results/figure2_confounding manuscript/generated/figure2_confounding
Rscript manuscript/scripts/figure2_confounding/plot.R
python3 manuscript/scripts/figure2_confounding/assemble.py
```

`run.sbatch` records the original Euler array submission. `resume.sbatch` runs
the same tasks in eight concurrent processes within one allocation, skipping
completed tasks, to avoid per-task scheduling delays. Adapt only the execution
paths/module setup for another installation. Both require the frozen WiNN 0.1.6
archive and library; the runner rejects a different version or archive hash.

The compact task archive contains synthetic sample allocations, per-feature
marginal effects and tests, seed metrics, design diagnostics, and completion
records. The aggregate exports include both the marginal primary endpoint and
the adjusted sensitivity analysis, with non-estimable boundary results recorded
as missing. All 190 tasks are required before aggregation. The zero-alignment
condition is fitted once per seed and displayed at the start of each curve.
Full batch and joint alignment have identical allocations and results.

Corrected, raw and ground-truth matrices are saved during execution but omitted
from the public archive; they can be regenerated. No phenotype labels or feature
truth enter WiNN. Plotting writes only Figure 2E/F and their source CSVs; assembly
refreshes Figure 2 and the bundles containing it. Panels A–D are reused verbatim.

`validate_outputs.R` independently checks every saved matrix's sample/feature
alignment, recomputes BH adjustments and compares 2,280 feature tests against
base R's `t.test`. `validate.sbatch` records that post-fit Euler check.
The user-requested standalone calibration comparison is retained separately in
Supplementary Methods S10 / Figure S11; it does not alter these WiNN fits.
