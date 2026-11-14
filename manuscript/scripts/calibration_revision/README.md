# Drift-gate calibration revision

This revision retains the Ljung–Box statistic, adaptive lag, log1p scale,
MAD preparation and smoothing. The default reference becomes 999 common,
uniformly sampled permutations of the injection columns within each segment.
For each permutation the maximum statistic across testable features is retained.
The adjusted p-value is `(1 + number of permuted maxima >= observed Q)/(999+1)`.
There is no second BH adjustment. Ties are counted conservatively.

At a fixed segment and prespecified threshold alpha this controls the probability
of any false selection if the joint matrix of null features is invariant under
column permutations. This is stronger than a null of zero autocorrelation.
It allows dependence between metabolites but does not cover arbitrary temporal
structure, scheduled sample-type differences, data-selected segmentation,
QC-selected thresholds, or the complete preprocessing/correction pipeline.
The unchanged ANOVA stage continues to use BH. The legacy `fdr_threshold`
argument denotes the maxT threshold for the permutation drift gate.

## Prespecified checks

- Reuse all 4,000 original Gaussian-null datasets (100 features; 1,000 seeds at
  each segment length 20, 50, 100, 200), before and after MAD preparation.
- Compare both algorithms on exactly the same profiles; verify all 8,000 legacy
  stage rows against the previous archived experiment, not just rounded means.
- Evaluate 100 seeds for each combination of four segment lengths, shared
  component weights 0 and 0.6, Gaussian/t3/skewed noise, and global null,
  20%-linear-drift, or periodic-QC scenarios. Periodic QC differences violate the
  permutation null and are a limitation check, not a nominal-error experiment.
- Evaluate 1,000-feature global-null families at lengths 50 and 200, shared
  component weights 0 and 0.6, over 100 seeds.
- Independently check 1,000 Gaussian-null datasets at length 50 using the same
  production default permutation seed in every fit.

The main design contains 11,600 datasets and 31,200 algorithm/stage rows.
`run_gate_validation.R` supports disjoint deterministic shards. Merge only after
`summarize_gate_validation.R`'s completeness checks pass. Report exact binomial 95% Monte Carlo
intervals for the probability of any false selection, mean false-discovery
proportions and sensitivity. All permutations and simulation seeds are recorded.
No parameter is chosen by the calibration results or biological endpoints.

Existing source snapshots remain historical evidence. New benchmark results must
use the frozen updated package, identical inputs, seeds, split assignments and
endpoint code. Comparator outputs and MAD/PQN-only stress results are unaffected.

The package source needed by the validation is included under
`source_snapshot/winn/` (development version 0.1.7.9999, commit
`ab4969792c98bd86208a28a4e83b7a69e204b9fc`). The source manifest records file
hashes and the full package archive hash. To run one of 40 disjoint shards:

```
Rscript scripts/calibration_revision/run_gate_validation.R \
  scripts/calibration_revision/source_snapshot/winn tmp/gate_validation 1 40
```

After running shard indices 1 through 40, merge and verify:

```
Rscript scripts/calibration_revision/summarize_gate_validation.R \
  tmp/gate_validation generated/review_validation/null_gate_verification_replicates.csv
```

The complete merged results, design and Monte Carlo summaries are stored in
`generated/calibration_revision/`. Shard hashes document the 40-part merge.
`plot_gate_validation.R` creates Supplementary Figure S11 from the validated
summary. `verify_package_scope.R` confirms that 37 existing function bodies,
including all preprocessing, smoothing, batch-correction and PQN functions,
are unchanged; the five changed existing bodies are the drift gate and its
call paths. The old and new snapshots make that comparison self-contained.

Metric interpretation: the replicate-level false-discovery fields compare
selection with the generating **instrumental-drift** indicator. Periodic QC
shifts violate the exchangeability null even though that indicator is false.
The summarizer therefore marks `null_model_applicable=FALSE` and leaves formal
`fwer`/`mean_fdp` fields undefined for that scenario. Its `any_selection_rate`
and `technical_false_selection_rate` describe failure to distinguish scheduled
sample-type structure from instrumental drift, not a failure of a valid null
permutation test. No scenario is removed based on its observed result.

## Additional checks and plots

From the manuscript root, rerun the independent fixed-seed experiment and
reproduce its interval and the deterministic scheduled-QC example:

```sh
Rscript scripts/calibration_revision/run_fixed_seed_check.R \
  scripts/calibration_revision/source_snapshot/winn \
  generated/calibration_revision/fixed_seed_null_replicates.csv
Rscript scripts/calibration_revision/summarize_additional_checks.R \
  scripts/calibration_revision/source_snapshot/winn generated/calibration_revision
Rscript scripts/calibration_revision/plot_gate_validation.R
```

The batch-gate check can be run in four disjoint shards (replace the shard index
`1` with each index from 1 through 4), then merged:

```sh
Rscript scripts/calibration_revision/run_batch_gate_check.R \
  scripts/calibration_revision/source_snapshot/winn tmp/batch_gate 1 4
Rscript scripts/calibration_revision/summarize_batch_gate.R tmp/batch_gate \
  generated/review_validation/batch_gate_calibration_replicates.csv \
  generated/calibration_revision
```

This check used local R 4.6.1 and mgcv 1.9.4, matching the historical batch-stress
environment. All 400 before-drift rows were unchanged. Running the same driver
against `scripts/targeted_validation/source_snapshot/winn` reproduced all 800
historical rows (four numeric metrics per row) exactly. Session records distinguish
these checks from the main Euler calibration experiments (R 4.5.1).
