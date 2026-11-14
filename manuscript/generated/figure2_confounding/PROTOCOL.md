# Figure 2E–F confounding design revision

Scope: change only the phenotype-allocation experiment and its evaluation for
Figure 2E–F. WiNN 0.1.6, technical-artifact generation, the ten existing CONF
master seeds, 1,000 features, 500 injections, five plates, QC placement, 100
responsive features and their balanced signed log effects of 0.35 remain fixed.
Figure 2A–D and every non-confounding benchmark remain unchanged.

The previous experiment retained an always-mixed middle plate and evaluated a
phenotype coefficient conditional on plate and acquisition position. That
coefficient can remain stable while correction removes the between-group
contrast. The revised experiment separates biological preservation from the
estimability of this conditional coefficient.

## Allocation, fixed before fitting

- Use 180 cases and 270 controls in every condition. This 40:60 ratio permits
  complete batch separation in the existing five equal plates without changing
  the technical matrices or discarding injections. At zero alignment each
  plate has 36 cases and 54 controls.
- Randomly choose two case-enriched plates once per seed, independently of
  feature values, drift and batch effects. Batch alignment interpolates plate
  case counts from (36,36,36,36,36) to a permutation of (90,90,0,0,0).
  Largest-remainder allocation preserves exactly 180 cases at every strength.
- Run-order alignment progressively exchanges late cases and early controls,
  retaining the assigned case count in each plate. At full alignment the two
  groups occupy nonoverlapping chronological blocks within each mixed plate.
- Joint alignment applies both operations. At strength one, batch-only and
  joint allocations necessarily coincide; this is an expected boundary, not
  an independently identifiable joint effect.
- Strengths: 0, 0.25, 0.50, 0.75, 0.90, 0.95, 1.00. Zero is represented once:
  19 conditions per seed, 190 fits total. Keep every planned condition.
- Center the biological log effect using the fixed prevalence 0.4. QC
  injections receive no phenotype effect. Correction sees no phenotype labels
  or response truth.

## Evaluation, fixed before fitting

Primary preservation: median signed case-minus-control mean log1p contrast
after correction divided by the same contrast in the artifact-free matrix,
over the 100 responsive features. Use the common log1p scale; do not divide by
the generating log-intensity effect. Also report contrast RMSE and sign
agreement. Neither group contrast is adjusted again for the variable with
which phenotype is deliberately confounded.

Discovery evaluation: two-sided Welch tests of case/control log1p values,
BH adjustment over all 1,000 features at 0.05, and known feature truth. Report
true-positive rate, false-discovery proportion (zero when there are no
discoveries), number of discoveries and the frequency of zero discoveries.
Figure 2F must show sensitivity alongside FDR: a low FDR with no discoveries
does not establish preservation. Analyze Raw, WiNN and artifact-free truth
using the same definitions. Truth is an evaluation control, never fit input.

Sensitivity analysis: retain the former phenotype+plate+within-plate-position
model, mark it non-estimable when rank deficient, and never replace a
non-estimable coefficient or FDP with zero. Report realized phenotype/batch
association, within-plate order association, overlap, design rank and condition
number. No algorithm, allocation, endpoint or seed is selected based on
whether WiNN deteriorates. This experiment describes this specified design;
it does not prove a general identifiability bound or infer true effects from
fully confounded observational data.

Means and Monte Carlo standard errors are across the same ten paired seeds.
The original adjusted experiment remains archived for comparison.
