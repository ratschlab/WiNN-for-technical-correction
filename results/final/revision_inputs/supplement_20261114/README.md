# Supplementary aggregate inputs: permutation-calibrated WiNN

These inputs come from the isolated 14 November 2026 rerun using WiNN development version 0.1.7.9999, package commit `ab4969792c98bd86208a28a4e83b7a69e204b9fc` and source archive SHA-256 `72c90896ac28f3ee4000a5c6d5287e7d4031e307c40b787d001755a6982c4c29`.

The rerun replaces the Ljung–Box drift gate's asymptotic/BH calibration with 999 joint-column permutations and maxT adjustment. The three WiNN modes, cumulative ablations, ten simulation seeds, repeated reference splits, partial confounding and downstream estimates were recomputed. Input matrices, candidate configurations, correction operators, endpoint definitions and resampling seeds were retained. The corrected auto-mode configuration count is 36 nominal combinations; the previous count of 72 in some metadata was a descriptive error, not a change to the search.

Raw and comparator rows retain their frozen values. As in the historical aggregate archive, QC-RFSC rows use the earlier adapter context; the supplementary table generator applies the unchanged corrected QC-RFSC exports from `manuscript/generated/review_validation`. Current publication inputs preserve those corrected comparator rows exactly.

The derivative tables use 2,000 bootstrap draws and, where applicable, 5,000 label permutations. New WiNN estimates are paired with the same saved comparator draws. The combined release validation records the new WiNN checks and the prior comparator validation. Detailed rerun code and provenance are in `analysis/config/calibration_revision/`.

Run `make tables` from the repository root. The table generator reads this archive, applies the corrected RFSC values and assembles the supplementary tables. No clinical identities, feature-level clinical data or matrices are included. `input_manifest.csv` records the exact uncompressed bytes and SHA-256 for every compressed CSV.
