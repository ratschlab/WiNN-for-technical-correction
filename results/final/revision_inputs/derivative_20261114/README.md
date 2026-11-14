# Derivative inputs: permutation-calibrated WiNN

These public aggregate and public-dataset inputs come from the 14 November 2026 rerun using frozen WiNN development version 0.1.7.9999. Package provenance is in `package/frozen_package.json`; execution and validation details are in `analysis/config/calibration_revision/`.

The three WiNN modes were recomputed from fresh corrected matrices using the same 2,000 bootstrap draws, 5,000 label permutations where applicable, endpoints and seeds. Raw and comparator rows retain their exact historical values after numerical verification. Public-dataset files contain BatchCorr identity matching, MTBLS79/WaveICA feature effects, Sacurine covariate effects and leave-one-batch-out coefficients. The clinical files contain aggregate bootstrap draws and method-level point estimates only, with no identities, feature identifiers or matrices.

Use this directory as `--derivative-root` for `manuscript/scripts/nearest_baseline/evaluate_baseline_biology.R`. Public-only reproduction uses `--datasets=batchcorr_set1,mtbls79,waveica,sacurine`; `--reuse-baselines-from` retains the unchanged baseline rows. The five `biological_point_metrics.csv` exports provide the independent point-estimate checks used by the revised evaluator. `input_manifest.csv` records exact uncompressed bytes and SHA-256 hashes.
