# Derivative inputs for the revised validation (WiNN 0.1.6, 11 November 2026)

These files have the same layout as `derivative/`, the 27 July 2026 inputs. They are exported from the WiNN 0.1.6 derivative run, and current baseline, retrieval and effect comparisons use them.

Only the WiNN rows changed. Raw and competitor rows come from the July derivative run, reused on the same frozen matrices and seeds. Every Raw and competitor row the 0.1.6 run recomputed reproduced July exactly before reuse.

The public-dataset inputs are:

- BatchCorr per-identity matching;
- MTBLS79 and WaveICA per-feature effects and leave-one-batch-out coefficients;
- Sacurine covariate associations and leave-one-batch-out coefficients.

The only clinical file contains aggregate bootstrap draws (bootstrap, method, metric, value). It has no clinical identity labels, feature identifiers or matrices. SHA-256 checksums in `input_manifest.csv` refer to the exact uncompressed bytes.

Use this directory as `--derivative-root` for `manuscript/scripts/nearest_baseline/evaluate_baseline_biology.R`. The public-only reproduction command uses `--datasets=batchcorr_set1,mtbls79,waveica,sacurine`. Add `--reuse-baselines-from=<previous output directory>` to reuse the WiNN-independent baseline rows instead of re-evaluating them.
