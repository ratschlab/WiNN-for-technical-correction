# Frozen derivative inputs used by the revised validation

These files are the exact, compressed inputs from the locally archived 27 July 2026 derivative analysis export. They were absent from public commit `21de7c4`; this revision publishes the required sources. They are historical inputs for reproducing and correcting the original uncertainty, not current publication estimates. Current query-weighted estimates are under `manuscript/generated/review_validation`.

Public-dataset inputs include BatchCorr per-identity matching and MTBLS79/WaveICA/Sacurine per-feature effects and leave-one-batch-out coefficients. The only clinical file contains aggregate bootstrap draws with method/metric/replicate/value columns; no clinical identity labels, feature identifiers or matrices are included. SHA-256 checksums refer to exact uncompressed source bytes.

Use this directory as `--derivative-root` for `manuscript/scripts/nearest_baseline/evaluate_baseline_biology.R`. That reader accepts `.csv.gz` files directly. The public-only reproduction command uses `--datasets=batchcorr_set1,mtbls79,waveica,sacurine`.
