# Frozen partial-confounding full-grid configuration

This directory was generated before any full-grid correction run.
It fixes the 320 scenario order, 18 analysis labels, 37 metrics,
hidden-reference seeds, and the material-degradation rule.
`portable_matrix_hashes.csv` freezes v3 SHA-256 identities after rounding
the complete pre-exponentiation raw/clean log matrices to 9 decimal
places. Exact bundle file hashes remain mandatory. Intensity object and
intensity-round9 hashes are platform diagnostics, not identity gates.
Exact reference-platform raw/clean log-object hashes are also frozen as
diagnostics; only the quantized log hashes are cross-platform gates.
The equivalence tolerance is one 1e-9 quantization cell; exact dimensions
and dimnames are also included in the serialized digest.
Every scenario has a new portable_execution_config_sha256_v3; v1 runs
and v2 runs must be archived and recomputed and are never silently reused.

Nominal strength 1.0 means maximal feasible confounding under the
balanced five-plate design, not complete global confounding. The primary
realized run-order axis is the weighted mean absolute within-plate
phenotype/order correlation; pooled correlation is diagnostic only.

Phenotype-aware ComBat is unavailable through the unchanged standard
wrapper and is not included in the phenotype-blind ranking.
