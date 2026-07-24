# Partial-confounding v2 archive and rerun protocol

Array 8385353 was held/cancelled after the v1 intensity-round10 gate
failed. The observed remote pre-v2 inventory was 33 canonical run
directories, 49 preflight directories, and 96 scheduler log files.
Before syncing v2, archive the entire runs/, validation/preflight/ and
validation/dry_run/ trees, the complete v1 config directory, and all
job-8385353 plus prior gate logs. Record archive paths and SHA-256 values.
Never delete or selectively overwrite failed/incomplete directories.

`pre_v2_attempt_ledger.csv` retains the task-level v1 attempt state and
the exact CONF03 expected/observed hashes and numerical comparison.
Every scenario, including the old task 1 result, requires a v2 run.

Release order:
1. Validate the regenerated config manifest and unchanged 320 bundle hashes.
2. On Euler, run CONF03_none_0p00 with --dry-run --force and retain proof
   that both authoritative log-round10 hashes pass while the raw-intensity
   round10 diagnostic mismatch remains non-blocking.
3. Run a full CONF01_none_0p00 v2 task with --force and all native checks.
4. Only after both gates pass, submit array tasks 2-320 with the v2 rerun
   flag so every pre-v2 directory is archived before recomputation.
5. Launch the strict v2 aggregate only with an afterok dependency on that
   array; v1 completion schemas are scientifically ineligible.
