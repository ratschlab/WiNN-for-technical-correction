# Partial-confounding v3 archive and rerun protocol

V2 gate job 8395516 completed task 1. Array 8396751 then produced
282 completed tasks, 17 raw-log round-10 portability-gate failures, and
20 tasks canceled when the defect was diagnosed. No correction method ran
inside the 17 failed tasks. A source-identical EPYC 9654 retry reproduced
the CONF09 failure.

An independent all-320 reconstruction comparison found 16 raw-log hash
mismatches at round-10, zero clean-log mismatches at round-10, and zero
raw or clean mismatches at round-9 and round-8. Round-9 is therefore the
narrowest tested portable decimal gate. Exact bundle hashes, dimensions,
dimnames, source-subset identity, and all scientific parameters remain
mandatory and unchanged.

Before syncing v3, move the complete v2 runs/, validation/, config/, and
scheduler logs into the checksum-manifested archive named in
`pre_v3_attempt_ledger.csv`. Never selectively reuse a v2 completion.

Release order:
1. Validate the v3 config manifest and unchanged 320 bundle hashes.
2. On an AMD compute node, run CONF09_none_0p00 with --dry-run --force
   and retain proof that both authoritative log-round9 hashes pass.
3. Run a full CONF01_none_0p00 v3 gate with --force and all native checks.
4. Submit tasks 2--320 with the v3 rerun flag only after both gates pass.
5. Launch the strict v3 aggregate with an afterok dependency on the array.
   V1 and v2 completion schemas are scientifically ineligible.
