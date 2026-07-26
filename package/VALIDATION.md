# Frozen WiNN build

The benchmark is frozen to WiNN 0.1.4 at Git commit
`a6130a123a5a56a089a0e3a0465e55bbf31e0cc6`.

The exact source archive is `winn_0.1.4.tar.gz`, with SHA-256:

```text
71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1
```

Validation completed before the benchmark freeze:

- 90 package tests passed.
- The run-order permutation-invariance regression test passed.
- The independent two-tail MAD regression tests passed.
- The `ljung_box_fitdf = 1` compatibility fixture was unchanged.
- `R CMD check --as-cran` on the exact source archive completed with 0 errors,
  0 warnings, and the expected single new-submission note.
- HTML and PDF manuals built successfully.
- GitHub Actions run 30178559865 passed on Ubuntu R-release, R-devel, and
  R-oldrel, Windows R-release, and macOS R-release.

Every production job must load the isolated Euler installation, assert version
0.1.4, verify that the installed package path is inside the frozen release
library, and verify the source-archive hash before analysis begins.
