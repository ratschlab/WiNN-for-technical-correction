# Frozen WiNN build

The benchmark is frozen to WiNN 0.1.6 at Git commit
`50b72a4ae7fc4093c4c624ebd0fb512bff58fa7a`.

The exact source archive is `winn_0.1.6.tar.gz`, built with
`git archive 50b72a4 | R CMD build --no-build-vignettes --no-manual winn`, with SHA-256:

```text
e40a2c3fd709f11a456ac3721379da3549e89c619e5a164743e203fe1e9be7ab
```

Version 0.1.6 fixes the drift-smoothing formulas. Versions 0.1.4 and 0.1.5 wrote the
smooth term as `mgcv::s()`, which mgcv's formula parser does not recognise; every spline
fit errored inside `tryCatch()`, so drift correction fell back to LOESS in the
conservative mode and to a linear fit in the standard mode. Version 0.1.6 uses `s()`.
No exported function or argument changed.

Validation of the 0.1.6 release:

- All 105 package test expectations passed, with no failures, warnings or skips,
  including new regression tests that require the documented smoother
  (`conservative_ps_reml`, `standard_gam_cr`) for corrected features; these tests fail
  on 0.1.5.
- `R CMD check --as-cran` on macOS (R 4.6.1) completed with 0 errors, 0 warnings and
  0 notes; examples (including `--run-donttest`), vignettes and the PDF and HTML
  manuals were checked. The GitHub Actions matrix is recorded in the package's
  `cran-comments.md`.
- CRAN accepted the 0.1.6 update on 11 November 2026.

Every production job must load the isolated release installation, assert version
0.1.6, verify that the installed package path is inside the frozen release library,
and verify the source-archive hash before analysis begins (`assert_frozen_winn.R`).

## Previous build

The WiNN 0.1.4 release (commit `a6130a123a5a56a089a0e3a0465e55bbf31e0cc6`, archive
`winn_0.1.4.tar.gz`, SHA-256
`71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1`) produced the
comparator outputs that are reused in the current results. Its archive is kept for
provenance. Before that freeze, 90 package tests and `R CMD check --as-cran` passed,
and GitHub Actions run 30178559865 passed on Ubuntu, Windows and macOS.
