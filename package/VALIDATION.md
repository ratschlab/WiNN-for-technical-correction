# Frozen permutation-calibration build

The revised WiNN results use development version 0.1.7.9999 at package commit
`ab4969792c98bd86208a28a4e83b7a69e204b9fc`. The source archive
`winn_0.1.7.9999.tar.gz` has SHA-256:

```text
72c90896ac28f3ee4000a5c6d5287e7d4031e307c40b787d001755a6982c4c29
```

The package retains the Ljung–Box statistic and calibrates its default drift gate
with 999 joint-column permutations, maxT adjustment and a plus-one correction.
Its default permutation seed is 104729; adjusted probabilities are compared directly
with 0.05 without another BH adjustment. The legacy asymptotic/BH rule remains an
explicit option. The batch ANOVA gate and other correction stages are unchanged.

The full package test suite passed 135 expectations. `R CMD check --no-manual` of
the source archive with built vignettes completed with 0 errors, 0 warnings and 0
notes. This development build is not the CRAN 0.1.6 release.

The same package commit passed [GitHub Actions](https://github.com/ratschlab/winn/actions/runs/43591609421); the run head was verified against the frozen commit.

Production uses the isolated Euler R 4.5.1 installation. Before analysis, each job
checks the version, installed-library location and tarball SHA-256 against
`frozen_package.json`. Independent gate validation also hashes the formals and
bodies of both calibration helpers to verify that its source matches production.
Bootstrap/PCA derivatives additionally pin `irlba` 2.3.7; its source URL and hash
are recorded in `analysis/config/calibration_revision/derivative_dependencies.csv`.

# Historical WiNN 0.1.6 build

The previous benchmark was frozen to WiNN 0.1.6 at Git commit
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
- CRAN published the 0.1.6 update on 10 November 2026.

Every historical 0.1.6 production job had to load the isolated release installation, assert version
0.1.6, verify that the installed package path is inside the frozen release library,
and verify the source-archive hash before analysis begins (`assert_frozen_winn.R`).

## Previous build

The WiNN 0.1.4 release (commit `a6130a123a5a56a089a0e3a0465e55bbf31e0cc6`, archive
`winn_0.1.4.tar.gz`, SHA-256
`71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1`) produced the
comparator outputs that are reused in the current results. Its archive is kept for
provenance. Before that freeze, 90 package tests and `R CMD check --as-cran` passed,
and GitHub Actions run 30178559865 passed on Ubuntu, Windows and macOS.
