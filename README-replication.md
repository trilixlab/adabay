# Code/ — adabay package and replication script

This folder holds the package source, the built tarball, the reference manual
PDF and the standalone replication script for the manuscript

> He Z, Yu F (2026). *adabay: An R package for rapid evaluation and
> calibration of Bayesian group sequential designs across common endpoint
> types.*

## What to install

The submission version is **adabay 0.1.0**:

* `adabay_0.1.0.tar.gz` — source tarball.
* `adabay/` — source tree (the contents of the tarball, unpacked).
* `adabay-manual.pdf` — PDF reference manual generated from the v0.1.0
  source.

To install from the bundled tarball:

```r
install.packages("adabay_0.1.0.tar.gz", repos = NULL, type = "source")
```

The manuscript, the supplemental material, the replication script, the
reference manual, and the rendered website (top-level `Website/` directory)
all correspond to **v0.1.0**.

## Replication

`replication.R` reproduces every numerical result reported in the manuscript.
It is controlled by a single `RUN_FULL_PRECISION` flag (smoke test at
`R = 2,000` versus the manuscript budgets at `R = 10,000` matched and
`R = 1,000,000` high precision) plus three environment-variable gates
(`RUN_ADAPTR`, `RUN_GSBDESIGN`, `RUN_BATSS`) for the comparator runs of
Section 7. See the script header for details.

## A note on the rendered website (`Website/`)

The top-level `Website/` directory of the supplementary materials is a
pkgdown-rendered HTML snapshot of v0.1.0 frozen at the build date listed
in `Website/pkgdown.yml`. It mirrors the source contents of
`adabay/` to the resolution of that build. Function-reference
pages (`Website/reference/evaluate_design.html`,
`Website/reference/build_cache.html`, etc.) show the v0.1.0 signatures as
they appear in `man/*.Rd`. The canonical source for the latest API
surface — including any post-build refinements such as the
`cross_core_reproducible` argument — is the package source (`R/` and
`man/`), the manual PDF (`adabay-manual.pdf`), the package vignettes,
`NEWS.md`, and the manuscript Sections 4.6 and 7.5. Reviewers comparing
the source against the rendered site should treat the source as
authoritative.
