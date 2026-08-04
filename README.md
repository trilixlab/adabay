# adabay <img src="man/figures/logo.png" align="right" height="139" alt="adabay logo" />

`adabay` is an R package for the rapid evaluation and calibration of Bayesian group sequential designs (GSDs) across continuous, binary, count and time-to-event endpoints. It implements the semi-simulation framework described in the companion manuscript: virtual trials are simulated by Monte Carlo, but per-look posteriors are computed in closed form by approximating any user-specified prior with a finite mixture of conjugate components.

## Installation

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("trilixlab/adabay")
```

or from the source tarball distributed with the manuscript's supplementary
materials:

```r
install.packages("adabay_0.1.0.tar.gz", repos = NULL, type = "source")
```

The optional `RBesT` package is recommended for fitting non-conjugate priors via `fit_mixture()`. The package depends on `Rcpp`, `parallel`, `graphics` and `stats`, and uses `ggplot2` and `ggsci` for plot methods if available.

## Workflow

A typical workflow has four steps: build a design, prior, and decision rule, then call `evaluate_design()`. The same function is the entry point for both the single-design fused path (pass an `adabay_design`) and the cached grid-evaluation path (pass an `adabay_cache` produced by `build_cache()`).

```r
library(adabay)

des <- set_design(
  endpoint     = "binary",
  n_per_look   = c(760, 1520, 2280, 3040, 3800),
  effect_scale = "risk_difference",
  alternative  = "less"
)

pri <- set_prior(
  endpoint = "binary",
  arms = list(c = list(family = "beta", a = 1, b = 1),
              t = list(family = "beta", a = 1, b = 1))
)

dec <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE
)

oc <- evaluate_design(des,
  prior    = pri,
  decision = dec,
  effect   = list(theta_c = 0.33, theta_t = 0.28),
  n_trials = 1e6,
  cores    = 8,
  seed     = 1L
)

print(oc)
plot(oc)
plot(oc, file = "adrenal.jpeg", width = 7, height = 4, dpi = 300)
```

## Supported endpoints

`endpoint` takes one of four values, given in the first column below.

| `endpoint`     | Outcome type  | Conjugate kernel              | Effect scales |
|----------------|---------------|-------------------------------|---------------|
| `"continuous"` | Continuous    | Normal, Normal--inverse-gamma | mean_difference, standardised_mean_difference |
| `"binary"`     | Binary        | Beta                          | risk_difference, risk_ratio, odds_ratio       |
| `"count"`      | Count         | Gamma                         | rate_difference, rate_ratio, log_rate_ratio   |
| `"tte"`        | Time-to-event | Gamma                         | hazard_ratio, log_hazard_ratio                |

Time-to-event designs are specified in pooled event counts (`d_total`, `d_per_look`) and require an accrual specification via `set_accrual()`, because the calendar-time cutoff at each event-driven look is driven by the simulated arrival times.

## Top-level functions (verb_object snake_case)

| Function            | Purpose |
|---------------------|---------|
| `set_design()`      | Specify the design skeleton (endpoint, look schedule, scale) |
| `set_prior()`       | Specify per-arm conjugate or mixture-of-conjugate priors |
| `set_decision()`    | Specify per-look efficacy and futility decision rules |
| `set_accrual()`     | Specify recruitment, follow-up and dropout |
| `fit_mixture()`     | Approximate a non-conjugate prior by a conjugate mixture |
| `build_cache()`     | Cache per-look tail probabilities for grid evaluation |
| `evaluate_design()` | Evaluate operating characteristics; polymorphic on an `adabay_design` (fused simulate-and-aggregate path) or an `adabay_cache` (cheap aggregation from a pre-built cache) |
| `calibrate_design()`| Search a threshold grid for target operating characteristics |
| `summarise_oc()`    | Tidy summary across designs in a grid |

## Worked examples

Reproducible scripts for the four case studies in the manuscript are provided under `inst/examples/`:

* `transform2.R` -- continuous endpoint (TRANSFORM-2 esketamine for treatment-resistant depression)
* `adrenal.R` -- binary endpoint (ADRENAL re-design)
* `bnt162b2.R` -- count endpoint (Pfizer--BioNTech BNT162b2 COVID-19 vaccine trial)
* `checkmate141.R` -- time-to-event endpoint (CheckMate-141 head and neck cancer trial)

## Performance

The slow per-trial posterior integrals on the fixed-node quadrature scales (`risk_difference`, `risk_ratio` and `odds_ratio` for binary, `rate_difference` for count) are implemented in C++ via Rcpp using a 64-node Gauss--Legendre quadrature; the continuous normal--inverse-gamma path evaluates its inner integral with `stats::integrate()` at the R level, and the remaining scales are closed-form. The trial-level loop is parallel-friendly through `parallel::mclapply` (POSIX) or `parallel::parLapply` (Windows), with per-worker streams from the L'Ecuyer combined multiple-recursive generator.

Results are bit-identical for any fixed `(seed, cores)` pair, so reproducing a published number requires matching both. The worked examples and the manuscript use `seed = 1L` and `cores = 8`. Bit-identity *across* worker counts is available by setting `cross_core_reproducible = TRUE` on `evaluate_design()` or `build_cache()`, which forces the simulator pass into a single driver process; the default is `FALSE`.

Cross-package wall-clock and per-virtual-trial timing comparisons are reported in Section 7 of the companion software paper (arXiv:2608.01068), across all four case studies. Comparator coverage differs by endpoint: gsbDesign provides the analytic reference for the continuous case only; adaptr covers the continuous and binary cases; BATSS covers the continuous, binary and count cases; and no faithful comparator is constructible for the time-to-event case, because no existing tool provides an equivalent survival-likelihood implementation under the conventions used here.

## Plot styling

`plot.adabay_oc()` uses `ggplot2::theme_bw()` and the BMJ palette via `ggsci::scale_fill_bmj()` (with a sensible fallback when `ggsci` is not installed). When the `file` argument is supplied, the plot is written via `ggplot2::ggsave()`; the file extension determines the device (`.jpeg`, `.jpg`, `.png`, `.pdf`).

## Citation

If you use `adabay` in published work, please cite both the software paper and
the companion methodology paper:

* He Z, Yu F (2026). *adabay: an R package for rapid evaluation and calibration
  of Bayesian group sequential designs across common endpoint types.*
  R package version 0.1.0. arXiv preprint arXiv:2608.01068.
  <https://arxiv.org/abs/2608.01068>
* He Z, Yu F, Cro S, Billot L (2026). *Rapid evaluation and calibration of
  Bayesian group sequential designs via conjugate-mixture semi-simulation.*
  arXiv preprint arXiv:2607.22900. <https://arxiv.org/abs/2607.22900>

`citation("adabay")` returns both entries; `toBibtex(citation("adabay"))` gives them in BibTeX form.

## Licence

MIT (see [LICENSE](LICENSE)).
