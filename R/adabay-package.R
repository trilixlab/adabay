#' adabay: rapid evaluation and calibration of Bayesian group sequential designs across common endpoint types
#'
#' \pkg{adabay} implements a semi-simulation framework for the rapid
#' evaluation and calibration of Bayesian group sequential designs across
#' continuous, binary, count and time-to-event endpoints. Virtual trials are
#' simulated by Monte Carlo, and per-look posteriors are computed in closed
#' form by approximating any user-specified prior with a finite mixture of
#' conjugate components.
#'
#' A typical workflow consists of building a design, prior and decision
#' specification with [set_design()], [set_prior()] and [set_decision()],
#' then passing them directly to [evaluate_design()] for a single-design
#' run, or to [build_cache()] followed by [evaluate_design()] (or
#' [calibrate_design()]) for grid-based threshold and look-time calibration.
#'
#' @section Top-level functions:
#' * [set_design()] -- specify the design skeleton
#' * [set_prior()] -- specify per-arm priors
#' * [set_decision()] -- specify per-look decision rules
#' * [set_accrual()] -- specify recruitment, follow-up and dropout
#' * [fit_mixture()] -- approximate a non-conjugate prior by a conjugate mixture
#' * [build_cache()] -- cache per-look tail probabilities over a grid
#' * [evaluate_design()] -- evaluate operating characteristics; polymorphic on
#'   a `adabay_design` (fused simulate-and-aggregate path) or a `adabay_cache`
#'   (cheap aggregation from a pre-built cache)
#' * [calibrate_design()] -- search a threshold grid for target operating characteristics
#' * [summarise_oc()] -- tidy summary across designs in a grid
#'
#' @author
#' Zhangyi He \email{zhe1@georgeinstitute.org.uk}
#'
#' Feng Yu \email{feng.yu@bristol.ac.uk} (maintainer)
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline barplot hist legend lines par
#' @importFrom parallel clusterSetRNGStream makeCluster mclapply parLapply stopCluster
#' @importFrom stats approx dbeta density dgamma dnorm dt integrate pbeta pgamma pnorm pt quantile rbinom rnorm rpois runif uniroot var
#' @useDynLib adabay, .registration = TRUE
"_PACKAGE"

## Quiet R CMD check NSE notes from ggplot2 aesthetics in plot.adabay_oc().
utils::globalVariables(c("look", "prob", "type", ".data"))
