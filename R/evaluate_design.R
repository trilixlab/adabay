#' Evaluate operating characteristics for a Bayesian group sequential design
#'
#' `evaluate_design()` is the single entry point for obtaining the
#' operating characteristics of a Bayesian group sequential design. It is
#' polymorphic in its first argument:
#'
#' \describe{
#'   \item{`adabay_design`}{Fused path: simulates `n_trials` virtual trials at
#'     the supplied data-generating values, computes per-look posterior tail
#'     probabilities of \eqn{\Delta} in closed form, and applies the supplied
#'     decision rule to obtain the per-look stopping probabilities, the
#'     overall type I error rate (or power under the alternative), and the
#'     expected sample size, all in one call. Use this for one-off evaluation
#'     of a single design.}
#'   \item{`adabay_cache`}{Cached path: looks up pre-computed per-look tail
#'     probabilities in a `adabay_cache` produced once by [build_cache()] and
#'     aggregates them under the supplied decision rule, without re-simulating
#'     trials. Use this for grid-style calibration when many decision rules
#'     share the same simulated trials.}
#' }
#'
#' Slow inner integrals are dispatched to C++ implementations via
#' \pkg{Rcpp}; the trial-level loop is embarrassingly parallel and is
#' parallelised across `cores` workers via [parallel::mclapply()] (POSIX) or
#' [parallel::parLapply()] (Windows).
#'
#' @param x A `adabay_design` object (fused path) or a `adabay_cache` object
#'   (cached path).
#' @param prior For the fused path only: a `adabay_prior` from [set_prior()].
#' @param decision A `adabay_decision` from [set_decision()].
#' @param effect For the fused path only: named list of data-generating
#'   values; see the worked examples in the package vignettes for endpoint-
#'   specific names.
#' @param accrual For the fused path only: a [set_accrual()] specification.
#'   **Required** for tte designs, which are always staggered: the
#'   recruitment rate fixes the Poisson enrolment process whose arrival times
#'   drive the calendar-time exposure cutoff. The tte endpoint currently
#'   supports only `model = "poisson"`. Optional for the other endpoints (continuous,
#'   binary, count), which have no simulated recruitment timeline, where it
#'   is used only to report an analytic expected duration.
#' @param n_trials For the fused path only: number of virtual trials.
#' @param cores For the fused path only: number of CPU cores. Defaults to 1.
#' @param seed For the fused path only: optional integer seed.
#' @param cross_core_reproducible For the fused path only: logical, default
#'   `FALSE`. When `TRUE`, the simulator is forced to run in a single driver
#'   process (the user-supplied `cores` is ignored for the simulation pass)
#'   so that the trial-level RNG draws are bit-identical across worker
#'   counts. The posterior and aggregation passes are core-independent in
#'   both modes. Set to `TRUE` for regulatory submissions or other settings
#'   where reproducibility across hardware configurations is required; leave
#'   `FALSE` for design exploration where wall-clock matters more than
#'   cross-`cores` bit-identity. The cost is the loss of simulator
#'   parallelism: small for binary, count and tte (where the simulator
#'   is a small fraction of total wall-clock); moderate for continuous.
#' @param look_subset For the cached path only: optional integer vector. If
#'   supplied, the design is evaluated at the subset of look indices given
#'   (in order). Defaults to all looks in the cached design.
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return An object of class `"adabay_oc"`, a list with elements:
#'   \describe{
#'     \item{`alpha`}{Overall type I error rate scalar (regulatory
#'       convention: under non-binding futility, computed as if futility
#'       were absent).}
#'     \item{`eff_prob`}{Overall efficacy probability scalar under binding
#'       trial conduct.}
#'     \item{`fut_prob`}{Overall futility probability scalar under binding
#'       trial conduct.}
#'     \item{`eff_prob_at_look`, `fut_prob_at_look`}{Per-look stopping
#'       probabilities under binding conduct; they partition the sample
#'       space.}
#'     \item{`alpha_at_look_binding`, `alpha_at_look_nonbinding`}{Alpha
#'       spent at each look under binding and non-binding semantics
#'       respectively.}
#'     \item{`expected_sample_size`}{Expected total subjects (continuous,
#'       binary), exposure (count) or events (tte) under binding
#'       conduct.}
#'     \item{`expected_duration`}{Expected calendar duration; present only
#'       when an [set_accrual()] specification is supplied (fused path) or
#'       was supplied to [build_cache()] (cached path), and absent otherwise.
#'       For tte this is the exact Monte Carlo mean of the calendar
#'       time realised (during simulation) at each trial's own stopping
#'       look; for continuous, binary and count -- which have no simulated
#'       recruitment timeline -- it is an analytic approximation.}
#'     \item{`tau`, `C`}{Per-trial stopping look and trial outcome code
#'       (`+1`/`-1`/`0` for efficacy/futility/no-cross).}
#'     \item{`n_trials`, `effect`, `design`, `decision`}{Echoes of the inputs
#'       for downstream reproducibility. The fused path additionally records
#'       `seed`, `call` and (when supplied) `accrual`; these are omitted on
#'       the cached path.}
#'   }
#'
#' @section Binding vs non-binding futility:
#'   The simulator always conducts virtual trials under binding semantics
#'   (a trial physically stops at the first crossing of either the efficacy
#'   or futility threshold). The per-look stopping probabilities
#'   `eff_prob_at_look[k]` and `fut_prob_at_look[k]` therefore always reflect
#'   actual binding conduct and partition the sample space:
#'   `sum(eff_prob_at_look) + sum(fut_prob_at_look) = 1` for both binding
#'   and non-binding settings of `futility_binding`.
#'
#'   The flag changes only the overall scalar type I error rate `alpha`:
#'   \describe{
#'     \item{Binding}{`alpha = sum(eff_prob_at_look)`: probability that the
#'       trial crosses efficacy before any futility crossing under the
#'       data-generating distribution.}
#'     \item{Non-binding}{`alpha = sum(alpha_at_look_nonbinding)`: probability
#'       that the trial would have crossed efficacy at some look had the
#'       futility rule been ignored (the standard regulatory convention).}
#'   }
#'
#' @examples
#' des <- set_design(endpoint = "binary",
#'                   n_per_look = c(760, 1520, 2280, 3040, 3800),
#'                   effect_scale = "risk_difference",
#'                   alternative = "less")
#' pri <- set_prior(endpoint = "binary",
#'                  arms = list(c = list(family = "beta", a = 1, b = 1),
#'                              t = list(family = "beta", a = 1, b = 1)))
#' dec <- set_decision(
#'   efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
#'   futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
#'   futility_binding = FALSE)
#' \donttest{
#' ## Fused path: evaluate a single design directly
#' oc <- evaluate_design(des, prior = pri, decision = dec,
#'                       effect = list(theta_c = 0.33, theta_t = 0.28),
#'                       n_trials = 200L, cores = 1L, seed = 1L)
#' print(oc)
#'
#' ## Cached path: build once, evaluate many decision rules cheaply
#' cache <- build_cache(des, pri,
#'                      effect         = list(theta_c = 0.33, theta_t = 0.28),
#'                      n_trials       = 200L,
#'                      seed           = 1L,
#'                      threshold_grid = list(efficacy = 0, futility = 0))
#' oc2 <- evaluate_design(cache, dec)
#' print(oc2)
#' }
#' @export
evaluate_design <- function(x, ...) UseMethod("evaluate_design")

#' @rdname evaluate_design
#' @export
evaluate_design.default <- function(x, ...) {
  stop("'x' must be a adabay_design object (fused path) or a adabay_cache ",
       "object (cached path), not an object of class ",
       sprintf("'%s'.", paste(class(x), collapse = "/")), call. = FALSE)
}

#' @rdname evaluate_design
#' @export
evaluate_design.adabay_design <- function(x, prior, decision, effect,
                                          accrual = NULL,
                                          n_trials = 10000L,
                                          cores = 1L, seed = NULL,
                                          cross_core_reproducible = FALSE,
                                          ...) {
  design <- x
  .check_inputs(design, prior, decision)
  .warn_quadrature_shape(design, prior)
  if (!is.null(accrual) && !inherits(accrual, "adabay_accrual"))
    stop("'accrual' must be a adabay_accrual object.", call. = FALSE)
  .check_tte_accrual(design, accrual)
  n_trials <- as.integer(n_trials)
  if (n_trials < 1L)
    stop("'n_trials' must be a positive integer.", call. = FALSE)

  sim   <- .simulate_data(design, effect, n_trials = n_trials,
                          seed = seed, cores = cores,
                          cross_core_reproducible = cross_core_reproducible,
                          accrual = accrual)
  tails <- .compute_tail_probs(design, prior, sim, decision, cores = cores)
  oc    <- .aggregate_oc(design, decision, tails, sim)

  if (!is.null(accrual))
    oc$expected_duration <- .expected_duration(design, oc, accrual, sim = sim)

  oc$call     <- match.call()
  oc$n_trials <- n_trials
  oc$effect   <- effect
  oc$accrual  <- accrual
  oc$seed     <- seed
  oc
}

#' @rdname evaluate_design
#' @export
evaluate_design.adabay_cache <- function(x, decision,
                                         look_subset = NULL, ...) {
  cache <- x
  if (!inherits(decision, "adabay_decision"))
    stop("'decision' must be a adabay_decision object.", call. = FALSE)
  n_looks_full <- cache$design$n_looks
  if (is.null(look_subset)) look_subset <- seq_len(n_looks_full)
  look_subset <- as.integer(look_subset)
  if (any(look_subset < 1L | look_subset > n_looks_full))
    stop("'look_subset' must lie within 1..n_looks of the cached design.",
         call. = FALSE)

  design_eff <- cache$design
  design_eff$n_looks <- length(look_subset)
  if (cache$design$endpoint %in% c("continuous", "binary")) {
    design_eff$schedule$n_c <- cache$design$schedule$n_c[look_subset]
    design_eff$schedule$n_t <- cache$design$schedule$n_t[look_subset]
  } else if (cache$design$endpoint == "count") {
    design_eff$schedule$E_c <- cache$design$schedule$E_c[look_subset]
    design_eff$schedule$E_t <- cache$design$schedule$E_t[look_subset]
  } else if (cache$design$endpoint == "tte") {
    design_eff$schedule$D <- cache$design$schedule$D[look_subset]
  }

  thresholds_full <- cache$thresholds
  thresholds_grid <- .threshold_grid(decision, n_looks = design_eff$n_looks)
  tails <- vector("list", design_eff$n_looks)
  for (kk in seq_along(look_subset)) {
    look <- look_subset[kk]
    full <- cache$tails[[look]]
    cols <- match(thresholds_grid[[kk]], thresholds_full)
    if (any(is.na(cols)))
      stop("Decision rule references thresholds not in the cache. Re-run build_cache() with the new thresholds.",
           call. = FALSE)
    tails[[kk]] <- full[, cols, drop = FALSE]
  }
  attr(tails, "thresholds") <- thresholds_grid

  sim <- .slice_sim(cache$sim, look_subset)
  oc  <- .aggregate_oc(design_eff, decision, tails, sim)
  oc$n_trials <- cache$n_trials
  oc$effect   <- cache$effect
  if (!is.null(cache$accrual)) {
    oc$accrual           <- cache$accrual
    oc$expected_duration <- .expected_duration(design_eff, oc, cache$accrual,
                                               sim = sim)
  }
  oc
}

## Slice the cached per-trial sim data down to a look_subset.
.slice_sim <- function(sim, idx) {
  e <- sim$endpoint
  if (e == "continuous") {
    sim$ybar_c <- sim$ybar_c[, idx, drop = FALSE]
    sim$ybar_t <- sim$ybar_t[, idx, drop = FALSE]
    sim$S2_c   <- sim$S2_c[, idx, drop = FALSE]
    sim$S2_t   <- sim$S2_t[, idx, drop = FALSE]
    sim$n_c    <- sim$n_c[idx]
    sim$n_t    <- sim$n_t[idx]
  } else if (e == "binary") {
    sim$S_c <- sim$S_c[, idx, drop = FALSE]
    sim$S_t <- sim$S_t[, idx, drop = FALSE]
    sim$n_c <- sim$n_c[idx]
    sim$n_t <- sim$n_t[idx]
  } else if (e == "count") {
    sim$S_c <- sim$S_c[, idx, drop = FALSE]
    sim$S_t <- sim$S_t[, idx, drop = FALSE]
    sim$E_c <- sim$E_c[idx]
    sim$E_t <- sim$E_t[idx]
  } else if (e == "tte") {
    sim$D_c <- sim$D_c[, idx, drop = FALSE]
    sim$D_t <- sim$D_t[, idx, drop = FALSE]
    sim$E_c <- sim$E_c[, idx, drop = FALSE]
    sim$E_t <- sim$E_t[, idx, drop = FALSE]
    sim$C_look <- sim$C_look[, idx, drop = FALSE]
    sim$d_per_look <- sim$d_per_look[idx]
  }
  sim
}
