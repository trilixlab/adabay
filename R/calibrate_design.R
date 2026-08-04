#' Calibrate design thresholds against operating-characteristic targets
#'
#' Sweeps a grid of efficacy and futility posterior-probability thresholds
#' against a precomputed cache, evaluating each candidate design and
#' returning the design that minimises the expected sample size under the
#' alternative subject to the type I error and power constraints.
#'
#' This is a convenience wrapper around [evaluate_design()] that constructs
#' the decision rule for each combination of efficacy threshold (\code{p})
#' and futility threshold (\code{q}) over the supplied grids, both applied
#' at every look at the null effect-size value.
#'
#' For more flexible calibration (per-look thresholds, dual-criterion rules,
#' alternative effect-size cuts), call [evaluate_design()] directly within a
#' user-defined search loop.
#'
#' @param cache An object of class \code{"adabay_cache"} cached under \eqn{H_0}
#'   (drives the type I error rate of every grid point).
#' @param cache_alt A second \code{"adabay_cache"} cached under \eqn{H_1}
#'   (drives the power and expected sample size of every grid point).
#'   Required: a calibration is a constrained search over both targets and
#'   cannot be carried out from a single cache.
#' @param alpha_target Required upper bound on the type I error rate.
#' @param power_target Required lower bound on the power (probability of crossing efficacy under the alternative effect).
#' @param efficacy_grid Numeric vector of efficacy posterior-probability
#'   thresholds to search.
#' @param futility_grid Numeric vector of futility posterior-probability
#'   thresholds to search. Use \code{c(1)} to disable futility.
#' @param futility_binding Logical. Defaults to \code{FALSE}.
#' @return An object of class \code{"adabay_calibration"}: a list with
#'   elements
#'   \describe{
#'     \item{\code{best}}{An object of class \code{"adabay_calibration_best"}
#'       (or \code{NULL}, with a warning, if no grid point meets both
#'       targets) for the \code{(p, q)} combination that meets the targets
#'       at minimum expected sample size under the alternative. Carries
#'       \code{p}, \code{q}, \code{type_I} (from the \eqn{H_0} cache),
#'       \code{power} and \code{expected_sample_size} (both from the
#'       \eqn{H_1} cache) as plain numeric fields, plus the full
#'       \code{adabay_oc} objects for each hypothesis as \code{oc_h0} /
#'       \code{oc_h1} for deeper inspection (e.g. \code{plot(cal$best$oc_h1)}
#'       for per-look stopping probabilities). Deliberately not of class
#'       \code{"adabay_oc"} itself: \code{oc_h1$alpha} is the *power* (it
#'       was evaluated under the alternative), and printing it through
#'       \code{print.adabay_oc()} would mislabel that as "Type I error".}
#'     \item{\code{grid}}{A data frame of every grid point evaluated, with
#'       columns \code{p}, \code{q}, \code{type_I}, \code{power}, \code{E_N}.}
#'   }
#' @examples
#' \donttest{
#'   des <- set_design(endpoint     = "binary",
#'                     n_per_look   = c(760, 1520, 2280),
#'                     effect_scale = "risk_difference",
#'                     alternative  = "less")
#'   pri <- set_prior(endpoint = "binary",
#'                    arms = list(c = list(family = "beta", a = 1, b = 1),
#'                                t = list(family = "beta", a = 1, b = 1)))
#'   cache_h0 <- build_cache(des, pri,
#'                           effect   = list(theta_c = 0.33, theta_t = 0.33),
#'                           n_trials = 200L, seed = 1L,
#'                           threshold_grid = list(efficacy = 0, futility = 0))
#'   cache_h1 <- build_cache(des, pri,
#'                           effect   = list(theta_c = 0.33, theta_t = 0.28),
#'                           n_trials = 200L, seed = 1L,
#'                           threshold_grid = list(efficacy = 0, futility = 0))
#'   cal <- calibrate_design(cache          = cache_h0,
#'                           cache_alt      = cache_h1,
#'                           alpha_target   = 0.05,
#'                           power_target   = 0.50,
#'                           efficacy_grid  = c(0.95, 0.99),
#'                           futility_grid  = c(0.80, 0.90))
#'   head(cal$grid)
#' }
#' @export
calibrate_design <- function(cache,
                             cache_alt,
                             alpha_target = 0.025,
                             power_target = 0.90,
                             efficacy_grid,
                             futility_grid = c(1),
                             futility_binding = FALSE) {
  if (!inherits(cache, "adabay_cache"))
    stop("'cache' must be a adabay_cache object.", call. = FALSE)
  if (missing(cache_alt) || is.null(cache_alt))
    stop("'cache_alt' is required: pass a adabay_cache built under the alternative effect.",
         call. = FALSE)
  if (!inherits(cache_alt, "adabay_cache"))
    stop("'cache_alt' must be a adabay_cache object.", call. = FALSE)
  d0 <- cache$design
  d1 <- cache_alt$design
  same_design <- identical(d0$endpoint, d1$endpoint) &&
    identical(d0$n_looks, d1$n_looks) &&
    identical(d0$schedule, d1$schedule) &&
    identical(d0$effect_scale, d1$effect_scale) &&
    identical(d0$delta_null, d1$delta_null)
  if (!same_design)
    stop("'cache' and 'cache_alt' must be built from the same design ",
         "(endpoint, look schedule, n_looks, effect scale and delta_null); ",
         "they currently describe different designs, so the type I error and ",
         "the power/expected sample size would refer to different trials.",
         call. = FALSE)

  delta_null <- cache$design$delta_null
  records <- list()
  best_en <- Inf
  best <- NULL

  for (p in efficacy_grid) {
    for (q in futility_grid) {
      dec <- if (q >= 1) {
        set_decision(
          efficacy = list(list(threshold_effect = delta_null,
                               threshold_prob = p)),
          futility = NULL,
          futility_binding = futility_binding)
      } else {
        set_decision(
          efficacy = list(list(threshold_effect = delta_null,
                               threshold_prob = p)),
          futility = list(list(threshold_effect = delta_null,
                               threshold_prob = q)),
          futility_binding = futility_binding)
      }
      oc_h0 <- evaluate_design(cache, dec)
      oc_h1 <- evaluate_design(cache_alt, dec)
      records[[length(records) + 1L]] <- data.frame(
        p = p, q = q,
        type_I = oc_h0$alpha,
        power  = oc_h1$alpha,
        E_N    = oc_h1$expected_sample_size)
      if (isTRUE(oc_h0$alpha <= alpha_target) &&
          isTRUE(oc_h1$alpha >= power_target) &&
          oc_h1$expected_sample_size < best_en) {
        best_en <- oc_h1$expected_sample_size
        best <- structure(
          list(p = p, q = q,
               type_I = oc_h0$alpha,
               power  = oc_h1$alpha,
               expected_sample_size = oc_h1$expected_sample_size,
               oc_h0 = oc_h0, oc_h1 = oc_h1),
          class = "adabay_calibration_best")
      }
    }
  }
  if (is.null(best))
    warning("calibrate_design(): no threshold combination met the alpha and ",
            "power targets; returning best = NULL. Inspect 'grid' and relax ",
            "the targets or widen the efficacy/futility grids.", call. = FALSE)
  structure(list(best = best, grid = do.call(rbind, records)),
            class = "adabay_calibration")
}

#' @noRd
#' @export
print.adabay_calibration <- function(x, ...) {
  cat("Threshold calibration\n")
  cat(sprintf("  Grid evaluated: %d threshold combination%s (see $grid)\n",
              nrow(x$grid), if (nrow(x$grid) == 1L) "" else "s"))
  if (is.null(x$best)) {
    cat("  No combination met the targets; $best is NULL.\n")
  } else {
    cat("\n")
    print(x$best)
  }
  invisible(x)
}

#' @noRd
#' @export
print.adabay_calibration_best <- function(x, ...) {
  cat("Calibrated design\n")
  cat(sprintf("  Efficacy threshold (p):  %.4f\n", x$p))
  cat(sprintf("  Futility threshold (q):  %s\n",
              if (x$q >= 1) "disabled" else sprintf("%.4f", x$q)))
  cat(sprintf("  Type I error (H0):       %.4f\n", x$type_I))
  cat(sprintf("  Power (H1):              %.4f\n", x$power))
  cat(sprintf("  Expected sample size (H1): %.1f\n", x$expected_sample_size))
  cat("  Note: Type I error and Power follow the calibration's futility\n")
  cat("        convention; under the default non-binding rule the efficacy\n")
  cat("        probability ignores futility crossings. See $oc_h0 / $oc_h1.\n")
  cat("  Full per-hypothesis operating characteristics: $oc_h0, $oc_h1\n")
  invisible(x)
}
