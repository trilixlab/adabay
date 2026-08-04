#' Set the Bayesian group sequential design skeleton
#'
#' Declares the endpoint type, the look schedule, and the pooled (both-arms)
#' sample size, exposure or event-count target at each look, plus the effect
#' scale and the direction of the alternative hypothesis.
#'
#' @param endpoint One of \code{"continuous"}, \code{"binary"}, \code{"count"}
#'   or \code{"tte"} (time-to-event).
#' @param n_per_look Cumulative sample size at each look, pooled across both
#'   arms (continuous and binary endpoints). The number of looks \eqn{K} is
#'   taken from \code{length(n_per_look)} -- there is no separate
#'   \code{n_looks} argument. Split into per-arm targets by
#'   \code{allocation_ratio}:
#'   \code{n_c[k] = round(n_per_look[k] / (1 + allocation_ratio))},
#'   \code{n_t[k] = n_per_look[k] - n_c[k]}, so \code{n_c[k] + n_t[k] ==
#'   n_per_look[k]} exactly at every look. With the default
#'   \code{allocation_ratio = 1} this is a 1:1 split (up to rounding).
#'   Errors if consecutive \code{n_per_look} values are too close together
#'   for the split to be strictly increasing in both arms at every look.
#' @param exposure_per_look Cumulative exposure at each look, pooled across
#'   both arms (count endpoint); \eqn{K} is taken from
#'   \code{length(exposure_per_look)}. Split into per-arm targets by
#'   \code{allocation_ratio} the same way as \code{n_per_look} (exactly,
#'   with no rounding, since exposure is not required to be an integer).
#' @param d_total Total target number of events for the tte endpoint,
#'   pooled across both arms; unaffected by \code{allocation_ratio}.
#'   Defaults to \code{NULL}, in which case the final entry of
#'   \code{d_per_look} is used.
#' @param d_per_look Cumulative number of events at each look (tte
#'   endpoint), pooled across both arms; \eqn{K} is taken from
#'   \code{length(d_per_look)}. Unaffected by \code{allocation_ratio}.
#' @param effect_scale Effect-size parameterisation. One of
#'   \code{"mean_difference"} or \code{"standardised_mean_difference"}
#'   (continuous); \code{"risk_difference"}, \code{"risk_ratio"} or
#'   \code{"odds_ratio"} (binary); \code{"rate_difference"},
#'   \code{"rate_ratio"} or \code{"log_rate_ratio"} (count); or
#'   \code{"hazard_ratio"} or \code{"log_hazard_ratio"} (time-to-event).
#' @param alternative Direction of the alternative hypothesis,
#'   \code{"greater"} (default) or \code{"less"}.
#' @param delta_null Null-hypothesis effect-size value (\eqn{\Delta_{H_0}} in
#'   the manuscript). Defaults to the no-effect value of the chosen
#'   \code{effect_scale}: 0 for the difference scales
#'   (\code{"mean_difference"}, \code{"standardised_mean_difference"},
#'   \code{"risk_difference"}, \code{"rate_difference"}) and for the log
#'   scales (\code{"log_rate_ratio"}, \code{"log_hazard_ratio"}), and 1 for
#'   the ratio scales (\code{"risk_ratio"}, \code{"odds_ratio"},
#'   \code{"rate_ratio"}, \code{"hazard_ratio"}).
#' @param sigma Common within-arm standard deviation of the simulated
#'   outcome data. Required for every continuous design, since it drives the
#'   Monte Carlo data generation, and additionally used by the known-variance
#'   normal prior.
#' @param allocation_ratio Treatment-to-control allocation ratio, with a
#'   uniform meaning across all four endpoints: continuous/binary split
#'   \code{n_per_look} into per-arm sample-size targets by this ratio; count
#'   splits \code{exposure_per_look} the same way; tte scales the
#'   treatment arm's Poisson recruitment rate (see \code{\link{set_accrual}})
#'   by this ratio relative to the control arm's, so the pooled
#'   \code{d_per_look}/\code{d_total} targets are reached under an unbalanced
#'   arm composition rather than a 1:1 split.
#' @return An object of class \code{"adabay_design"}.
#' @examples
#' set_design(endpoint = "binary",
#'            n_per_look = c(760, 1520, 2280, 3040, 3800),
#'            effect_scale = "risk_difference",
#'            alternative = "less")
#' @export
set_design <- function(endpoint,
                       n_per_look = NULL,
                       exposure_per_look = NULL,
                       d_total = NULL,
                       d_per_look = NULL,
                       effect_scale = NULL,
                       alternative = c("greater", "less"),
                       delta_null = NULL,
                       sigma = NULL,
                       allocation_ratio = 1) {
  endpoint    <- .match_choice(endpoint, .adabay_endpoints)
  alternative <- .match_choice(alternative, c("greater", "less"))

  .assert_numeric(allocation_ratio, len = 1L, positive = TRUE)

  if (is.null(effect_scale))
    effect_scale <- .adabay_scales[[endpoint]][1]
  effect_scale <- .match_choice(effect_scale, .adabay_scales[[endpoint]])

  if (is.null(delta_null))
    delta_null <- .adabay_null_values[[effect_scale]]
  .assert_numeric(delta_null, len = 1L)

  ## The number of looks K is not a separate argument: it is taken from the
  ## length of whichever pooled per-look vector the endpoint requires
  ## (n_per_look / exposure_per_look / d_per_look), so there is nothing for
  ## a stray 'n_looks' value to disagree with.
  schedule <- switch(endpoint,
                     continuous = ,
                     binary = {
                       if (is.null(n_per_look) || length(n_per_look) < 1L)
                         stop("Provide a non-empty 'n_per_look' for continuous and binary endpoints.",
                              call. = FALSE)
                       .assert_numeric(n_per_look, positive = TRUE, integer = TRUE)
                       if (any(diff(n_per_look) <= 0))
                         stop("'n_per_look' must be strictly increasing.",
                              call. = FALSE)
                       split <- .split_per_look(n_per_look, allocation_ratio,
                                                integer = TRUE,
                                                varname = "n_per_look")
                       list(n_c = split$c, n_t = split$t)
                     },
                     count = {
                       if (is.null(exposure_per_look) || length(exposure_per_look) < 1L)
                         stop("Provide a non-empty 'exposure_per_look' for the count endpoint.",
                              call. = FALSE)
                       .assert_numeric(exposure_per_look, positive = TRUE)
                       if (any(diff(exposure_per_look) <= 0))
                         stop("'exposure_per_look' must be strictly increasing.",
                              call. = FALSE)
                       split <- .split_per_look(exposure_per_look, allocation_ratio,
                                                integer = FALSE,
                                                varname = "exposure_per_look")
                       list(E_c = split$c, E_t = split$t)
                     },
                     tte = {
                       if (is.null(d_per_look) || length(d_per_look) < 1L)
                         stop("Provide a non-empty 'd_per_look' for the tte endpoint.",
                              call. = FALSE)
                       .assert_numeric(d_per_look, positive = TRUE, integer = TRUE)
                       if (any(diff(d_per_look) <= 0))
                         stop("'d_per_look' must be strictly increasing.",
                              call. = FALSE)
                       if (is.null(d_total)) {
                         d_total <- max(d_per_look)
                       } else {
                         .assert_numeric(d_total, len = 1L, positive = TRUE,
                                         integer = TRUE)
                         if (d_total < max(d_per_look))
                           stop("'d_total' must be at least max(d_per_look).",
                                call. = FALSE)
                       }
                       list(D = as.integer(d_per_look),
                            D_total = as.integer(d_total))
                     })
  n_looks <- switch(endpoint,
                    continuous = , binary = length(n_per_look),
                    count      = length(exposure_per_look),
                    tte        = length(d_per_look))

  if (!is.null(sigma)) .assert_numeric(sigma, len = 1L, positive = TRUE)

  structure(
    list(endpoint         = endpoint,
         n_looks          = n_looks,
         schedule         = schedule,
         effect_scale     = effect_scale,
         alternative      = alternative,
         delta_null       = delta_null,
         sigma            = sigma,
         allocation_ratio = allocation_ratio),
    class = "adabay_design")
}

## Split a pooled (both-arms) per-look schedule into per-arm sequences,
## exact by construction: arm_c[k] + arm_t[k] == pooled[k] at every look.
## For integer schedules (continuous/binary sample size) arm_c is rounded
## and arm_t takes the remainder, so the pooled total is matched exactly
## even though the individual arms are rounded; for continuous schedules
## (count exposure) both are exact fractions of pooled, with no rounding.
## Errors if the derived per-arm sequences are not strictly increasing at
## every look -- e.g. consecutive pooled values too close together to
## guarantee each arm gains at least one whole unit at every look under the
## given allocation_ratio.
.split_per_look <- function(pooled, allocation_ratio, integer, varname) {
  if (integer) {
    ## Note: round() uses IEC 60559 round-half-to-even ("banker's rounding"),
    ## so a pooled value whose control-arm share is an exact half-integer at
    ## consecutive looks (e.g. n_per_look = c(5, 7, 9) at allocation_ratio = 1)
    ## can leave arm_c non-monotone and be rejected below even though the
    ## pooled increments are adequate. Realistic schedules use increments
    ## large enough that the half-integer boundary is not hit.
    arm_c <- as.integer(round(pooled / (1 + allocation_ratio)))
    arm_t <- as.integer(pooled) - arm_c
  } else {
    arm_c <- as.numeric(pooled) / (1 + allocation_ratio)
    arm_t <- as.numeric(pooled) - arm_c
  }
  if (any(arm_c <= 0) || any(arm_t <= 0))
    stop(sprintf(
      "'%s' too small to give both arms a positive size at every look at allocation_ratio = %g (a derived per-arm value rounds to zero).",
      varname, allocation_ratio), call. = FALSE)
  if (any(diff(arm_c) <= 0) || any(diff(arm_t) <= 0))
    stop(sprintf(
      "'%s' increments too small to guarantee a strictly increasing per-arm schedule at allocation_ratio = %g.",
      varname, allocation_ratio), call. = FALSE)
  list(c = arm_c, t = arm_t)
}
