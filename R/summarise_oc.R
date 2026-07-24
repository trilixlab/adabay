#' Summarise operating characteristics in tidy form
#'
#' Returns a single-row \code{data.frame} summarising the operating
#' characteristics of one \code{adabay_oc} object, or a multi-row
#' \code{data.frame} when supplied with a list of \code{adabay_oc} objects.
#'
#' \code{alpha} is the overall type I error rate (regulatory convention:
#' under non-binding futility, computed as if futility were absent).
#' \code{eff_prob} is the overall probability of stopping for efficacy
#' under actual (binding) trial conduct -- the power when the simulation
#' was run under the alternative effect, and a coincident value with
#' \code{alpha} under binding futility. \code{fut_prob} is the overall
#' probability of stopping for futility under actual binding conduct.
#'
#' @param x A \code{adabay_oc} object or a list of them.
#' @return A \code{data.frame} with columns \code{alpha} (overall type I error),
#'   \code{eff_prob} (overall efficacy probability), \code{fut_prob} (overall
#'   futility probability), \code{E_N} (expected sample size / exposure / events
#'   under binding trial conduct), \code{E_T} (expected calendar duration when
#'   an \code{\link{set_accrual}} specification is supplied, otherwise \code{NA}),
#'   \code{n_trials}, \code{n_looks}, and the list-columns
#'   \code{eff_prob_at_look} and \code{fut_prob_at_look} carrying the stage
#'   stopping probability vectors.
#' @examples
#' \donttest{
#'   des <- set_design(endpoint     = "binary",
#'                     n_per_look   = c(760, 1520, 2280),
#'                     effect_scale = "risk_difference",
#'                     alternative  = "less")
#'   pri <- set_prior(endpoint = "binary",
#'                    arms = list(c = list(family = "beta", a = 1, b = 1),
#'                                t = list(family = "beta", a = 1, b = 1)))
#'   dec <- set_decision(
#'     efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
#'     futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
#'     futility_binding = TRUE)
#'   oc <- evaluate_design(des, pri, dec,
#'                         effect   = list(theta_c = 0.33, theta_t = 0.28),
#'                         n_trials = 200L, seed = 1L)
#'   summarise_oc(oc)
#' }
#' @export
summarise_oc <- function(x) {
  if (inherits(x, "adabay_oc")) return(.summary_one(x))
  if (is.list(x) && all(vapply(x, inherits, logical(1), what = "adabay_oc")))
    return(do.call(rbind, lapply(x, .summary_one)))
  stop("'x' must be a adabay_oc object or a list of them.", call. = FALSE)
}

.summary_one <- function(oc) {
  data.frame(
    alpha    = oc$alpha,
    eff_prob = oc$eff_prob,
    fut_prob = oc$fut_prob,
    E_N      = oc$expected_sample_size,
    E_T      = oc$expected_duration %||% NA_real_,
    n_trials = oc$n_trials %||% NA_integer_,
    n_looks  = oc$design$n_looks,
    eff_prob_at_look = I(list(oc$eff_prob_at_look)),
    fut_prob_at_look = I(list(oc$fut_prob_at_look)),
    row.names = NULL
  )
}
