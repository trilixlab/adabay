#' Set per-look decision rules
#'
#' Defines the efficacy and futility stopping rules for a Bayesian group
#' sequential design. A decision rule is a list of one or more criteria, each
#' a list with elements \code{threshold_effect} (effect-size threshold
#' \eqn{e_{k,u}} or \eqn{f_{k,v}}) and \code{threshold_prob}
#' (posterior-probability threshold \eqn{p_{k,u}} or \eqn{q_{k,v}}). All
#' criteria must be simultaneously met for the rule to fire.
#'
#' Rules can vary per look by passing a list of lists with one element per
#' look (length \eqn{K}, the number of looks). Supplying a single list applies
#' the same rule at every look.
#'
#' @param efficacy List of efficacy criteria (each with elements
#'   \code{threshold_effect} and \code{threshold_prob}), or a list of such
#'   lists, one per look.
#' @param futility List of futility criteria, with the same structure as
#'   \code{efficacy}. Pass \code{NULL} or an empty list to disable futility
#'   stopping.
#' @param futility_binding Logical. The simulator always conducts the trial
#'   under binding behaviour (trials stop on either efficacy or futility).
#'   This flag changes only the reported type I error rate. When
#'   \code{TRUE}, the type I error rate counts only efficacy crossings that
#'   occur before any futility crossing under the data-generating
#'   distribution. When \code{FALSE} (the default), the type I error rate
#'   is reported as the rate that would be observed if the futility rule
#'   had been ignored, which is identical to a design with no futility rule
#'   at all. The expected sample size, the futility probability, the per-look
#'   stopping probabilities and the trial-level outcomes \code{tau} and
#'   \code{C} are unchanged by the flag.
#' @return An object of class \code{"adabay_decision"}.
#' @examples
#' set_decision(
#'   efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
#'   futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
#'   futility_binding = FALSE)
#' @export
set_decision <- function(efficacy,
                         futility = NULL,
                         futility_binding = FALSE) {
  if (!is.logical(futility_binding) || length(futility_binding) != 1L)
    stop("'futility_binding' must be a single logical.", call. = FALSE)

  efficacy_per_look <- .normalise_rule_list(efficacy, "efficacy")
  if (is.null(futility) ||
      (is.list(futility) && length(futility) == 0L)) {
    futility_per_look <- list()
  } else {
    futility_per_look <- .normalise_rule_list(futility, "futility")
  }

  structure(
    list(efficacy         = efficacy_per_look,
         futility         = futility_per_look,
         futility_binding = futility_binding),
    class = "adabay_decision")
}

.normalise_rule_list <- function(rule, name) {
  if (!is.list(rule))
    stop(sprintf("'%s' must be a list.", name), call. = FALSE)
  if (.is_criterion(rule)) rule <- list(rule)
  if (all(vapply(rule, .is_criterion, logical(1)))) {
    return(list(per_look = FALSE,
                criteria = lapply(rule, .check_criterion, name = name)))
  }
  per_look <- lapply(rule, function(crit_list) {
    if (!is.list(crit_list))
      stop(sprintf("Each per-look entry in '%s' must be a list of criteria.",
                   name), call. = FALSE)
    lapply(crit_list, .check_criterion, name = name)
  })
  list(per_look = TRUE, criteria = per_look)
}

.is_criterion <- function(x) {
  is.list(x) && !is.null(x$threshold_effect) && !is.null(x$threshold_prob)
}

.check_criterion <- function(crit, name) {
  if (!.is_criterion(crit))
    stop(sprintf("Each %s criterion must be a list with elements 'threshold_effect' and 'threshold_prob'.",
                 name), call. = FALSE)
  .assert_numeric(crit$threshold_effect, name = "threshold_effect", len = 1L)
  .assert_numeric(crit$threshold_prob, name = "threshold_prob", len = 1L)
  if (crit$threshold_prob <= 0 || crit$threshold_prob > 1)
    stop("'threshold_prob' must lie in (0, 1].", call. = FALSE)
  crit
}

.criteria_at_look <- function(rule, k) {
  if (length(rule$criteria) == 0L) return(list())
  if (isTRUE(rule$per_look)) {
    if (k > length(rule$criteria))
      stop(sprintf("Per-look rule list shorter than the number of looks (k = %d).", k),
           call. = FALSE)
    rule$criteria[[k]]
  } else {
    rule$criteria
  }
}
