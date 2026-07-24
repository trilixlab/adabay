#' Set per-arm priors for a Bayesian group sequential design
#'
#' Constructs a prior specification matched to the endpoint of the design.
#' Each arm is given a (possibly mixture) conjugate prior. Pure conjugate
#' priors are passed in directly; non-conjugate priors are first approximated
#' by a finite mixture of conjugate components via [fit_mixture()] and then
#' passed in as the resulting object.
#'
#' Supported conjugate kernels per endpoint:
#'
#' \describe{
#'   \item{continuous}{\code{family = "normal"} with \code{mean} and
#'         \code{sd} (known variance), or
#'         \code{family = "normal_inverse_gamma"} with \code{nu, kappa,
#'         alpha, beta} (unknown variance).}
#'   \item{binary}{\code{family = "beta"} with \code{a} and \code{b}.}
#'   \item{count, tte}{\code{family = "gamma"} with \code{a} (shape) and
#'         \code{b} (rate).}
#' }
#'
#' Each per-arm specification can also be a multi-component mixture supplied
#' as a list with element \code{components} (list of conjugate components)
#' and \code{weights} (numeric vector summing to 1).
#'
#' On the effect scales evaluated by fixed-node Gauss--Legendre quadrature
#' (\code{"risk_difference"}, \code{"risk_ratio"} and \code{"odds_ratio"} for
#' binary endpoints, and \code{"rate_difference"} for count endpoints), the
#' posterior tail probability loses accuracy when a posterior shape parameter
#' falls below 1, as can happen at early looks under a prior shape below 1 (for
#' example a Jeffreys prior); [evaluate_design()] and [build_cache()] warn in
#' that case. Prefer a prior whose shape parameters are all at least 1, or, for
#' count endpoints, a relative-effect scale (\code{"rate_ratio"} or
#' \code{"log_rate_ratio"}), which uses an exact closed form at any shape.
#'
#' @param endpoint One of \code{"continuous"}, \code{"binary"}, \code{"count"}
#'   or \code{"tte"} (time-to-event).
#' @param arms A list with named elements \code{c} (control) and \code{t}
#'   (treatment), each a per-arm prior specification (see Details).
#' @return An object of class \code{"adabay_prior"}.
#' @examples
#' set_prior(endpoint = "binary",
#'           arms = list(c = list(family = "beta", a = 1, b = 1),
#'                       t = list(family = "beta", a = 1, b = 1)))
#' @export
set_prior <- function(endpoint, arms) {
  endpoint <- .match_choice(endpoint, .adabay_endpoints)
  if (!is.list(arms) || !all(c("c", "t") %in% names(arms)))
    stop("'arms' must be a list with named elements 'c' and 't'.",
         call. = FALSE)

  arm_c <- .normalise_arm_prior(arms$c, endpoint, arm_label = "c")
  arm_t <- .normalise_arm_prior(arms$t, endpoint, arm_label = "t")

  structure(
    list(endpoint = endpoint,
         arms = list(c = arm_c, t = arm_t)),
    class = "adabay_prior")
}

## Validate and canonicalise a per-arm prior specification.
.normalise_arm_prior <- function(spec, endpoint, arm_label) {
  if (!is.list(spec))
    stop(sprintf("Per-arm prior for '%s' must be a list.", arm_label),
         call. = FALSE)
  family <- spec$family
  if (is.null(family))
    stop(sprintf("Per-arm prior for '%s' must include 'family'.", arm_label),
         call. = FALSE)

  ok <- switch(endpoint,
               continuous = family %in% c("normal", "normal_inverse_gamma"),
               binary     = family %in% c("beta"),
               count      = family %in% c("gamma"),
               tte        = family %in% c("gamma"))
  if (!ok)
    stop(sprintf("Family '%s' is not a conjugate kernel for endpoint '%s'.",
                 family, endpoint), call. = FALSE)

  if (!is.null(spec$components) || !is.null(spec$weights)) {
    if (is.null(spec$components) || is.null(spec$weights))
      stop("A mixture prior must supply both 'components' and 'weights'.",
           call. = FALSE)
    if (length(spec$components) != length(spec$weights))
      stop("Length of 'components' and 'weights' must match.",
           call. = FALSE)
    .assert_numeric(spec$weights, positive = TRUE)
    if (abs(sum(spec$weights) - 1) > 1e-8)
      stop("Mixture weights must sum to 1.", call. = FALSE)
    components <- lapply(spec$components, .check_component, family = family)
    out <- list(family = family,
                weights = spec$weights,
                components = components)
    ## Preserve diagnostics and class on a fit_mixture() output so users can
    ## still inspect the per-quantile tail-error table and plot the overlay
    ## via the adabay_prior_arm methods.
    if (inherits(spec, "adabay_prior_arm")) {
      diag_attr <- attr(spec, "diagnostics")
      arm_attr  <- attr(spec, "arm")
      if (!is.null(diag_attr)) attr(out, "diagnostics") <- diag_attr
      if (!is.null(arm_attr))  attr(out, "arm")         <- arm_attr
      class(out) <- c("adabay_prior_arm", "list")
    }
    return(out)
  }
  comp <- .check_component(spec, family = family)
  list(family = family, weights = 1, components = list(comp))
}

.check_component <- function(comp, family) {
  if (!is.list(comp))
    stop("Each conjugate component must be a list.", call. = FALSE)
  switch(family,
         normal = {
           if (is.null(comp$mean) || is.null(comp$sd))
             stop("Normal component requires 'mean' and 'sd'.",
                  call. = FALSE)
           .assert_numeric(comp$mean, len = 1L)
           .assert_numeric(comp$sd, len = 1L, positive = TRUE)
           list(mean = comp$mean, sd = comp$sd)
         },
         normal_inverse_gamma = {
           need <- c("nu", "kappa", "alpha", "beta")
           miss <- setdiff(need, names(comp))
           if (length(miss))
             stop(sprintf("Normal--inverse-gamma component requires %s.",
                          paste0("'", need, "'", collapse = ", ")),
                  call. = FALSE)
           .assert_numeric(comp$nu,    len = 1L)
           .assert_numeric(comp$kappa, len = 1L, positive = TRUE)
           .assert_numeric(comp$alpha, len = 1L, positive = TRUE)
           .assert_numeric(comp$beta,  len = 1L, positive = TRUE)
           list(nu = comp$nu, kappa = comp$kappa,
                alpha = comp$alpha, beta = comp$beta)
         },
         beta = {
           if (is.null(comp$a) || is.null(comp$b))
             stop("Beta component requires 'a' and 'b'.", call. = FALSE)
           .assert_numeric(comp$a, len = 1L, positive = TRUE)
           .assert_numeric(comp$b, len = 1L, positive = TRUE)
           list(a = comp$a, b = comp$b)
         },
         gamma = {
           if (is.null(comp$a) || is.null(comp$b))
             stop("Gamma component requires 'a' (shape) and 'b' (rate).",
                  call. = FALSE)
           .assert_numeric(comp$a, len = 1L, positive = TRUE)
           .assert_numeric(comp$b, len = 1L, positive = TRUE)
           list(a = comp$a, b = comp$b)
         })
}
