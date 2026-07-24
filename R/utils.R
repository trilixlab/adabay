## Internal helpers shared across the package.
## These functions are not exported.

## Supported endpoint types.
.adabay_endpoints <- c("continuous", "binary", "count", "tte")

## Supported effect scales for each endpoint.
.adabay_scales <- list(
  continuous = c("mean_difference", "standardised_mean_difference"),
  binary     = c("risk_difference", "risk_ratio", "odds_ratio"),
  count      = c("rate_difference", "rate_ratio", "log_rate_ratio"),
  tte        = c("hazard_ratio", "log_hazard_ratio")
)

## Default null values for each effect scale.
.adabay_null_values <- list(
  mean_difference              = 0,
  standardised_mean_difference = 0,
  risk_difference              = 0,
  risk_ratio                   = 1,
  odds_ratio                   = 1,
  rate_difference              = 0,
  rate_ratio                   = 1,
  log_rate_ratio               = 0,
  hazard_ratio                 = 1,
  log_hazard_ratio             = 0
)

## Default conjugate kernel for each endpoint.
.adabay_kernels <- list(
  continuous = c("normal", "normal_inverse_gamma"),
  binary     = "beta",
  count      = "gamma",
  tte        = "gamma"
)

## Effect scales whose posterior tail probability is evaluated by fixed-node
## Gauss-Legendre quadrature rather than a closed-form CDF. On these scales the
## quadrature loses accuracy when a posterior shape parameter falls below 1,
## because it cannot resolve the integrable density-endpoint singularity. The
## relative-effect scales for count and time-to-event endpoints use an exact closed
## form (a beta CDF) and are unaffected at any shape.
.adabay_quadrature_scales <- c("risk_difference", "risk_ratio", "odds_ratio",
                               "rate_difference")

## Smallest prior shape parameter that governs an integrable density-endpoint
## singularity on a quadrature scale: (a, b) for beta kernels (either endpoint),
## a (shape) for gamma kernels (lower endpoint). Returns Inf when no relevant
## shape applies (e.g. continuous kernels).
.min_prior_shape <- function(prior) {
  shapes <- unlist(lapply(prior$arms, function(arm) {
    lapply(arm$components, function(comp) {
      if (identical(prior$endpoint, "binary"))
        c(comp$a, comp$b)   # beta: singularity at 0 if a<1, at 1 if b<1
      else if (!is.null(comp$a))
        comp$a              # gamma: singularity at 0 if shape<1
      else NA_real_
    })
  }), use.names = FALSE)
  shapes <- shapes[is.finite(shapes)]
  if (!length(shapes)) return(Inf)
  min(shapes)
}

## Warn once per entry-point call when a fixed-node quadrature scale is paired
## with a prior that can drive a posterior shape below 1, in which case the
## per-look tail-probability estimate may be inaccurate at looks with few
## events. Closed-form (relative-effect) scales never trigger this.
.warn_quadrature_shape <- function(design, prior) {
  scale <- design$effect_scale
  if (is.null(scale) || !(scale %in% .adabay_quadrature_scales))
    return(invisible())
  if (.min_prior_shape(prior) >= 1)
    return(invisible())
  alt <- if (identical(design$endpoint, "count"))
    paste0(" Use a relative-effect scale (\"rate_ratio\" or ",
           "\"log_rate_ratio\"), which uses an exact closed form, or a prior ",
           "with all shape parameters at least 1.")
  else
    " Use a prior whose shape parameters (a, b) are all at least 1."
  warning(sprintf(
    paste0("Effect scale '%s' evaluates posterior tail probabilities by ",
           "fixed-node Gauss-Legendre quadrature, which loses accuracy when a ",
           "posterior shape parameter falls below 1 (as can occur at early ",
           "looks under a prior shape < 1, e.g. a Jeffreys prior).%s"),
    scale, alt), call. = FALSE)
}

## Validate that x is one of choices and return it.
##
## Accepts the match.arg() convention: if x is a vector identical to choices
## (the default-argument case), the first element is selected. Otherwise x
## must be a single value matching one of choices.
.match_choice <- function(x, choices, name = deparse(substitute(x))) {
  if (length(x) > 1L && identical(as.character(x), as.character(choices))) {
    return(as.character(choices[1]))
  }
  if (length(x) != 1L || is.na(x))
    stop(sprintf("'%s' must be a single non-missing value.", name), call. = FALSE)
  match.arg(as.character(x), choices)
}

## Validate that x is a finite numeric vector of given length.
.assert_numeric <- function(x, name = deparse(substitute(x)),
                            len = NULL, positive = FALSE,
                            integer = FALSE) {
  if (!is.numeric(x) || any(!is.finite(x)))
    stop(sprintf("'%s' must be a finite numeric vector.", name), call. = FALSE)
  if (!is.null(len) && length(x) != len)
    stop(sprintf("'%s' must have length %d (got %d).",
                 name, len, length(x)), call. = FALSE)
  if (positive && any(x <= 0))
    stop(sprintf("'%s' must be strictly positive.", name), call. = FALSE)
  if (integer && any(x != round(x)))
    stop(sprintf("'%s' must be a vector of integers.", name), call. = FALSE)
  invisible(x)
}

## Cross-platform parallel apply with reproducible RNG streams.
##
## Uses parallel::mclapply on Unix when cores > 1, and parLapply on Windows.
## Pass cores = 1 to force sequential execution.
##
## The L'Ecuyer-CMRG generator is required for independent per-worker streams
## under mclapply. We always switch the calling thread to it for the duration
## of any parallel call (POSIX) and restore the user's original RNGkind on
## exit. When 'seed' is NULL we still switch generators -- otherwise mclapply
## forks would inherit a shared Mersenne-Twister state and produce identical
## per-worker streams -- but we seed the parent stream from R's current state
## rather than from a fixed value, so the run remains non-deterministic.
.par_apply <- function(X, FUN, cores = 1L, seed = NULL, ...) {
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) cores <- 1L
  if (cores == 1L) {
    if (!is.null(seed)) set.seed(seed)
    return(lapply(X, FUN, ...))
  }
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    if (!is.null(seed))
      parallel::clusterSetRNGStream(cl, iseed = seed)
    parallel::parLapply(cl, X, FUN, ...)
  } else {
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old_kind[1L]), add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    parallel::mclapply(X, FUN, ..., mc.cores = cores)
  }
}

## Apply the alternative direction to a posterior tail probability.
##
## For "greater" the rule is P(Delta > e) > p; for "less" it is P(Delta < e) > p.
## .tail_prob_eff() always returns the probability of the side that triggers
## efficacy under the user-supplied alternative.
##
## prob_greater is P(Delta > e) for the marginal posterior.
.tail_prob_eff <- function(prob_greater, alternative) {
  if (alternative == "greater") prob_greater else 1 - prob_greater
}

.tail_prob_fut <- function(prob_greater, alternative) {
  ## Futility uses the opposite side from efficacy.
  if (alternative == "greater") 1 - prob_greater else prob_greater
}

## Normalise effect = list(...) into a named list with elements psi_c and psi_t,
## using the per-endpoint convention from the manuscript:
##   continuous: mu_c, mu_t        -> psi_c = mu_c,     psi_t = mu_t
##   binary:     theta_c, theta_t  -> psi_c = theta_c,  psi_t = theta_t
##   count, tte: lambda_c, lambda_t -> psi_c = lambda_c, psi_t = lambda_t
.normalise_effect <- function(effect, endpoint) {
  if (!is.list(effect) || is.null(names(effect)))
    stop("'effect' must be a named list.", call. = FALSE)
  pick <- function(c_keys, t_keys) {
    cn <- intersect(c_keys, names(effect))
    tn <- intersect(t_keys, names(effect))
    if (length(cn) == 0L || length(tn) == 0L)
      stop(sprintf("'effect' must supply control and treatment values; expected one of %s and %s.",
                   paste0("'", c_keys, "'", collapse = "/"),
                   paste0("'", t_keys, "'", collapse = "/")),
           call. = FALSE)
    list(psi_c = effect[[cn[1]]], psi_t = effect[[tn[1]]])
  }
  out <- switch(endpoint,
         continuous = pick(c("mu_c", "psi_c"),         c("mu_t", "psi_t")),
         binary     = pick(c("theta_c", "p_c", "psi_c"),  c("theta_t", "p_t", "psi_t")),
         count      = pick(c("lambda_c", "psi_c"),     c("lambda_t", "psi_t")),
         tte        = pick(c("lambda_c", "psi_c"),     c("lambda_t", "psi_t")),
         stop("Unknown endpoint", call. = FALSE))
  ## Domain-validate the data-generating values so invalid parameters fail fast
  ## with an informative error rather than producing NA or silently wrong
  ## operating characteristics downstream (e.g. via rbinom/rpois NAs).
  if (!is.numeric(out$psi_c) || !is.numeric(out$psi_t) ||
      length(out$psi_c) != 1L || length(out$psi_t) != 1L ||
      !is.finite(out$psi_c) || !is.finite(out$psi_t))
    stop("'effect' control and treatment values must each be a single finite number.",
         call. = FALSE)
  if (endpoint == "binary" &&
      (out$psi_c <= 0 || out$psi_c >= 1 || out$psi_t <= 0 || out$psi_t >= 1))
    stop("Binary 'effect' response rates (theta_c, theta_t) must lie strictly in (0, 1).",
         call. = FALSE)
  if (endpoint %in% c("count", "tte") &&
      (out$psi_c <= 0 || out$psi_t <= 0))
    stop("Count/tte 'effect' rates (lambda_c, lambda_t) must be strictly positive.",
         call. = FALSE)
  out
}
