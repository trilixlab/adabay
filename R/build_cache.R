#' Build a precomputation cache for grid evaluation
#'
#' Implements the precomputation strategy of Section 6.1 of the manuscript.
#' Simulates \code{n_trials} virtual trials once at the union of all
#' candidate look times, computes per-look posterior tail probabilities at
#' all candidate effect-size thresholds, and stores the result. Subsequent
#' calls to [evaluate_design()] or [calibrate_design()] then evaluate any
#' candidate \eqn{(K, \boldsymbol{n}, p_{k,u}, q_{k,v})} configuration
#' without further simulation cost.
#'
#' @param design A \code{adabay_design} object whose schedule represents the
#'   union of all candidate look times.
#' @param prior A \code{adabay_prior} object.
#' @param effect Data-generating values, as in [evaluate_design()].
#' @param accrual Optional [set_accrual()] specification. **Required** for
#'   tte designs, which are always staggered: the recruitment rate fixes
#'   the Poisson enrolment process whose arrival times drive the calendar-time
#'   exposure cutoff. The tte endpoint currently supports only `model = "poisson"`.
#'   When supplied (any endpoint) it is stored in the cache so that
#'   [evaluate_design()] on the cache reports the expected duration (exact,
#'   from simulated calendar times, for tte; analytic for the other
#'   endpoints, which have no simulated recruitment timeline).
#' @param n_trials Number of virtual trials.
#' @param cores Number of CPU cores. Defaults to 1.
#' @param seed Optional integer seed.
#' @param threshold_grid Named list with elements \code{efficacy} and
#'   \code{futility}, each a numeric vector of effect-size thresholds to
#'   cache. Defaults to \code{list(efficacy = design$delta_null, futility =
#'   design$delta_null)}.
#' @param cross_core_reproducible Logical, default `FALSE`. When `TRUE`,
#'   the simulator is forced to run in a single driver process (the user-
#'   supplied `cores` is ignored for the simulation pass) so that the
#'   per-trial RNG draws stored in the cache are bit-identical across
#'   worker counts. The flag belongs on the simulating call only: the
#'   cached path ([evaluate_design()] on an `adabay_cache`) and
#'   [calibrate_design()] re-use the stored draws and are already
#'   core-independent, so setting it there has no effect.
#'   See [evaluate_design()] for the cost discussion.
#' @return An object of class \code{"adabay_cache"}.
#' @examples
#' \donttest{
#'   des <- set_design(endpoint     = "binary",
#'                     n_per_look   = c(760, 1520, 2280),
#'                     effect_scale = "risk_difference",
#'                     alternative  = "less")
#'   pri <- set_prior(endpoint = "binary",
#'                    arms = list(c = list(family = "beta", a = 1, b = 1),
#'                                t = list(family = "beta", a = 1, b = 1)))
#'   cache <- build_cache(des, pri,
#'                        effect         = list(theta_c = 0.33, theta_t = 0.28),
#'                        n_trials       = 200L,
#'                        seed           = 1L,
#'                        threshold_grid = list(efficacy = 0, futility = 0))
#'   print(cache)
#' }
#' @export
build_cache <- function(design, prior, effect,
                        accrual = NULL,
                        n_trials = 10000L,
                        cores = 1L, seed = NULL,
                        threshold_grid = NULL,
                        cross_core_reproducible = FALSE) {
  if (!inherits(design, "adabay_design"))
    stop("'design' must be a adabay_design object.", call. = FALSE)
  if (!inherits(prior, "adabay_prior"))
    stop("'prior' must be a adabay_prior object.", call. = FALSE)
  if (design$endpoint != prior$endpoint)
    stop("Design and prior endpoints differ.", call. = FALSE)
  if (!is.null(accrual) && !inherits(accrual, "adabay_accrual"))
    stop("'accrual' must be a adabay_accrual object.", call. = FALSE)
  .check_tte_accrual(design, accrual)
  .warn_quadrature_shape(design, prior)
  if (is.null(threshold_grid)) {
    threshold_grid <- list(efficacy = design$delta_null,
                           futility = design$delta_null)
  }
  n_trials <- as.integer(n_trials)
  if (n_trials < 1L)
    stop("'n_trials' must be a positive integer.", call. = FALSE)

  thresholds <- sort(unique(c(threshold_grid$efficacy,
                              threshold_grid$futility)))
  if (length(thresholds) == 0L)
    stop("'threshold_grid' must include at least one threshold.",
         call. = FALSE)

  sim <- .simulate_data(design, effect, n_trials = n_trials,
                        seed = seed, cores = cores,
                        cross_core_reproducible = cross_core_reproducible,
                        accrual = accrual)

  n_looks <- design$n_looks
  tails <- vector("list", n_looks)
  for (k in seq_len(n_looks))
    tails[[k]] <- .compute_tail_probs_at_look(design, prior, sim, k,
                                              thresholds)

  structure(
    list(design     = design,
         prior      = prior,
         effect     = effect,
         accrual    = accrual,
         thresholds = thresholds,
         n_trials   = n_trials,
         tails      = tails,
         sim        = sim,
         seed       = seed),
    class = "adabay_cache")
}
