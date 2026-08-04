## adabay -- trial-level simulation engine.
##
## This file contains the internal trial-simulation helpers that back the
## public `evaluate_design()` (fused path) and `build_cache()` entry points.
## It exports nothing on its own.

## Solve E{D_c(t)} + E{D_t(t)} = d for the calendar time t under independent
## per-arm homogeneous Poisson recruitment (rates rate_c/rate_t, allowing an
## unbalanced allocation_ratio) and per-arm exponential hazards
## (lambda_c/lambda_t), using the per-arm expected-event-count relation
## summed across arms:
##   E{D_a(t)} = rate_a t - (rate_a / lambda_a) (1 - e^{-lambda_a t}).
## Used only to size each arm's enrolment pool in .sim_tte() with a
## safety margin -- the expected trial duration reported to the user is
## computed directly from simulated calendar times (see
## .expected_duration()), not from this deterministic approximation.
## Returns NA when no positive root exists (e.g. a non-positive hazard).
.tte_tstar <- function(rate_c, rate_t, lambda_c, lambda_t, d) {
  if (!is.finite(lambda_c) || lambda_c <= 0 ||
      !is.finite(lambda_t) || lambda_t <= 0) return(NA_real_)
  f <- function(t) {
    (rate_c * t - rate_c / lambda_c * (1 - exp(-lambda_c * t))) +
      (rate_t * t - rate_t / lambda_t * (1 - exp(-lambda_t * t))) - d
  }
  tryCatch(stats::uniroot(f, c(1e-6, 1e6))$root,
           error = function(e) NA_real_)
}

## Compute the expected calendar duration E(T) under the supplied accrual
## specification. For continuous, binary and count there is no simulated
## recruitment timeline (these endpoints draw sufficient statistics directly
## at fixed cumulative sample-size/exposure targets, see .sim_continuous() /
## .sim_binary() / .sim_count()), so E(T) is a closed-form approximation for
## homogeneous Poisson recruitment at the design's pooled rate, combined with
## the simulated per-look stopping probabilities. accrual$rate is always the
## pooled (both-arms) rate (see set_accrual()), matching how n_per_look /
## exposure_per_look are pooled targets in set_design() -- so no per-arm
## split is needed here, the pooled sample size / exposure at each look is
## simply divided by the pooled rate directly.
##
## For tte, the per-trial calendar time at every look is already
## computed during simulation (sim_tte_cpp(); see .sim_tte()), so
## E(T) is instead the exact Monte Carlo mean of the calendar time realised
## at each trial's own stopping look -- no analytic approximation and no
## extra simulation cost. (Solving the deterministic mean-field equation
## E{D(t)} = D_k here instead would be systematically biased relative to the
## true expected hitting time of the underlying stochastic process, with the
## bias growing with the target event count.)
##
## Returns NA for unsupported configurations.
.expected_duration <- function(design, oc, accrual, sim = NULL) {
  if (accrual$model != "poisson") return(NA_real_)
  rate <- accrual$rate
  if (!is.numeric(rate) || length(rate) != 1L) return(NA_real_)
  follow_up <- if (is.numeric(accrual$follow_up)) accrual$follow_up else 0
  endpoint <- design$endpoint
  if (endpoint %in% c("continuous", "binary")) {
    n_at_look <- design$schedule$n_c + design$schedule$n_t
    t_at_look <- n_at_look / rate + follow_up
  } else if (endpoint == "count") {
    e_at_look <- design$schedule$E_c + design$schedule$E_t
    t_at_look <- e_at_look / rate
  } else if (endpoint == "tte") {
    if (is.null(sim) || is.null(sim$C_look)) return(NA_real_)
    return(mean(sim$C_look[cbind(seq_len(nrow(sim$C_look)), oc$tau)]))
  } else {
    return(NA_real_)
  }
  ## P(tau = k) = alpha_at_look_binding + fut_prob_at_look in either flag,
  ## because the simulator always stops trials at futility.
  p_tau <- oc$alpha_at_look_binding + oc$fut_prob_at_look
  sum(t_at_look * p_tau)
}

## Validate top-level inputs.
.check_inputs <- function(design, prior, decision) {
  if (!inherits(design, "adabay_design"))
    stop("'design' must be a adabay_design object.", call. = FALSE)
  if (!inherits(prior, "adabay_prior"))
    stop("'prior' must be a adabay_prior object.", call. = FALSE)
  if (!inherits(decision, "adabay_decision"))
    stop("'decision' must be a adabay_decision object.", call. = FALSE)
  if (design$endpoint != prior$endpoint)
    stop(sprintf("Design endpoint '%s' does not match prior endpoint '%s'.",
                 design$endpoint, prior$endpoint), call. = FALSE)
}

## Time-to-event (tte) designs are always staggered: a homogeneous-Poisson
## recruitment specification is mandatory, because the arrival times drive
## the calendar-time exposure cutoff. Piecewise / callback recruitment with
## a tte endpoint is not yet supported and errors here rather than being
## silently ignored. No-op for the other endpoints.
.check_tte_accrual <- function(design, accrual) {
  if (design$endpoint != "tte") return(invisible(NULL))
  if (is.null(accrual))
    stop("Time-to-event designs require an accrual specification; supply ",
         "accrual = set_accrual(model = \"poisson\", rate = ...).",
         call. = FALSE)
  if (accrual$model != "poisson")
    stop("Time-to-event designs currently support only Poisson recruitment; ",
         "supply accrual = set_accrual(model = \"poisson\", rate = ...).",
         call. = FALSE)
  invisible(NULL)
}

## Simulate cumulative per-look summary statistics under data-generating
## values, parallelised across cores by chunking the n_trials.
##
## When cross_core_reproducible = TRUE, the simulator is forced to run in a
## single driver process (cores = 1L is used internally regardless of the
## user-supplied 'cores'), so the trial-level RNG draws are bit-identical at
## any worker count. The posterior and aggregation passes are already core-
## independent. The wall-clock cost of cross_core_reproducible = TRUE is
## exactly the loss of simulator parallelism (small for binary/count/
## tte, where the simulator is a small fraction of total wall-clock;
## moderate for continuous, where the rnorm pass is larger). The default is
## FALSE, so the simulator uses all requested cores; its output is then
## bit-identical only at a fixed (seed, cores) pair.
.simulate_data <- function(design, effect, n_trials,
                           seed = NULL, cores = 1L,
                           cross_core_reproducible = FALSE,
                           accrual = NULL) {
  endpoint <- design$endpoint
  psi      <- .normalise_effect(effect, endpoint)
  if (isTRUE(cross_core_reproducible)) cores <- 1L
  if (cores <= 1L) {
    if (!is.null(seed)) set.seed(seed)
    return(switch(endpoint,
                  continuous = .sim_continuous(design, psi, n_trials),
                  binary     = .sim_binary(design, psi, n_trials),
                  count      = .sim_count(design, psi, n_trials),
                  tte        = .sim_tte(design, psi, n_trials, accrual)))
  }
  chunks  <- .chunk_indices(n_trials, cores)
  results <- .par_apply(seq_along(chunks), function(b) {
    chunk_size <- length(chunks[[b]])
    switch(endpoint,
           continuous = .sim_continuous(design, psi, chunk_size),
           binary     = .sim_binary(design, psi, chunk_size),
           count      = .sim_count(design, psi, chunk_size),
           tte        = .sim_tte(design, psi, chunk_size, accrual))
  }, cores = cores, seed = seed)
  .merge_sims(results, endpoint)
}

.chunk_indices <- function(n_trials, cores) {
  cores <- max(1L, min(as.integer(cores), n_trials))
  split(seq_len(n_trials),
        cut(seq_len(n_trials), cores, labels = FALSE))
}

.merge_sims <- function(results, endpoint) {
  out <- results[[1]]
  if (endpoint == "continuous") {
    for (nm in c("ybar_c", "ybar_t", "S2_c", "S2_t"))
      out[[nm]] <- do.call(rbind, lapply(results, `[[`, nm))
  } else if (endpoint == "binary") {
    for (nm in c("S_c", "S_t"))
      out[[nm]] <- do.call(rbind, lapply(results, `[[`, nm))
  } else if (endpoint == "count") {
    for (nm in c("S_c", "S_t"))
      out[[nm]] <- do.call(rbind, lapply(results, `[[`, nm))
  } else if (endpoint == "tte") {
    for (nm in c("D_c", "D_t", "E_c", "E_t", "C_look"))
      out[[nm]] <- do.call(rbind, lapply(results, `[[`, nm))
  }
  out
}

.sim_continuous <- function(design, psi, n_trials) {
  n_looks <- design$n_looks
  n_c     <- design$schedule$n_c
  n_t     <- design$schedule$n_t
  sigma   <- design$sigma
  if (is.null(sigma))
    stop("For continuous endpoints, supply 'sigma' to set_design().",
         call. = FALSE)
  inc_c <- diff(c(0, n_c))
  inc_t <- diff(c(0, n_t))
  cum_sum_c <- numeric(n_trials); cum_sum_t <- numeric(n_trials)
  cum_ss_c  <- numeric(n_trials); cum_ss_t  <- numeric(n_trials)
  ybar_c <- ybar_t <- S2_c <- S2_t <-
    matrix(0, nrow = n_trials, ncol = n_looks)
  for (k in seq_len(n_looks)) {
    if (inc_c[k] > 0) {
      x <- matrix(stats::rnorm(n_trials * inc_c[k],
                               mean = psi$psi_c, sd = sigma),
                  nrow = n_trials)
      cum_sum_c <- cum_sum_c + .rowSums(x, n_trials, inc_c[k])
      cum_ss_c  <- cum_ss_c  + .rowSums(x^2, n_trials, inc_c[k])
    }
    if (inc_t[k] > 0) {
      x <- matrix(stats::rnorm(n_trials * inc_t[k],
                               mean = psi$psi_t, sd = sigma),
                  nrow = n_trials)
      cum_sum_t <- cum_sum_t + .rowSums(x, n_trials, inc_t[k])
      cum_ss_t  <- cum_ss_t  + .rowSums(x^2, n_trials, inc_t[k])
    }
    ybar_c[, k] <- cum_sum_c / n_c[k]
    ybar_t[, k] <- cum_sum_t / n_t[k]
    S2_c[, k]   <- cum_ss_c - n_c[k] * ybar_c[, k]^2
    S2_t[, k]   <- cum_ss_t - n_t[k] * ybar_t[, k]^2
  }
  list(endpoint = "continuous",
       ybar_c = ybar_c, ybar_t = ybar_t,
       S2_c = S2_c, S2_t = S2_t,
       n_c = n_c, n_t = n_t, psi = psi)
}

.rowSums <- function(m, n, p) base::.rowSums(m, n, p, na.rm = FALSE)

.sim_binary <- function(design, psi, n_trials) {
  n_looks <- design$n_looks
  n_c     <- design$schedule$n_c
  n_t     <- design$schedule$n_t
  inc_c <- diff(c(0, n_c))
  inc_t <- diff(c(0, n_t))
  cum_S_c <- integer(n_trials); cum_S_t <- integer(n_trials)
  S_c <- S_t <- matrix(0L, nrow = n_trials, ncol = n_looks)
  for (k in seq_len(n_looks)) {
    if (inc_c[k] > 0)
      cum_S_c <- cum_S_c + stats::rbinom(n_trials, inc_c[k], psi$psi_c)
    if (inc_t[k] > 0)
      cum_S_t <- cum_S_t + stats::rbinom(n_trials, inc_t[k], psi$psi_t)
    S_c[, k] <- cum_S_c
    S_t[, k] <- cum_S_t
  }
  list(endpoint = "binary",
       S_c = S_c, S_t = S_t,
       n_c = n_c, n_t = n_t, psi = psi)
}

.sim_count <- function(design, psi, n_trials) {
  n_looks <- design$n_looks
  E_c <- design$schedule$E_c
  E_t <- design$schedule$E_t
  inc_E_c <- diff(c(0, E_c))
  inc_E_t <- diff(c(0, E_t))
  cum_S_c <- integer(n_trials); cum_S_t <- integer(n_trials)
  S_c <- S_t <- matrix(0L, nrow = n_trials, ncol = n_looks)
  for (k in seq_len(n_looks)) {
    if (inc_E_c[k] > 0)
      cum_S_c <- cum_S_c + stats::rpois(n_trials, inc_E_c[k] * psi$psi_c)
    if (inc_E_t[k] > 0)
      cum_S_t <- cum_S_t + stats::rpois(n_trials, inc_E_t[k] * psi$psi_t)
    S_c[, k] <- cum_S_c
    S_t[, k] <- cum_S_t
  }
  list(endpoint = "count",
       S_c = S_c, S_t = S_t,
       E_c = E_c, E_t = E_t, psi = psi)
}

## Simulate the tte sufficient statistics (D_a, E_a) under staggered
## homogeneous-Poisson enrolment. accrual$rate is the pooled (both-arms)
## recruitment rate -- matching how n_per_look/exposure_per_look are pooled
## targets in set_design() -- split into per-arm rates by
## design$allocation_ratio exactly as those are split into per-arm sample
## size/exposure targets, so allocation_ratio has a uniform meaning across
## all four endpoints.
##
## Recruitment is mandatory for tte: the per-arm rate fixes the
## enrolment process whose arrival times feed the calendar-time exposure
## cutoff (see sim_tte_cpp). Each arm's enrolment pool is sized
## separately from the analytic calendar time at which the (possibly
## asymmetric) pooled E{D(t)} reaches the final look (the per-arm
## expected-event-count relation of .tte_tstar(), summed per arm), inflated
## by an 8-SD margin so that the last-enrolled subject in
## each arm almost surely arrives after the final look -- this sizing is an
## internal safety margin only; it is not used to report the expected trial
## duration (see .expected_duration()).
.sim_tte <- function(design, psi, n_trials, accrual) {
  if (is.null(accrual) || !inherits(accrual, "adabay_accrual") ||
      accrual$model != "poisson")
    stop("Time-to-event simulation requires a Poisson accrual specification ",
         "(set_accrual(model = \"poisson\", rate = ...)).", call. = FALSE)
  d_per_look   <- design$schedule$D
  d_total      <- design$schedule$D_total
  alloc        <- design$allocation_ratio %||% 1
  rate_c       <- accrual$rate / (1 + alloc)
  rate_t       <- rate_c * alloc
  t_star       <- .tte_tstar(rate_c, rate_t, psi$psi_c, psi$psi_t, d_total)
  mu_c         <- if (is.finite(t_star)) rate_c * t_star else 0
  mu_t         <- if (is.finite(t_star)) rate_t * t_star else 0
  n_pool_c     <- max(ceiling(mu_c + 8 * sqrt(mu_c) + 50), 4L * d_total, 100L)
  n_pool_t     <- max(ceiling(mu_t + 8 * sqrt(mu_t) + 50), 4L * d_total, 100L)
  out <- sim_tte_cpp(R = n_trials,
                     n_pool_c = n_pool_c, n_pool_t = n_pool_t,
                     lambda_c = psi$psi_c, lambda_t = psi$psi_t,
                     rate_c = rate_c, rate_t = rate_t,
                     D_per_look = as.integer(d_per_look))
  list(endpoint = "tte",
       D_c = out$D_c, D_t = out$D_t,
       E_c = out$E_c, E_t = out$E_t,
       C_look = out$C,
       d_per_look = d_per_look, psi = psi)
}

.compute_tail_probs <- function(design, prior, sim, decision, cores = 1L) {
  n_looks    <- design$n_looks
  thresholds <- .threshold_grid(decision, n_looks)
  out <- vector("list", n_looks)
  for (k in seq_len(n_looks))
    out[[k]] <- .compute_tail_probs_at_look(design, prior, sim, k,
                                            thresholds[[k]])
  attr(out, "thresholds") <- thresholds
  out
}

.threshold_grid <- function(decision, n_looks) {
  lapply(seq_len(n_looks), function(k) {
    eff <- .criteria_at_look(list(criteria = decision$efficacy$criteria,
                                  per_look = decision$efficacy$per_look), k)
    fut <- if (length(decision$futility$criteria))
      .criteria_at_look(list(criteria = decision$futility$criteria,
                             per_look = decision$futility$per_look), k)
      else list()
    e_eff <- vapply(eff, function(c) c$threshold_effect, numeric(1))
    e_fut <- if (length(fut))
      vapply(fut, function(c) c$threshold_effect, numeric(1))
    else
      numeric(0)
    sort(unique(c(e_eff, e_fut)))
  })
}

.compute_tail_probs_at_look <- function(design, prior, sim, k, thresholds) {
  endpoint <- design$endpoint
  if (length(thresholds) == 0L)
    return(matrix(NA_real_, nrow = .R_of(sim), ncol = 0L))
  switch(endpoint,
         continuous = .tail_prob_continuous(prior, design$sigma,
                                            scale = design$effect_scale,
                                            n_c = sim$n_c[k], n_t = sim$n_t[k],
                                            ybar_c = sim$ybar_c[, k],
                                            ybar_t = sim$ybar_t[, k],
                                            S2_c = sim$S2_c[, k],
                                            S2_t = sim$S2_t[, k],
                                            thresholds = thresholds),
         binary     = .tail_prob_binary(prior, scale = design$effect_scale,
                                        n_c = sim$n_c[k], n_t = sim$n_t[k],
                                        S_c = sim$S_c[, k], S_t = sim$S_t[, k],
                                        thresholds = thresholds),
         count      = .tail_prob_gamma(prior, scale = design$effect_scale,
                                       E_c = sim$E_c[k], E_t = sim$E_t[k],
                                       S_c = sim$S_c[, k], S_t = sim$S_t[, k],
                                       thresholds = thresholds),
         tte        = .tail_prob_gamma(prior, scale = design$effect_scale,
                                       E_c = sim$E_c[, k], E_t = sim$E_t[, k],
                                       S_c = sim$D_c[, k], S_t = sim$D_t[, k],
                                       thresholds = thresholds))
}

.R_of <- function(sim) {
  if (sim$endpoint == "continuous") return(nrow(sim$ybar_c))
  if (sim$endpoint == "binary")     return(nrow(sim$S_c))
  if (sim$endpoint == "count")      return(nrow(sim$S_c))
  if (sim$endpoint == "tte")        return(nrow(sim$D_c))
  stop("Unknown endpoint", call. = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## Apply decision rules and aggregate operating characteristics.
##
## The simulator always conducts the trial under the assumption that
## crossing futility terminates it (binding behaviour). The actual stopping
## time \code{tau}, the trial result \code{C}, the expected sample size
## \code{expected_sample_size}, the per-look stopping probabilities
## \code{eff_prob_at_look} / \code{fut_prob_at_look}, and the futility probability
## \code{fut_prob} therefore reflect this binding conduct irrespective of the
## \code{futility_binding} flag, and \code{eff_prob_at_look + fut_prob_at_look}
## partitions the sample space (sum = 1).
##
## The flag changes only how the overall scalar type I error rate
## (\code{alpha}) is reported:
## (i) binding     -> alpha = sum(alpha_at_look_binding)
##                  = probability that the trial crosses efficacy before any
##                    futility crossing under the data-generating distribution.
## (ii) non-binding -> alpha = sum(alpha_at_look_nonbinding)
##                  = probability that the trial would have crossed efficacy
##                    at some look had the futility rule been ignored
##                    (= type I error of the design with no futility rule),
##                    which is the standard regulatory convention for type I
##                    error control. Under non-binding,
##                    alpha != sum(eff_prob_at_look) by design.
##
## To support both reports without re-simulating, the aggregator records the
## per-look efficacy crossings under two separate active-set policies:
## (i) binding (active until first eff or fut crossing), used for
## \code{alpha_at_look_binding}, \code{eff_prob_at_look}, and \code{tau} / \code{C};
## (ii) eff-only (active until first eff crossing, ignoring fut), used for
## \code{alpha_at_look_nonbinding} and feeding the non-binding overall alpha.
.aggregate_oc <- function(design, decision, tails, sim) {
  n_looks    <- design$n_looks
  n_trials   <- .R_of(sim)
  thresholds <- attr(tails, "thresholds")
  binding    <- isTRUE(decision$futility_binding)

  tau    <- rep(n_looks, n_trials)
  C      <- integer(n_trials)
  active_b <- rep(TRUE, n_trials)  # for binding: stops on eff or fut
  active_e <- rep(TRUE, n_trials)  # for non-binding alpha: stops on eff only
  per_look_eff_b <- per_look_fut_b <- integer(n_looks)
  per_look_eff_n <- integer(n_looks)

  for (k in seq_len(n_looks)) {
    eff_crit <- .criteria_at_look(list(criteria = decision$efficacy$criteria,
                                       per_look = decision$efficacy$per_look),
                                  k)
    fut_crit <- if (length(decision$futility$criteria))
      .criteria_at_look(list(criteria = decision$futility$criteria,
                             per_look = decision$futility$per_look), k)
      else list()
    e_grid   <- thresholds[[k]]
    tail_mat <- tails[[k]]

    eff_fire <- rep(TRUE, n_trials)
    for (crit in eff_crit) {
      idx <- match(crit$threshold_effect, e_grid)
      pg  <- tail_mat[, idx]
      eff_fire <- eff_fire &
        (.tail_prob_eff(pg, design$alternative) > crit$threshold_prob)
    }
    fut_fire <- if (length(fut_crit) == 0L) {
      rep(FALSE, n_trials)
    } else {
      ff <- rep(TRUE, n_trials)
      for (crit in fut_crit) {
        idx <- match(crit$threshold_effect, e_grid)
        pg  <- tail_mat[, idx]
        ff <- ff &
          (.tail_prob_fut(pg, design$alternative) > crit$threshold_prob)
      }
      ff
    }
    if (k == n_looks) fut_fire <- rep(FALSE, n_trials)

    ## Binding active-set update (governs trial conduct, E(N), futility rate).
    eff_now_b <- active_b & eff_fire
    if (k < n_looks) {
      fut_now_b <- active_b & fut_fire & !eff_now_b
      per_look_eff_b[k] <- sum(eff_now_b)
      per_look_fut_b[k] <- sum(fut_now_b)
      tau[eff_now_b | fut_now_b] <- k
      C[eff_now_b] <-  1L
      C[fut_now_b] <- -1L
      active_b[eff_now_b | fut_now_b] <- FALSE
    } else {
      per_look_eff_b[k] <- sum(active_b & eff_fire)
      per_look_fut_b[k] <- sum(active_b & !eff_fire)
      tau[active_b] <- n_looks
      C[active_b &  eff_fire] <-  1L
      C[active_b & !eff_fire] <- -1L
      active_b[] <- FALSE
    }

    ## Non-binding active-set update (governs alpha as if no futility).
    eff_now_n <- active_e & eff_fire
    per_look_eff_n[k] <- sum(eff_now_n)
    active_e[eff_now_n] <- FALSE
  }

  eff_prob_at_look_b <- per_look_eff_b / n_trials
  eff_prob_at_look_n <- per_look_eff_n / n_trials
  fut_prob_at_look    <- per_look_fut_b / n_trials

  ## Per-look stopping probabilities ALWAYS use binding-conduct
  ## (actual-trial first crossings), so eff_prob_at_look + fut_prob_at_look = 1
  ## by partition. The overall scalar alpha switches on the binding flag:
  ##   binding     -> sum(eff_prob_at_look_b)  (= sum of eff_prob_at_look)
  ##   non-binding -> sum(eff_prob_at_look_n)  (the "as if futility absent"
  ##                                           type I error rate, the
  ##                                           standard regulatory
  ##                                           convention; differs from
  ##                                           sum(eff_prob_at_look) by the
  ##                                           probability of crossing
  ##                                           futility-then-efficacy).
  eff_prob_at_look <- eff_prob_at_look_b
  alpha_overall  <- if (binding) sum(eff_prob_at_look_b) else sum(eff_prob_at_look_n)

  sample_at_look <- switch(design$endpoint,
                           continuous = design$schedule$n_c + design$schedule$n_t,
                           binary     = design$schedule$n_c + design$schedule$n_t,
                           count      = design$schedule$E_c + design$schedule$E_t,
                           tte        = design$schedule$D)
  ## E(N) reflects actual trial conduct: trials always stop at futility.
  ## P(tau = k) = eff_prob_at_look_b[k] + fut_prob_at_look[k] for all k.
  e_n <- sum(sample_at_look * (eff_prob_at_look_b + fut_prob_at_look))

  structure(
    list(design = design, decision = decision,
         eff_prob_at_look          = eff_prob_at_look,
         fut_prob_at_look          = fut_prob_at_look,
         alpha_at_look_binding    = eff_prob_at_look_b,
         alpha_at_look_nonbinding = eff_prob_at_look_n,
         alpha    = alpha_overall,
         eff_prob = sum(eff_prob_at_look),
         fut_prob = sum(fut_prob_at_look),
         expected_sample_size = e_n,
         tau = tau, C = C),
    class = "adabay_oc")
}
