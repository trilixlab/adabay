## Tests for the staggered-Poisson enrolment tte simulator.
##
## The tte OC simulator draws Poisson enrolment (arrival) times A_{a,i},
## draws per-subject time-to-event T_{a,i} ~ Exp(lambda_a), forms the calendar
## event time X = A + T, and builds the per-arm sufficient statistic
## (D_a, E_a) = (events, exposure) at each event-driven look using a
## CALENDAR-time cutoff c_k (the D_k-th smallest X pooled across arms):
##
##   D_{a,k} = #{ i : X_{a,i} <= c_k }
##   E_{a,k} = sum_{ i : A_{a,i} <= c_k } min(T_{a,i}, c_k - A_{a,i}).
##
## These tests pin (i) the mandatory-Poisson-accrual contract, (ii) the
## structural event-count identity D_c + D_t == D_k, and (iii) the
## operational-time identity E[lambda_a E_a] = E[D_a] that any faithful
## exposure construction must satisfy (lambda_a E_a | D_a ~ Gamma(D_a, 1) by
## exponential memorylessness, independent of the recruitment rate).

surv_design <- function(d_total = 60L, d_per_look = c(20L, 40L, 60L)) {
  set_design(endpoint = "tte",
             d_total = d_total, d_per_look = d_per_look,
             effect_scale = "log_hazard_ratio", alternative = "less")
}
surv_prior <- function() {
  set_prior(endpoint = "tte",
            arms = list(c = list(family = "gamma", a = 1, b = 1),
                        t = list(family = "gamma", a = 1, b = 1)))
}
surv_decision <- function() {
  set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)))
}
surv_effect <- function() list(lambda_c = log(2) / 12, lambda_t = log(2) / 24)

test_that("tte evaluate_design requires an accrual specification", {
  expect_error(
    evaluate_design(surv_design(), surv_prior(), surv_decision(),
                    effect = surv_effect(),
                    n_trials = 50, cores = 1, seed = 1L),
    "accrual")
})

test_that("tte evaluate_design rejects non-Poisson accrual", {
  acc <- set_accrual(model = "piecewise", rate = c(10, 20),
                     breakpoints = 5)
  expect_error(
    evaluate_design(surv_design(), surv_prior(), surv_decision(),
                    effect = surv_effect(), accrual = acc,
                    n_trials = 50, cores = 1, seed = 1L),
    "Poisson")
})

test_that("tte build_cache requires a Poisson accrual specification", {
  expect_error(
    build_cache(surv_design(), surv_prior(), effect = surv_effect(),
                n_trials = 50, seed = 1L,
                threshold_grid = list(efficacy = 0, futility = 0)),
    "accrual")
  acc <- set_accrual(model = "callback",
                     callback = function(n, ...) sort(stats::runif(n)))
  expect_error(
    build_cache(surv_design(), surv_prior(), effect = surv_effect(),
                accrual = acc, n_trials = 50, seed = 1L,
                threshold_grid = list(efficacy = 0, futility = 0)),
    "Poisson")
})

test_that("tte OC populates a finite expected duration", {
  acc <- set_accrual(model = "poisson", rate = 20)
  oc <- evaluate_design(surv_design(), surv_prior(), surv_decision(),
                        effect = surv_effect(), accrual = acc,
                        n_trials = 200, cores = 1, seed = 20260512L)
  expect_true(is.finite(oc$expected_duration))
  expect_true(oc$expected_duration > 0)
})

test_that("cached tte path also returns a finite expected duration", {
  acc <- set_accrual(model = "poisson", rate = 20)
  cache <- build_cache(surv_design(), surv_prior(), effect = surv_effect(),
                       accrual = acc, n_trials = 200, seed = 20260512L,
                       threshold_grid = list(efficacy = 0, futility = 0))
  oc <- evaluate_design(cache, surv_decision())
  expect_true(is.finite(oc$expected_duration))
  expect_true(oc$expected_duration > 0)
})

test_that("tte expected_duration is the exact simulated mean, not an analytic approximation", {
  ## Recompute the expected duration independently from the same simulated
  ## trials (same seed => same RNG draws, since no RNG-consuming code runs
  ## between set.seed() and .sim_tte() -- see .simulate_data()) and
  ## check it matches oc$expected_duration exactly. This locks in that
  ## expected_duration is genuinely derived from sim_tte_cpp()'s
  ## per-trial calendar times (sim$C_look), not a deterministic mean-field
  ## formula.
  acc    <- set_accrual(model = "poisson", rate = 20)
  des    <- surv_design()
  effect <- surv_effect()
  seed   <- 20260512L
  n      <- 2000L

  oc <- evaluate_design(des, surv_prior(), surv_decision(),
                        effect = effect, accrual = acc,
                        n_trials = n, cores = 1, seed = seed)

  psi <- adabay:::.normalise_effect(effect, "tte")
  set.seed(seed)
  sim   <- adabay:::.sim_tte(des, psi, n_trials = n, accrual = acc)
  tails <- adabay:::.compute_tail_probs(des, surv_prior(), sim, surv_decision(), cores = 1)
  oc2   <- adabay:::.aggregate_oc(des, surv_decision(), tails, sim)

  expect_true(is.matrix(sim$C_look))
  expect_identical(dim(sim$C_look), c(n, des$n_looks))
  true_duration <- mean(sim$C_look[cbind(seq_len(n), oc2$tau)])
  expect_equal(oc$expected_duration, true_duration, tolerance = 1e-8)
})

test_that("tte allocation_ratio scales the treatment arm's recruitment rate", {
  ## Regression guard for the bug where allocation_ratio was silently
  ## ignored for the tte endpoint (set_design() stored it but
  ## .sim_tte() never read it, so both arms always recruited at the
  ## accrual's single rate regardless of allocation_ratio).
  acc  <- set_accrual(model = "poisson", rate = 20)
  des1 <- surv_design()
  des2 <- set_design(endpoint = "tte", d_total = 60,
                     d_per_look = c(20, 40, 60),
                     effect_scale = "log_hazard_ratio", alternative = "less",
                     allocation_ratio = 2)
  expect_equal(des1$allocation_ratio, 1)
  expect_equal(des2$allocation_ratio, 2)

  psi <- list(psi_c = log(2) / 12, psi_t = log(2) / 24)
  n <- 20000L
  set.seed(1L)
  sim1 <- adabay:::.sim_tte(des1, psi, n_trials = n, accrual = acc)
  set.seed(1L)
  sim2 <- adabay:::.sim_tte(des2, psi, n_trials = n, accrual = acc)

  ## allocation_ratio must actually change the simulated data now (was a
  ## no-op before the fix).
  expect_false(identical(sim1$D_t, sim2$D_t))

  ## The pooled target event count at each look is invariant to allocation
  ## (D_c + D_t == d_per_look[k] regardless of the arm composition)...
  d_per_look <- des1$schedule$D
  for (k in seq_along(d_per_look)) {
    expect_true(all(sim1$D_c[, k] + sim1$D_t[, k] == d_per_look[k]))
    expect_true(all(sim2$D_c[, k] + sim2$D_t[, k] == d_per_look[k]))
  }
  ## ...but a treatment arm that recruits twice as fast (relative to
  ## control) picks up a larger share of those pooled events. accrual$rate
  ## is the pooled (both-arms) rate, so the *total* patient influx rate is
  ## the same (20) under both allocation ratios -- only the composition
  ## shifts. Here lambda_t is half of lambda_c (treatment survives longer,
  ## i.e. generates events more slowly per patient), so funnelling a larger
  ## share of the same total patient stream into the slower-event arm means
  ## the pooled event target is reached *later* in calendar time, not
  ## sooner -- the opposite of what a naive "faster recruitment -> shorter
  ## trial" intuition would suggest, and a useful check that the fix
  ## genuinely holds the total rate fixed rather than inflating it with
  ## allocation_ratio (which is what the pre-fix bug effectively did).
  final <- length(d_per_look)
  expect_gt(mean(sim2$D_t[, final]), mean(sim1$D_t[, final]))
  expect_gt(mean(sim2$C_look[, final]), mean(sim1$C_look[, final]))
})

test_that("simulated per-look event counts sum to the planned look sizes", {
  acc <- set_accrual(model = "poisson", rate = 20)
  psi <- list(psi_c = log(2) / 12, psi_t = log(2) / 24)
  set.seed(20260512L)
  sim <- adabay:::.sim_tte(surv_design(), psi,
                                   n_trials = 500L, accrual = acc)
  d_per_look <- surv_design()$schedule$D
  for (k in seq_along(d_per_look))
    expect_true(all(sim$D_c[, k] + sim$D_t[, k] == d_per_look[k]))
  ## Exposure is strictly positive and, with calendar cutoffs increasing in k,
  ## both events and exposure are non-decreasing across looks.
  expect_true(all(sim$E_c > 0))
  expect_true(all(sim$E_t > 0))
  for (k in 2:length(d_per_look)) {
    expect_true(all(sim$D_c[, k] >= sim$D_c[, k - 1]))
    expect_true(all(sim$D_t[, k] >= sim$D_t[, k - 1]))
    expect_true(all(sim$E_c[, k] >= sim$E_c[, k - 1]))
    expect_true(all(sim$E_t[, k] >= sim$E_t[, k - 1]))
  }
})

test_that("exposure obeys the operational-time identity E[lambda*E] = E[D]", {
  ## lambda_a E_a | D_a ~ Gamma(D_a, 1) by exponential memorylessness, so a
  ## faithful arrival-time + calendar-cutoff construction must satisfy
  ## E[lambda_a E_a] = E[D_a] at every look, for any recruitment rate. A bug
  ## that dropped the arrival offset (or summed exposure over not-yet-enrolled
  ## subjects) would inflate E and break this identity.
  psi <- list(psi_c = log(2) / 12, psi_t = log(2) / 24)
  d_per_look <- surv_design()$schedule$D
  for (rate in c(10, 40)) {
    acc <- set_accrual(model = "poisson", rate = rate)
    set.seed(rate)
    sim <- adabay:::.sim_tte(surv_design(), psi,
                                     n_trials = 20000L, accrual = acc)
    for (k in seq_along(d_per_look)) {
      expect_equal(mean(psi$psi_c * sim$E_c[, k]), mean(sim$D_c[, k]),
                   tolerance = 0.02)
      expect_equal(mean(psi$psi_t * sim$E_t[, k]), mean(sim$D_t[, k]),
                   tolerance = 0.02)
    }
  }
})
