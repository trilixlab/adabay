# Regression tests for the three blocking correctness bugs identified in review.
#   B1 - finite-mixture posterior must reweight component weights by the
#        per-trial marginal likelihood (supplement Eq. S11/S12), not use the
#        prior weights.
#   B2 - count + rate_difference must not over-read a scalar exposure vector in
#        the C++ tail routine (it returned NA/NaN operating characteristics).
#   B3 - standardised_mean_difference must standardise by sigma, i.e. it must
#        not be identical to mean_difference at a non-zero threshold.

test_that("B1: mixture posterior reweights component weights by marginal likelihood", {
  # Control prior: 0.5 * Beta(2,8) + 0.5 * Beta(20,5). Data S_c = 80/100 strongly
  # favour the high-rate component, so the *posterior* weights are far from 0.5.
  mk <- function(cw, ccomp)
    list(endpoint = "binary",
         arms = list(c = list(family = "beta", weights = cw, components = ccomp),
                     t = list(family = "beta", weights = 1,
                              components = list(list(a = 1, b = 1)))))
  comp1 <- list(a = 2, b = 8); comp2 <- list(a = 20, b = 5)
  n_c <- 100; S_c <- 80; n_t <- 100; S_t <- 50; e <- 0

  pkg <- adabay:::.tail_prob_binary(mk(c(0.5, 0.5), list(comp1, comp2)),
                           "risk_difference", n_c, n_t, S_c, S_t, e)[1, 1]
  P1  <- adabay:::.tail_prob_binary(mk(1, list(comp1)),
                           "risk_difference", n_c, n_t, S_c, S_t, e)[1, 1]
  P2  <- adabay:::.tail_prob_binary(mk(1, list(comp2)),
                           "risk_difference", n_c, n_t, S_c, S_t, e)[1, 1]

  bb <- function(a, b) exp(lbeta(a + S_c, b + n_c - S_c) - lbeta(a, b))
  ml <- c(bb(comp1$a, comp1$b), bb(comp2$a, comp2$b))
  wstar <- c(0.5, 0.5) * ml / sum(c(0.5, 0.5) * ml)   # correct posterior weights
  correct  <- wstar[1] * P1 + wstar[2] * P2
  priorwt  <- 0.5 * P1 + 0.5 * P2                      # the (buggy) prior-weighted value

  # Guard: the scenario must actually distinguish the two formulas.
  expect_false(isTRUE(all.equal(correct, priorwt)))
  expect_equal(pkg, correct, tolerance = 1e-8)
})

test_that("B1: single-component priors are unchanged (posterior weight == 1)", {
  pr <- list(endpoint = "binary",
             arms = list(c = list(family = "beta", weights = 1,
                                  components = list(list(a = 1, b = 1))),
                         t = list(family = "beta", weights = 1,
                                  components = list(list(a = 1, b = 1)))))
  v <- adabay:::.tail_prob_binary(pr, "risk_difference",
                         n_c = 50, n_t = 50, S_c = 20, S_t = 15, thresholds = 0)[1, 1]
  expect_true(is.finite(v) && v >= 0 && v <= 1)
})

test_that("B2: count + rate_difference returns finite operating characteristics", {
  des <- set_design(endpoint = "count",
                    exposure_per_look = c(2540, 5080),
                    effect_scale = "rate_difference")
  ## Shapes >= 1 keep the fixed-node quadrature accurate (see
  ## test-quadrature-shape.R); this test targets the B2 scalar-exposure
  ## recycling bug, not sub-unit-shape accuracy.
  pri <- set_prior(endpoint = "count",
                   arms = list(c = list(family = "gamma", a = 1, b = 1),
                               t = list(family = "gamma", a = 1, b = 1)))
  dec <- set_decision(efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
                      futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))
  oc <- evaluate_design(des, pri, dec,
                        effect = list(lambda_c = 0.02, lambda_t = 0.012),
                        n_trials = 500L, cores = 1L, seed = 1L)
  expect_true(is.finite(oc$alpha))
  expect_true(oc$alpha >= 0 && oc$alpha <= 1)
  expect_true(is.finite(oc$eff_prob))
})

test_that("B3: standardised_mean_difference standardises by sigma", {
  pri <- set_prior(endpoint = "continuous",
                   arms = list(c = list(family = "normal", mean = 0, sd = 1e3),
                               t = list(family = "normal", mean = 0, sd = 1e3)))
  run <- function(scale, thr, dnull) {
    des <- set_design(endpoint = "continuous", n_per_look = 400,
                      effect_scale = scale, alternative = "less",
                      sigma = 10, delta_null = dnull)
    dec <- set_decision(
      efficacy = list(list(threshold_effect = thr, threshold_prob = 0.975)),
      futility = list(list(threshold_effect = thr, threshold_prob = 0.999)))
    evaluate_design(des, pri, dec, effect = list(mu_c = 0, mu_t = -3.0),
                    n_trials = 3000L, cores = 1L, seed = 7L)
  }
  smd <- run("standardised_mean_difference", -0.3, -0.3)
  raw_same_thr <- run("mean_difference", -0.3, -0.3)
  # At a non-zero threshold the standardised scale must differ from the raw scale.
  expect_false(isTRUE(all.equal(smd$eff_prob, raw_same_thr$eff_prob)))
  # Correctness: standardised threshold e with sigma equals raw threshold e*sigma.
  raw_scaled <- run("mean_difference", -3.0, -3.0)
  expect_equal(smd$eff_prob, raw_scaled$eff_prob, tolerance = 1e-9)
})
