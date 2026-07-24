test_that("binary semi-simulation runs end-to-end and gives sane output", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))

  oc_h0 <- evaluate_design(des, pri, dec,
                           effect = list(theta_c = 0.30, theta_t = 0.30),
                           n_trials = 200, cores = 1, seed = 42L)
  oc_h1 <- evaluate_design(des, pri, dec,
                           effect = list(theta_c = 0.30, theta_t = 0.10),
                           n_trials = 200, cores = 1, seed = 42L)

  expect_s3_class(oc_h0, "adabay_oc")
  expect_true(is.finite(oc_h0$alpha))
  expect_true(is.finite(oc_h1$alpha))
  expect_true(oc_h1$alpha > oc_h0$alpha)
  expect_true(oc_h0$alpha < 0.5)
  expect_true(oc_h0$expected_sample_size <= 300)
})

test_that("continuous semi-simulation runs and respects direction", {
  des <- set_design(endpoint = "continuous",
                    n_per_look = c(40, 80, 160),
                    effect_scale = "mean_difference",
                    alternative  = "greater",
                    sigma = 1)
  pri <- set_prior(endpoint = "continuous",
                   arms = list(c = list(family = "normal", mean = 0, sd = 5),
                               t = list(family = "normal", mean = 0, sd = 5)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.95)))
  oc_h1 <- evaluate_design(des, pri, dec,
                           effect = list(mu_c = 0, mu_t = 0.6),
                           n_trials = 200, cores = 1, seed = 1L)
  expect_s3_class(oc_h1, "adabay_oc")
  expect_true(oc_h1$alpha > 0)
})

test_that("count and tte paths run end-to-end", {
  des_c <- set_design(endpoint = "count",
                      exposure_per_look = c(100, 200),
                      effect_scale = "log_rate_ratio",
                      alternative = "less")
  pri_c <- set_prior(endpoint = "count",
                     arms = list(c = list(family = "gamma", a = 1, b = 1),
                                 t = list(family = "gamma", a = 1, b = 1)))
  dec_c <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)))
  oc_c <- evaluate_design(des_c, pri_c, dec_c,
                          effect = list(lambda_c = 1.0, lambda_t = 0.5),
                          n_trials = 100, cores = 1, seed = 1L)
  expect_s3_class(oc_c, "adabay_oc")

  des_s <- set_design(endpoint = "tte",
                      d_total = 50, d_per_look = c(25, 50),
                      effect_scale = "log_hazard_ratio",
                      alternative = "less")
  pri_s <- set_prior(endpoint = "tte",
                     arms = list(c = list(family = "gamma", a = 1, b = 1),
                                 t = list(family = "gamma", a = 1, b = 1)))
  dec_s <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)))
  acc_s <- set_accrual(model = "poisson", rate = 20)
  oc_s <- evaluate_design(des_s, pri_s, dec_s,
                          effect = list(lambda_c = log(2)/12,
                                        lambda_t = log(2)/24),
                          accrual = acc_s,
                          n_trials = 50, cores = 1, seed = 1L)
  expect_s3_class(oc_s, "adabay_oc")
})

test_that("cached path (build_cache + evaluate_design) matches fused path (evaluate_design on a adabay_design)", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))

  oc1 <- evaluate_design(des, pri, dec,
                         effect = list(theta_c = 0.30, theta_t = 0.20),
                         n_trials = 200, cores = 1, seed = 7L)

  cache <- build_cache(des, pri,
                       effect = list(theta_c = 0.30, theta_t = 0.20),
                       n_trials = 200, cores = 1, seed = 7L,
                       threshold_grid = list(efficacy = 0, futility = 0))
  oc2 <- evaluate_design(cache, dec)

  expect_equal(oc1$alpha, oc2$alpha)
  expect_equal(oc1$fut_prob,  oc2$fut_prob)
})

test_that("cross_core_reproducible = TRUE gives bit-identical OCs at cores 1, 2, 4 (fused path)", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
    futility_binding = TRUE)
  args <- list(prior = pri, decision = dec,
               effect = list(theta_c = 0.30, theta_t = 0.20),
               n_trials = 200L, seed = 13L,
               cross_core_reproducible = TRUE)
  oc_1 <- do.call(evaluate_design, c(list(x = des, cores = 1L), args))
  oc_2 <- do.call(evaluate_design, c(list(x = des, cores = 2L), args))
  oc_4 <- do.call(evaluate_design, c(list(x = des, cores = 4L), args))
  expect_identical(oc_1$alpha, oc_2$alpha)
  expect_identical(oc_1$alpha, oc_4$alpha)
  expect_identical(oc_1$eff_prob, oc_2$eff_prob)
  expect_identical(oc_1$eff_prob, oc_4$eff_prob)
  expect_identical(oc_1$eff_prob_at_look, oc_2$eff_prob_at_look)
  expect_identical(oc_1$eff_prob_at_look, oc_4$eff_prob_at_look)
  expect_identical(oc_1$fut_prob_at_look, oc_2$fut_prob_at_look)
  expect_identical(oc_1$fut_prob_at_look, oc_4$fut_prob_at_look)
  expect_identical(oc_1$expected_sample_size, oc_2$expected_sample_size)
  expect_identical(oc_1$expected_sample_size, oc_4$expected_sample_size)
})

test_that("cross_core_reproducible = TRUE gives bit-identical caches at cores 1, 2, 4 (cached path)", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
    futility_binding = TRUE)
  build_args <- list(prior = pri,
                     effect = list(theta_c = 0.30, theta_t = 0.20),
                     n_trials = 200L, seed = 17L,
                     threshold_grid = list(efficacy = 0, futility = 0),
                     cross_core_reproducible = TRUE)
  cache_1 <- do.call(build_cache, c(list(design = des, cores = 1L), build_args))
  cache_2 <- do.call(build_cache, c(list(design = des, cores = 2L), build_args))
  cache_4 <- do.call(build_cache, c(list(design = des, cores = 4L), build_args))
  expect_identical(cache_1$tails, cache_2$tails)
  expect_identical(cache_1$tails, cache_4$tails)
  oc_1 <- evaluate_design(cache_1, dec)
  oc_2 <- evaluate_design(cache_2, dec)
  oc_4 <- evaluate_design(cache_4, dec)
  expect_identical(oc_1$alpha, oc_2$alpha)
  expect_identical(oc_1$alpha, oc_4$alpha)
})
