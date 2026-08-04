## Regression tests: pin the headline operating characteristics returned
## by evaluate_design() at a fixed (small) seed and budget, separately
## for each of the four endpoint types. A mismatch flags a change in the
## simulation pipeline, the per-look posterior path, or the aggregator,
## even when the unit tests still pass.
##
## The pinned numbers were generated under R 4.x with the L'Ecuyer-CMRG
## streams used internally by .par_apply; bit-identical reproducibility
## is therefore expected on the same R/platform. Tolerances allow for
## small numerical drift across BLAS or platform changes.

test_that("binary semi-simulation matches stored operating characteristics", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))

  oc <- evaluate_design(des, pri, dec,
                        effect = list(theta_c = 0.30, theta_t = 0.15),
                        n_trials = 500, cores = 1, seed = 20260512L)

  expect_equal(oc$alpha, 0.958, tolerance = 1e-6)
  expect_equal(oc$fut_prob,  0.042, tolerance = 1e-6)
  expect_equal(oc$expected_sample_size, 158.2, tolerance = 1e-6)
  expect_equal(oc$eff_prob_at_look, c(0.542, 0.330, 0.086), tolerance = 1e-6)
  expect_equal(oc$fut_prob_at_look,  c(0.002, 0.000, 0.040), tolerance = 1e-6)
})

test_that("continuous semi-simulation matches stored operating characteristics", {
  des <- set_design(endpoint = "continuous",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "mean_difference",
                    alternative  = "greater",
                    sigma = 1)
  pri <- set_prior(endpoint = "continuous",
                   arms = list(c = list(family = "normal", mean = 0, sd = 10),
                               t = list(family = "normal", mean = 0, sd = 10)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))

  oc <- evaluate_design(des, pri, dec,
                        effect = list(mu_c = 0, mu_t = 0.3),
                        n_trials = 500, cores = 1, seed = 20260512L)

  expect_equal(oc$alpha, 0.88, tolerance = 1e-6)
  expect_equal(oc$fut_prob,  0.12, tolerance = 1e-6)
  expect_equal(oc$expected_sample_size, 184.4, tolerance = 1e-6)
})

test_that("normal-inverse-gamma continuous prior works end-to-end", {
  ## Guards against accidental clobbering of the `beta` hyperparameter of
  ## the NIG prior. The partition `eff_prob + fut_prob == 1` and finite
  ## non-NA outputs are the minimum invariants to enforce.
  des <- set_design(endpoint = "continuous",
                    n_per_look = c(60, 120, 180),
                    effect_scale = "mean_difference",
                    alternative  = "greater",
                    sigma = 1)
  pri <- set_prior(endpoint = "continuous",
                   arms = list(c = list(family = "normal_inverse_gamma",
                                        nu = 0, kappa = 1, alpha = 2, beta = 1),
                               t = list(family = "normal_inverse_gamma",
                                        nu = 0, kappa = 1, alpha = 2, beta = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
    futility_binding = TRUE)
  oc <- evaluate_design(des, pri, dec,
                        effect = list(mu_c = 0, mu_t = 0.5),
                        n_trials = 500L, cores = 1L, seed = 1L)
  expect_equal(sum(oc$eff_prob_at_look) + sum(oc$fut_prob_at_look), 1)
  expect_false(any(is.na(c(oc$alpha, oc$eff_prob, oc$fut_prob,
                            oc$eff_prob_at_look, oc$fut_prob_at_look))))
  expect_true(oc$alpha > 0 && oc$alpha < 1)
  expect_true(oc$expected_sample_size > 0)
})

test_that("count semi-simulation matches stored operating characteristics", {
  des <- set_design(endpoint = "count",
                    exposure_per_look = c(40, 80, 120),
                    effect_scale = "log_rate_ratio",
                    alternative  = "less")
  pri <- set_prior(endpoint = "count",
                   arms = list(c = list(family = "gamma", a = 1, b = 1),
                               t = list(family = "gamma", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)))

  oc <- evaluate_design(des, pri, dec,
                        effect = list(lambda_c = 1.0, lambda_t = 0.5),
                        n_trials = 500, cores = 1, seed = 20260512L)

  expect_equal(oc$alpha, 0.952, tolerance = 1e-6)
  expect_equal(oc$fut_prob,  0.048, tolerance = 1e-6)
  expect_equal(oc$expected_sample_size, 60.48, tolerance = 1e-6)
})

test_that("tte semi-simulation matches stored operating characteristics", {
  des <- set_design(endpoint = "tte",
                    d_total = 60, d_per_look = c(20, 40, 60),
                    effect_scale = "log_hazard_ratio",
                    alternative  = "less")
  pri <- set_prior(endpoint = "tte",
                   arms = list(c = list(family = "gamma", a = 1, b = 1),
                               t = list(family = "gamma", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)))
  acc <- set_accrual(model = "poisson", rate = 40)

  oc <- evaluate_design(des, pri, dec,
                        effect = list(lambda_c = log(2) / 12,
                                      lambda_t = log(2) / 24),
                        accrual = acc,
                        n_trials = 200, cores = 1, seed = 20260512L)

  ## accrual$rate is the pooled (both-arms) rate, so rate = 40 here is
  ## 20 per arm at the 1:1 allocation. The pinned values below are tied to
  ## that convention, because .sim_tte()'s enrolment-pool sizing depends on
  ## the per-arm rates derived from the pooled rate.
  expect_equal(oc$alpha, 0.850, tolerance = 1e-6)
  expect_equal(oc$fut_prob,  0.150, tolerance = 1e-6)
})

test_that("repeated calls at the same seed and cores give bit-identical results", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))

  oc1 <- evaluate_design(des, pri, dec,
                         effect = list(theta_c = 0.30, theta_t = 0.20),
                         n_trials = 200, cores = 1, seed = 20260512L)
  oc2 <- evaluate_design(des, pri, dec,
                         effect = list(theta_c = 0.30, theta_t = 0.20),
                         n_trials = 200, cores = 1, seed = 20260512L)

  expect_equal(oc1$alpha, oc2$alpha)
  expect_equal(oc1$fut_prob,  oc2$fut_prob)
  expect_equal(oc1$eff_prob_at_look, oc2$eff_prob_at_look)
  expect_equal(oc1$fut_prob_at_look,  oc2$fut_prob_at_look)
  expect_equal(oc1$expected_sample_size, oc2$expected_sample_size)
})

test_that("calibrate_design requires cache_alt and matches direct sweep", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200, 300),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  grid <- list(efficacy = 0, futility = 0)
  cache_h0 <- build_cache(des, pri,
                          effect = list(theta_c = 0.30, theta_t = 0.30),
                          n_trials = 500, cores = 1, seed = 20260512L,
                          threshold_grid = grid)
  cache_h1 <- build_cache(des, pri,
                          effect = list(theta_c = 0.30, theta_t = 0.15),
                          n_trials = 500, cores = 1, seed = 20260512L,
                          threshold_grid = grid)
  expect_error(
    calibrate_design(cache_h0,
                     alpha_target = 0.05, power_target = 0.5,
                     efficacy_grid = c(0.90, 0.95),
                     futility_grid = c(1)),
    "cache_alt"
  )
  cal <- calibrate_design(cache_h0, cache_h1,
                          alpha_target = 0.5, power_target = 0.1,
                          efficacy_grid = c(0.50, 0.90, 0.95),
                          futility_grid = c(1))
  expect_s3_class(cal, "adabay_calibration")
  expect_s3_class(cal$best, "adabay_calibration_best")
  expect_true(all(c("p", "q", "type_I", "power", "E_N") %in% names(cal$grid)))
  expect_true(all(cal$grid$type_I >= 0 & cal$grid$type_I <= 1))
  expect_true(all(cal$grid$power  >= 0 & cal$grid$power  <= 1))

  ## Guard: cal$best$type_I / cal$best$power must come from the H0 / H1
  ## caches respectively. This is why $best carries its own class rather
  ## than being an adabay_oc: an adabay_oc's $alpha field evaluated under
  ## H1 is the power, which the generic print method would label
  ## "Type I error".
  best_row <- cal$grid[cal$grid$p == cal$best$p & cal$grid$q == cal$best$q, ]
  expect_equal(cal$best$type_I, best_row$type_I)
  expect_equal(cal$best$power, best_row$power)
  expect_equal(cal$best$expected_sample_size, best_row$E_N)
  expect_s3_class(cal$best$oc_h0, "adabay_oc")
  expect_s3_class(cal$best$oc_h1, "adabay_oc")
  expect_equal(cal$best$type_I, cal$best$oc_h0$alpha)
  expect_equal(cal$best$power, cal$best$oc_h1$alpha)
})

test_that("summary fields are present and self-consistent", {
  des <- set_design(endpoint = "binary",
                    n_per_look = c(100, 200),
                    effect_scale = "risk_difference",
                    alternative  = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.95)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))
  oc <- evaluate_design(des, pri, dec,
                        effect = list(theta_c = 0.30, theta_t = 0.15),
                        n_trials = 200, cores = 1, seed = 20260512L)
  s <- summary(oc)
  expect_named(s, c("alpha", "eff_prob", "fut_prob", "ess", "edur",
                    "n_trials", "eff_prob_at_look", "fut_prob_at_look"))
  expect_equal(s$alpha, oc$alpha)
  expect_equal(s$eff_prob, oc$eff_prob)
  expect_equal(s$fut_prob, oc$fut_prob)
  expect_equal(s$ess, oc$expected_sample_size)
  expect_equal(s$n_trials, oc$n_trials)

  tab <- summarise_oc(oc)
  expect_true(all(c("alpha", "eff_prob", "fut_prob", "E_N",
                    "E_T", "n_trials", "n_looks") %in% names(tab)))
  expect_equal(tab$alpha,    oc$alpha)
  expect_equal(tab$eff_prob, oc$eff_prob)
  expect_equal(tab$fut_prob, oc$fut_prob)
})
