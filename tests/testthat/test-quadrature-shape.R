# Guard for the fixed-node Gauss-Legendre quadrature used on difference and
# binary relative-effect scales: it loses accuracy when a posterior shape
# parameter falls below 1 (integrable density-endpoint singularity). The entry
# points warn in that case; the closed-form relative-effect scales never do.

test_that("evaluate_design warns for a sub-unit prior shape on a quadrature scale", {
  des <- set_design(endpoint = "binary", n_per_look = c(100, 200),
                    effect_scale = "risk_difference", alternative = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 0.5, b = 0.5),
                               t = list(family = "beta", a = 0.5, b = 0.5)))
  dec <- set_decision(efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
                      futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))
  expect_warning(
    evaluate_design(des, pri, dec, effect = list(theta_c = 0.3, theta_t = 0.2),
                    n_trials = 200L, cores = 1L, seed = 1L),
    "fixed-node Gauss-Legendre quadrature"
  )
})

test_that("no quadrature warning when all prior shapes are at least 1", {
  des <- set_design(endpoint = "binary", n_per_look = c(100, 200),
                    effect_scale = "risk_difference", alternative = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
                      futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))
  expect_no_warning(
    evaluate_design(des, pri, dec, effect = list(theta_c = 0.3, theta_t = 0.2),
                    n_trials = 200L, cores = 1L, seed = 1L)
  )
})

test_that("closed-form relative-effect scale never warns, even at shape < 1", {
  des <- set_design(endpoint = "count", exposure_per_look = c(2540, 5080),
                    effect_scale = "log_rate_ratio", alternative = "less",
                    delta_null = log(0.70))
  pri <- set_prior(endpoint = "count",
                   arms = list(c = list(family = "gamma", a = 0.5, b = 1),
                               t = list(family = "gamma", a = 0.5, b = 1)))
  dec <- set_decision(
    efficacy = list(list(threshold_effect = log(0.70), threshold_prob = 0.99)),
    futility = list(list(threshold_effect = log(0.70), threshold_prob = 0.90)))
  expect_no_warning(
    evaluate_design(des, pri, dec, effect = list(lambda_c = 0.02, lambda_t = 0.012),
                    n_trials = 200L, cores = 1L, seed = 1L)
  )
})

test_that("build_cache warns for a sub-unit prior shape on a quadrature scale", {
  des <- set_design(endpoint = "count", exposure_per_look = c(2540, 5080),
                    effect_scale = "rate_difference")
  pri <- set_prior(endpoint = "count",
                   arms = list(c = list(family = "gamma", a = 0.5, b = 1),
                               t = list(family = "gamma", a = 0.5, b = 1)))
  expect_warning(
    build_cache(des, pri, effect = list(lambda_c = 0.02, lambda_t = 0.012),
                n_trials = 200L, cores = 1L, seed = 1L,
                threshold_grid = list(efficacy = 0, futility = 0)),
    "fixed-node Gauss-Legendre quadrature"
  )
})
