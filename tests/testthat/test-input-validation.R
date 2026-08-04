# Input-domain validation hardening.
# Each bad input must fail fast with an informative error rather than
# producing NA / silently-wrong operating characteristics downstream.

test_that("M1: invalid data-generating effect values error informatively", {
  pri_b <- set_prior("binary", arms = list(c = list(family="beta",a=1,b=1),
                                            t = list(family="beta",a=1,b=1)))
  dec   <- set_decision(efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)))
  des_b <- set_design("binary", n_per_look = c(100,200),
                      effect_scale = "risk_difference")
  expect_error(evaluate_design(des_b, pri_b, dec,
                 effect = list(theta_c = 1.5, theta_t = 0.2), n_trials = 100L),
               "\\(0, 1\\)")
  expect_error(evaluate_design(des_b, pri_b, dec,
                 effect = list(theta_c = 0.3, theta_t = 0), n_trials = 100L),
               "\\(0, 1\\)")

  des_c <- set_design("count", exposure_per_look = c(200,400),
                      effect_scale = "log_rate_ratio")
  pri_c <- set_prior("count", arms = list(c = list(family="gamma",a=0.7,b=1),
                                          t = list(family="gamma",a=0.7,b=1)))
  expect_error(evaluate_design(des_c, pri_c, dec,
                 effect = list(lambda_c = -0.5, lambda_t = 0.3), n_trials = 100L),
               "positive")
  expect_error(evaluate_design(des_c, pri_c, dec,
                 effect = list(lambda_c = 0.1, lambda_t = 0), n_trials = 100L),
               "positive")
})

test_that("m1: non-integer / non-numeric n_per_look errors informatively", {
  expect_error(set_design("binary", n_per_look = c(20.5, 40)),
               "integer")
  expect_error(set_design("binary", n_per_look = c("twenty", "forty")),
               "numeric")
})

test_that("m2: tte d_total below max(d_per_look) is rejected", {
  expect_error(set_design("tte", d_per_look = c(50,100),
                          d_total = 80),
               "d_total")
})

test_that("m3: threshold_prob = NA errors informatively", {
  expect_error(set_decision(efficacy = list(list(threshold_effect = 0,
                                                  threshold_prob = NA))),
               "threshold_prob")
})

test_that("m5: evaluate_design on a non-adabay object errors informatively", {
  expect_error(evaluate_design(list(a = 1)), "adabay_design")
})

test_that("M2/M3: calibrate_design checks cache consistency and feasibility", {
  mk_cache <- function(np, tt) {
    des <- set_design("binary", n_per_look = np,
                      effect_scale = "risk_difference", alternative = "less")
    pri <- set_prior("binary", arms = list(c = list(family="beta",a=1,b=1),
                                           t = list(family="beta",a=1,b=1)))
    build_cache(des, pri, effect = list(theta_c = 0.33, theta_t = tt),
                n_trials = 300L, seed = 1L,
                threshold_grid = list(efficacy = 0, futility = 0))
  }
  c_h0  <- mk_cache(c(120,240,360), 0.33)
  c_h1  <- mk_cache(c(120,240,360), 0.28)
  c_bad <- mk_cache(c(80,160,240), 0.28)   # different look schedule

  expect_error(
    calibrate_design(c_h0, c_bad, alpha_target = 0.05, power_target = 0.5,
                     efficacy_grid = c(0.95, 0.99), futility_grid = c(0.9)),
    "same design")
  expect_warning(
    calibrate_design(c_h0, c_h1, alpha_target = 1e-4, power_target = 0.9999,
                     efficacy_grid = c(0.95, 0.99), futility_grid = c(0.9)),
    "target")
})
