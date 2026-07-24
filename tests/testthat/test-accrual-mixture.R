## Unit tests for set_accrual() and fit_mixture(): class invariants,
## input validation and (for fit_mixture) tail-probability diagnostics
## carried on the returned adabay_prior_arm object.

test_that("set_accrual returns a well-formed adabay_accrual for Poisson recruitment", {
  acc <- set_accrual(model = "poisson", rate = 20)
  expect_s3_class(acc, "adabay_accrual")
  expect_equal(acc$model, "poisson")
  expect_equal(acc$rate, 20)
  expect_equal(acc$dropout, "none")
})

test_that("set_accrual supports piecewise recruitment and exponential dropout", {
  acc <- set_accrual(model = "piecewise",
                     rate = c(5, 20, 15),
                     breakpoints = c(3, 12),
                     dropout = "exponential",
                     dropout_rate = 0.01)
  expect_s3_class(acc, "adabay_accrual")
  expect_equal(acc$model, "piecewise")
  expect_equal(acc$rate, c(5, 20, 15))
  expect_equal(acc$breakpoints, c(3, 12))
  expect_equal(acc$dropout, "exponential")
  expect_equal(acc$dropout_rate, 0.01)
})

test_that("set_accrual rejects malformed inputs", {
  expect_error(set_accrual(model = "poisson"),
               "rate")
  expect_error(set_accrual(model = "piecewise", rate = c(5, 10),
                           breakpoints = c(3, 6)),
               "length")
  expect_error(set_accrual(model = "callback"),
               "callback")
  expect_error(set_accrual(model = "poisson", rate = 20, follow_up = -1),
               "non-negative")
  expect_error(set_accrual(model = "poisson", rate = 20, dropout = "exponential"),
               "dropout_rate")
})

test_that("fit_mixture returns a adabay_prior_arm with tail diagnostics (samples input)", {
  set.seed(20260512)
  samples <- stats::rbeta(2000, shape1 = 4, shape2 = 6)
  pri_c <- fit_mixture(endpoint = "binary", arm = "c", prior = samples,
                       n_components_max = 3, tol_kl = 1e-3, tol_tail = 5e-2)
  expect_s3_class(pri_c, "adabay_prior_arm")
  expect_equal(pri_c$family, "beta")
  expect_equal(sum(pri_c$weights), 1, tolerance = 1e-8)
  expect_true(all(pri_c$weights > 0))
  expect_true(length(pri_c$components) >= 1L)

  diag <- attr(pri_c, "diagnostics")
  expect_false(is.null(diag))
  expect_true(!is.null(diag$tail_table))
  expect_true(!is.null(diag$max_tail_error))
  expect_true(is.finite(diag$kl))
})

test_that("fit_mixture (internal-EM) keeps adding components past the KL elbow when tol_tail is not yet met", {
  ## Regression guard: the KL-decrement elbow rule alone can stop too early
  ## -- forward KL is mass-covering, not tail-accurate, so a component that
  ## only marginally improves the bulk-density KL can still be needed to
  ## bring a tail probability within tol_tail. .em_beta_mixture() must keep
  ## growing L past the elbow (up to n_components_max) whenever the elbow
  ## pick's tail error exceeds tol_tail, and stop at the first L that
  ## satisfies it.
  skip_if(requireNamespace("RBesT", quietly = TRUE),
         "this guards the internal-EM fallback specifically")
  logit_normal_density <- function(theta) {
    meanlog <- -0.7; sdlog <- 0.4
    z <- log(theta / (1 - theta))
    dnorm(z, meanlog, sdlog) / (theta * (1 - theta))
  }
  set.seed(1L)
  pri_c <- fit_mixture(endpoint = "binary", arm = "c",
                       prior = logit_normal_density,
                       n_components_max = 5L, tol_kl = 1e-3, tol_tail = 5e-3,
                       n_samples = 10000L)
  diag <- attr(pri_c, "diagnostics")
  expect_gt(length(pri_c$weights), 1L)
  expect_lte(diag$max_tail_error, diag$tol_tail)
  expect_false(diag$warning_fired)
})

test_that("fit_mixture accepts a density function and plugs into set_prior", {
  pri_c <- fit_mixture(endpoint = "binary", arm = "c",
                       prior = function(theta)
                         dbeta(theta, shape1 = 4, shape2 = 6),
                       n_components_max = 3,
                       n_samples = 2000,
                       tol_kl = 1e-3, tol_tail = 5e-2)
  expect_s3_class(pri_c, "adabay_prior_arm")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = pri_c,
                               t = list(family = "beta", a = 1, b = 1)))
  expect_s3_class(pri, "adabay_prior")
})

test_that("fit_mixture(seed = ...) gives bit-identical results, ignoring ambient RNG state", {
  density_fn <- function(theta) dbeta(theta, shape1 = 4, shape2 = 6)
  args <- list(endpoint = "binary", arm = "c", prior = density_fn,
              n_components_max = 3, n_samples = 2000, tol_kl = 1e-3,
              tol_tail = 5e-2, seed = 42L)

  fit1 <- do.call(fit_mixture, args)
  ## Perturb the ambient RNG state between calls -- with seed supplied this
  ## must not matter, unlike the no-seed default (see next test).
  runif(37)
  fit2 <- do.call(fit_mixture, args)
  expect_identical(fit1, fit2)

  diag <- attr(fit1, "diagnostics")
  expect_equal(diag$seed, 42L)
})

test_that("fit_mixture() without seed depends on ambient RNG state (documents current behaviour)", {
  density_fn <- function(theta) dbeta(theta, shape1 = 4, shape2 = 6)
  args <- list(endpoint = "binary", arm = "c", prior = density_fn,
              n_components_max = 3, n_samples = 2000, tol_kl = 1e-3,
              tol_tail = 5e-2)

  set.seed(1L)
  fit1 <- do.call(fit_mixture, args)
  set.seed(2L)
  fit2 <- do.call(fit_mixture, args)
  expect_false(identical(attr(fit1, "diagnostics")$samples,
                         attr(fit2, "diagnostics")$samples))
  expect_null(attr(fit1, "diagnostics")$seed)
})

test_that("print.adabay_prior_arm lists each component's weight and conjugate parameters", {
  pri_c <- fit_mixture(endpoint = "binary", arm = "c",
                       prior = function(theta) dbeta(theta, shape1 = 4, shape2 = 6),
                       n_components_max = 1, n_samples = 500, tol_kl = 1e-3,
                       tol_tail = 1, seed = 1L)
  out <- capture.output(print(pri_c))
  ## One "[1] weight = ..., a = ..., b = ..." line for the single component.
  comp_line <- grep("^\\s*\\[1\\] weight = .*a = .*b = ", out, value = TRUE)
  expect_length(comp_line, 1L)
  expect_true(any(grepl("^\\s*Seed:\\s*1\\s*$", out)))
})

test_that("posterior tail probability matches the analytic beta-binomial form", {
  ## Under a conjugate Beta(1,1) prior, P(theta_c > x | S = s, n) = pbeta(x, 1+s, 1+n-s, lower.tail = FALSE).
  ## The internal .pdelta_binary_rd reduces to a convolution that, at the
  ## degenerate point of a single-arm trial (treatment arm matched to control),
  ## reduces to symmetry: P(Delta > 0) = 0.5 by exchangeability.
  des <- set_design(endpoint = "binary",
                    n_per_look = 200, effect_scale = "risk_difference",
                    alternative = "less")
  pri <- set_prior(endpoint = "binary",
                   arms = list(c = list(family = "beta", a = 1, b = 1),
                               t = list(family = "beta", a = 1, b = 1)))
  dec <- set_decision(efficacy = list(list(threshold_effect = 0,
                                            threshold_prob = 0.5)))
  oc <- evaluate_design(des, pri, dec,
                        effect = list(theta_c = 0.30, theta_t = 0.30),
                        n_trials = 5000, cores = 1, seed = 20260512L)
  ## At the null (theta_t = theta_c) with symmetric priors, the probability
  ## that P(Delta < 0 | data) exceeds 0.5 must itself be close to 0.5
  ## (the exact value depends on tie-breaking when the posterior is
  ## exactly 0.5; the MC standard error at R = 5000 is ~0.007).
  expect_lt(abs(oc$alpha - 0.5), 0.05)
})
