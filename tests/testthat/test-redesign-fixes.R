# Regression tests for two correctness fixes in the pooled-API redesign.

test_that("set_design rejects a pooled schedule that rounds an arm to zero", {
  # allocation_ratio = 9 rounds the control arm to 0 at the first look:
  # n_c = round(c(5,30,60)/10) = c(0, 3, 6). The per-arm split is still
  # strictly increasing, so only the positivity guard catches it.
  expect_error(
    set_design("continuous", n_per_look = c(5, 30, 60),
               allocation_ratio = 9, sigma = 1),
    "positive size"
  )
  # Reachable at the default allocation_ratio = 1 too:
  # n_c = round(c(1,4,6)/2) = c(0, 2, 3).
  expect_error(
    set_design("binary", n_per_look = c(1, 4, 6), alternative = "less"),
    "positive size"
  )
  # A schedule that gives both arms a positive, strictly increasing size is fine.
  des <- set_design("binary", n_per_look = c(760, 1520, 2280, 3040, 3800),
                    effect_scale = "risk_difference", alternative = "less")
  expect_s3_class(des, "adabay_design")
})

test_that(".em_beta_mixture returns the least-tail-error fit, not blindly the largest L", {
  # With an unsatisfiable tol_tail the fitter must fall back to the fit with
  # the smallest maximum tail error among those searched, NOT fits[[L_max]]
  # (more components is a fresh EM fit and does not monotonically reduce the
  # tail error). For this reproducible sample the L = 2 fit beats L = 5.
  set.seed(1)
  x  <- stats::rbeta(3000, 0.6, 4)
  tq <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
  ncm <- 5L
  max_tail <- function(fit)
    max(abs(adabay:::.tail_error_table(x, "beta", fit$weights,
                                       fit$components, tq)$error))

  fit      <- adabay:::.em_beta_mixture(x, n_components_max = ncm, tol_kl = 1e-3,
                                        tail_quantiles = tq, tol_tail = 0)
  ret_err  <- max_tail(fit)
  last_err <- max_tail(adabay:::.em_beta(x, L = ncm))

  # The buggy version returned fits[[ncm]] (== last_err); the fix returns a
  # strictly better fit here, and in general never a worse one.
  expect_lt(ret_err, last_err)
})
