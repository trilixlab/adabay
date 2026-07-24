test_that("set_prior accepts conjugate beta priors", {
  p <- set_prior(endpoint = "binary",
                 arms = list(c = list(family = "beta", a = 2, b = 5),
                             t = list(family = "beta", a = 1, b = 1)))
  expect_s3_class(p, "adabay_prior")
  expect_equal(p$arms$c$weights, 1)
  expect_equal(p$arms$c$components[[1]]$a, 2)
})

test_that("set_prior accepts beta mixtures with weights summing to 1", {
  p <- set_prior(endpoint = "binary",
                 arms = list(
                   c = list(family = "beta",
                            components = list(list(a = 5, b = 10),
                                              list(a = 1, b = 1)),
                            weights = c(0.7, 0.3)),
                   t = list(family = "beta", a = 1, b = 1)))
  expect_equal(length(p$arms$c$weights), 2)
  expect_equal(sum(p$arms$c$weights), 1)
})

test_that("set_prior rejects mismatched endpoint and family", {
  expect_error(
    set_prior(endpoint = "binary",
              arms = list(c = list(family = "gamma", a = 1, b = 1),
                          t = list(family = "beta",  a = 1, b = 1))))
})

test_that("set_decision validates rule structure", {
  d <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
    futility = list(list(threshold_effect = 0, threshold_prob = 0.90)))
  expect_s3_class(d, "adabay_decision")
  expect_false(d$futility_binding)

  expect_error(set_decision(efficacy = list(list(threshold_effect = 0))))
  expect_error(set_decision(efficacy = list(list(threshold_effect = 0,
                                                 threshold_prob = 1.5))))
})

test_that("set_decision supports dual-criterion efficacy", {
  d <- set_decision(
    efficacy = list(list(threshold_effect =  0,    threshold_prob = 0.99),
                    list(threshold_effect = -0.05, threshold_prob = 0.50)),
    futility = list(list(threshold_effect =  0,    threshold_prob = 0.90)))
  expect_equal(length(d$efficacy$criteria), 2)
})
