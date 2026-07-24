test_that("set_design rejects bad inputs", {
  expect_error(set_design(endpoint = "binary"))
  expect_error(set_design(endpoint = "binary", n_per_look = integer(0)))
  expect_error(set_design(endpoint = "binary",
                          n_per_look = c(10, 5, 7)))
})

test_that("set_design returns a well-formed object for binary", {
  d <- set_design(endpoint = "binary",
                  n_per_look = c(200, 400, 600),
                  effect_scale = "risk_difference",
                  alternative = "less")
  expect_s3_class(d, "adabay_design")
  expect_equal(d$schedule$n_c, c(100L, 200L, 300L))
  expect_equal(d$schedule$n_t, c(100L, 200L, 300L))
  expect_equal(d$delta_null, 0)
  expect_equal(d$alternative, "less")
  expect_equal(d$n_looks, 3L)
})

test_that("set_design infers n_looks from the length of the per-look vector", {
  ## No separate n_looks argument exists: the number of looks is always
  ## length(n_per_look) / length(exposure_per_look) / length(d_per_look).
  expect_equal(
    set_design(endpoint = "binary", n_per_look = 200,
              effect_scale = "risk_difference")$n_looks,
    1L)
  expect_equal(
    set_design(endpoint = "count", exposure_per_look = c(50, 100, 150, 200))$n_looks,
    4L)
})

test_that("set_design supports continuous, count and tte endpoints", {
  expect_s3_class(
    set_design(endpoint = "continuous",
               n_per_look = c(100, 200), sigma = 1),
    "adabay_design")
  expect_s3_class(
    set_design(endpoint = "count",
               exposure_per_look = c(200, 400)),
    "adabay_design")
  expect_s3_class(
    set_design(endpoint = "tte",
               d_total = 100, d_per_look = c(33, 66, 100)),
    "adabay_design")
})

test_that("set_design errors when allocation_ratio makes the per-arm split non-increasing", {
  expect_error(
    set_design(endpoint = "binary",
              n_per_look = c(100, 101),
              effect_scale = "risk_difference", alternative = "less"),
    "increments too small")
})
