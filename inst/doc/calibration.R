## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----design, message = FALSE--------------------------------------------------
library(adabay)

des <- set_design(
  endpoint     = "binary",
  n_per_look   = c(760, 1520, 2280, 3040, 3800),
  effect_scale = "risk_difference",
  alternative  = "less"
)

R_CACHE <- 5000L

threshold_grid <- list(efficacy = 0, futility = 0)

## ----conj-prior, message = FALSE----------------------------------------------
pri_conj <- set_prior(
  endpoint = "binary",
  arms = list(c = list(family = "beta", a = 1, b = 1),
              t = list(family = "beta", a = 1, b = 1))
)

cache_conj_h0 <- build_cache(
  des, pri_conj,
  effect   = list(theta_c = 0.33, theta_t = 0.33),
  n_trials = R_CACHE, cores = 1, seed = 1L,
  threshold_grid = threshold_grid
)
cache_conj_h1 <- build_cache(
  des, pri_conj,
  effect   = list(theta_c = 0.33, theta_t = 0.28),
  n_trials = R_CACHE, cores = 1, seed = 1L,
  threshold_grid = threshold_grid
)

## ----conj-calibrate-v1--------------------------------------------------------
cal_conj_90 <- calibrate_design(
  cache         = cache_conj_h0,
  cache_alt     = cache_conj_h1,
  alpha_target  = 0.025,
  power_target  = 0.90,
  efficacy_grid = seq(0.95, 0.999, by = 0.005),
  futility_grid = seq(0.80, 0.95,  by = 0.05)
)
is.null(cal_conj_90$best)

## ----conj-calibrate-----------------------------------------------------------
cal_conj <- calibrate_design(
  cache         = cache_conj_h0,
  cache_alt     = cache_conj_h1,
  alpha_target  = 0.025,
  power_target  = 0.80,
  efficacy_grid = seq(0.95, 0.999, by = 0.005),
  futility_grid = seq(0.80, 0.95,  by = 0.05)
)
head(cal_conj$grid)

## ----conj-best----------------------------------------------------------------
if (!is.null(cal_conj$best)) {
  cal_conj$best
} else {
  cat("No (p, q) in the searched grid met both alpha and power targets",
      "at R =", R_CACHE, "per cache; raise R or widen the grid.\n")
}

## ----conj-stage---------------------------------------------------------------
p_stage <- c(0.998, 0.998, 0.998, 0.998, 0.9775)
q_stage <- 0.90

dec_stage <- set_decision(
  efficacy = lapply(p_stage, function(p)
    list(list(threshold_effect = 0, threshold_prob = p))),
  futility = list(list(threshold_effect = 0, threshold_prob = q_stage)),
  futility_binding = TRUE
)

oc_conj_stage_h0 <- evaluate_design(cache_conj_h0, dec_stage)
oc_conj_stage_h1 <- evaluate_design(cache_conj_h1, dec_stage)

cat(sprintf("Conjugate prior, stage-specific p_1:5 = (0.998, 0.998, 0.998, 0.998, 0.9775):\n"))
cat(sprintf("  Type I error rate  = %.2f%%\n", 100 * oc_conj_stage_h0$alpha))
cat(sprintf("  Power              = %.2f%%\n", 100 * oc_conj_stage_h1$alpha))
cat(sprintf("  E(N | H1)          = %.0f patients\n",
            oc_conj_stage_h1$expected_sample_size))

## ----conj-plot, fig.alt = "Per-look stopping probabilities for the stage-specific calibrated design under the conjugate prior."----
if (requireNamespace("ggplot2", quietly = TRUE)) plot(oc_conj_stage_h1)

## ----nonconj-prior, message = FALSE, warning = FALSE--------------------------
logit_normal_density <- function(theta) {
  meanlog <- -0.7
  sdlog   <-  0.4
  z <- log(theta / (1 - theta))
  dnorm(z, meanlog, sdlog) / (theta * (1 - theta))
}

pri_c_arm <- fit_mixture(
  endpoint         = "binary",
  arm              = "c",
  prior            = logit_normal_density,
  n_components_max = 5L,
  tol_kl           = 1e-3,
  tol_tail         = 5e-3,
  n_samples        = 10000L,
  seed             = 1L
)
pri_c_arm

pri_nc <- set_prior(
  endpoint = "binary",
  arms     = list(c = pri_c_arm,
                  t = list(family = "beta", a = 1, b = 1))
)

## ----nonconj-caches, message = FALSE------------------------------------------
cache_nc_h0 <- build_cache(
  des, pri_nc,
  effect   = list(theta_c = 0.33, theta_t = 0.33),
  n_trials = R_CACHE, cores = 1, seed = 1L,
  threshold_grid = threshold_grid
)
cache_nc_h1 <- build_cache(
  des, pri_nc,
  effect   = list(theta_c = 0.33, theta_t = 0.28),
  n_trials = R_CACHE, cores = 1, seed = 1L,
  threshold_grid = threshold_grid
)

## ----nonconj-calibrate--------------------------------------------------------
cal_nc <- calibrate_design(
  cache         = cache_nc_h0,
  cache_alt     = cache_nc_h1,
  alpha_target  = 0.025,
  power_target  = 0.80,
  efficacy_grid = seq(0.95, 0.999, by = 0.005),
  futility_grid = seq(0.80, 0.95,  by = 0.05)
)
head(cal_nc$grid)
cal_nc$best

## ----nonconj-stage------------------------------------------------------------
oc_nc_stage_h0 <- evaluate_design(cache_nc_h0, dec_stage)
oc_nc_stage_h1 <- evaluate_design(cache_nc_h1, dec_stage)

cat(sprintf("Non-conjugate prior, stage-specific p_1:5 = (0.998, 0.998, 0.998, 0.998, 0.9775):\n"))
cat(sprintf("  Type I error rate  = %.2f%%\n", 100 * oc_nc_stage_h0$alpha))
cat(sprintf("  Power              = %.2f%%\n", 100 * oc_nc_stage_h1$alpha))
cat(sprintf("  E(N | H1)          = %.0f patients\n",
            oc_nc_stage_h1$expected_sample_size))

## ----compare-stage------------------------------------------------------------
data.frame(
  prior     = c("Beta(1,1) flat", "logit-normal mixture"),
  type_I_pp = c(100 * oc_conj_stage_h0$alpha, 100 * oc_nc_stage_h0$alpha),
  power_pp  = c(100 * oc_conj_stage_h1$alpha, 100 * oc_nc_stage_h1$alpha),
  E_N_H1    = c(oc_conj_stage_h1$expected_sample_size,
                oc_nc_stage_h1$expected_sample_size)
)

