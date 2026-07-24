## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)
set.seed(1L)

## ----design, message = FALSE--------------------------------------------------
library(adabay)

des <- set_design(
  endpoint     = "binary",
  n_per_look   = c(760, 1520, 2280, 3040, 3800),
  effect_scale = "risk_difference",
  alternative  = "less"
)

dec <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = FALSE
)

## ----conj-prior---------------------------------------------------------------
pri_conj <- set_prior(
  endpoint = "binary",
  arms = list(c = list(family = "beta", a = 1, b = 1),
              t = list(family = "beta", a = 1, b = 1))
)
pri_conj

## ----conj-eval, message = FALSE-----------------------------------------------
oc_conj_h0 <- evaluate_design(
  des, pri_conj, dec,
  effect   = list(theta_c = 0.33, theta_t = 0.33),
  n_trials = 2000, cores = 1, seed = 1L
)
oc_conj_h1 <- evaluate_design(
  des, pri_conj, dec,
  effect   = list(theta_c = 0.33, theta_t = 0.28),
  n_trials = 2000, cores = 1, seed = 1L
)

cat(sprintf("Conjugate Beta(1,1) prior:\n"))
cat(sprintf("  Type I error rate  = %.2f%%\n", 100 * oc_conj_h0$alpha))
cat(sprintf("  Power              = %.2f%%\n", 100 * oc_conj_h1$eff_prob))
cat(sprintf("  E(N | H1)          = %.0f patients\n",
            oc_conj_h1$expected_sample_size))

## ----conj-plot, fig.alt = "Per-look efficacy and futility stopping probabilities under H1 for the conjugate-prior design."----
if (requireNamespace("ggplot2", quietly = TRUE)) plot(oc_conj_h1)

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

## ----nonconj-overlay, fig.alt = "Empirical logit-normal prior density overlaid on the fitted beta-mixture density, with tail-quantile rugs."----
if (requireNamespace("ggplot2", quietly = TRUE)) plot(pri_c_arm)

## ----nonconj-eval, message = FALSE--------------------------------------------
oc_nc_h0 <- evaluate_design(
  des, pri_nc, dec,
  effect   = list(theta_c = 0.33, theta_t = 0.33),
  n_trials = 2000, cores = 1, seed = 1L
)
oc_nc_h1 <- evaluate_design(
  des, pri_nc, dec,
  effect   = list(theta_c = 0.33, theta_t = 0.28),
  n_trials = 2000, cores = 1, seed = 1L
)

cat(sprintf("Non-conjugate (logit-normal) prior via beta mixture:\n"))
cat(sprintf("  Type I error rate  = %.2f%%\n", 100 * oc_nc_h0$alpha))
cat(sprintf("  Power              = %.2f%%\n", 100 * oc_nc_h1$eff_prob))
cat(sprintf("  E(N | H1)          = %.0f patients\n",
            oc_nc_h1$expected_sample_size))

## ----compare------------------------------------------------------------------
data.frame(
  prior     = c("Beta(1,1) flat", "logit-normal mixture"),
  type_I_pp = c(100 * oc_conj_h0$alpha, 100 * oc_nc_h0$alpha),
  power_pp  = c(100 * oc_conj_h1$eff_prob, 100 * oc_nc_h1$eff_prob),
  E_N_H1    = c(oc_conj_h1$expected_sample_size,
                oc_nc_h1$expected_sample_size)
)

