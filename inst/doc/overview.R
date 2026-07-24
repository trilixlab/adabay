## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)
set.seed(1L)

## ----mwe, message = FALSE-----------------------------------------------------
library(adabay)

des <- set_design(
  endpoint     = "binary",
  n_per_look   = c(760, 1520, 2280, 3040, 3800),
  effect_scale = "risk_difference",
  alternative  = "less"
)

pri <- set_prior(
  endpoint = "binary",
  arms     = list(c = list(family = "beta", a = 1, b = 1),
                  t = list(family = "beta", a = 1, b = 1))
)

dec <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = FALSE
)

oc_h1 <- evaluate_design(
  des, pri, dec,
  effect   = list(theta_c = 0.33, theta_t = 0.28),
  n_trials = 2000, cores = 1, seed = 1L
)
oc_h1

## ----plot-oc, fig.alt = "Per-look stopping probabilities for the ADRENAL re-design under H1."----
if (requireNamespace("ggplot2", quietly = TRUE)) plot(oc_h1)

