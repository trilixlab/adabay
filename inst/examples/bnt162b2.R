## Pfizer-BioNTech BNT162b2 COVID-19 vaccine trial (count endpoint).
## Reproduces the count-endpoint case study of Section 5.3 of the manuscript.
## Original trial: Polack FP, Thomas SJ, Kitchin N, et al. NEJM 2020;383:2603-2615.
##
## n_trials = 1e6 (the high-precision manuscript budget) with cores = 8
## (the manuscript reference workstation has eight logical threads) so the
## example reproduces the values in Table 6 of the manuscript and
## Table S4 of the supplement. The run takes a few seconds on the
## reference workstation.

library(adabay)

des <- set_design(endpoint          = "count",
                  exposure_per_look = c(2540, 4920, 7300, 9520, 13020),
                  effect_scale      = "log_rate_ratio",
                  alternative       = "less",
                  delta_null        = log(0.70))

pri <- set_prior(endpoint = "count",
                 arms = list(c = list(family = "gamma", a = 1, b = 1),
                             t = list(family = "gamma", a = 1, b = 1)))

dec <- set_decision(
  efficacy = list(list(threshold_effect = log(0.70), threshold_prob = 0.99)),
  futility = list(list(threshold_effect = log(0.70), threshold_prob = 0.90)),
  futility_binding = TRUE)

## H0: rate ratio at the protocol superiority threshold (VE = 30%).
## H1: rate ratio at the designed-to-detect alternative   (VE = 60%).
oc_h0 <- evaluate_design(des, pri, dec,
                         effect = list(lambda_c = 0.018, lambda_t = 0.0126),
                         n_trials = 1e6, cores = 8, seed = 1L)
oc_h1 <- evaluate_design(des, pri, dec,
                         effect = list(lambda_c = 0.018, lambda_t = 0.0072),
                         n_trials = 1e6, cores = 8, seed = 1L)

print(oc_h0)
print(oc_h1)
