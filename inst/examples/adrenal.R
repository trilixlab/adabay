## ADRENAL re-design (binary endpoint).
## Reproduces the binary-endpoint case study of Section 5.2 of the manuscript.
##
## n_trials = 1e6 (the high-precision manuscript budget) with cores = 4
## (the manuscript reference workstation has four physical cores) so the
## example reproduces the values in Tables M.7.2 and S.3.2 of the
## manuscript. The run takes a few minutes on the reference workstation.

library(adabay)

des <- set_design(endpoint     = "binary",
                  n_per_look   = c(760, 1520, 2280, 3040, 3800),
                  effect_scale = "risk_difference",
                  alternative  = "less")

pri <- set_prior(endpoint = "binary",
                 arms = list(c = list(family = "beta", a = 1, b = 1),
                             t = list(family = "beta", a = 1, b = 1)))

dec <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE)

oc_h0 <- evaluate_design(des, pri, dec,
                         effect = list(theta_c = 0.33, theta_t = 0.33),
                         n_trials = 1e6, cores = 4, seed = 1L)
oc_h1 <- evaluate_design(des, pri, dec,
                         effect = list(theta_c = 0.33, theta_t = 0.28),
                         n_trials = 1e6, cores = 4, seed = 1L)

print(oc_h0)
print(oc_h1)
