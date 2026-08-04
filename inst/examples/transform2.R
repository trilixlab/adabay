## TRANSFORM-2: esketamine nasal spray plus a newly initiated oral
## antidepressant vs placebo nasal spray plus a newly initiated oral
## antidepressant in adults with treatment-resistant depression.
## Continuous (MADRS change at Day 28) endpoint.
## Reproduces the continuous-endpoint case study of Section 5.1 of the manuscript.
## Original trial: Popova V, Daly EJ, Trivedi M, et al.
##                 Am J Psychiatry 2019;176:428-438.
##
## n_trials = 1e6 (the high-precision manuscript budget) with cores = 8
## (the manuscript reference workstation has eight logical threads) so the
## example reproduces the values in Table 4 of the manuscript and
## Table S2 of the supplement. The run takes under a minute on the
## reference workstation.

library(adabay)

des <- set_design(endpoint     = "continuous",
                  n_per_look   = c(90, 178, 266, 354, 444),
                  effect_scale = "mean_difference",
                  alternative  = "less",
                  sigma        = 13)

pri <- set_prior(endpoint = "continuous",
                 arms = list(c = list(family = "normal", mean = -10, sd = 1e3),
                             t = list(family = "normal", mean = -10, sd = 1e3)))

dec <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE)

## Control: placebo nasal spray + oral antidepressant
##   (TRANSFORM-2 placebo MADRS change at Day 28: -15.8).
## Treatment: esketamine nasal spray + oral antidepressant
##   (TRANSFORM-2 esketamine MADRS change at Day 28: -19.8).
oc_h0 <- evaluate_design(des, pri, dec,
                         effect = list(mu_c = -15.8, mu_t = -15.8),
                         n_trials = 1e6, cores = 8, seed = 1L)
oc_h1 <- evaluate_design(des, pri, dec,
                         effect = list(mu_c = -15.8, mu_t = -19.8),
                         n_trials = 1e6, cores = 8, seed = 1L)

print(oc_h0)
print(oc_h1)
