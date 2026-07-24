## CheckMate-141: nivolumab vs investigator's choice chemotherapy in
## recurrent/metastatic head and neck squamous cell carcinoma (HNSCC).
## Time-to-event (tte; overall survival) endpoint.
## Reproduces the time-to-event case study of Section 5.4 of the manuscript.
## Original trial: Ferris RL, Blumenschein G, Fayette J, et al.
##                 NEJM 2016;375:1856-1867.
##
## n_trials = 1e6 (the high-precision manuscript budget) with cores = 8
## (the manuscript reference workstation has eight logical threads) so the
## example reproduces the values in Tables M.7.4 and S.3.4 of the
## manuscript. The run takes a couple of minutes on the reference workstation.
## Time-to-event (tte) designs are always staggered: an accrual specification is
## mandatory, because the Poisson arrival times drive the calendar-time
## exposure cutoff (see set_accrual() below).
##
## CheckMate-141 was 2:1 randomised (nivolumab:chemo); for design-stage
## exposition we re-implement under a 1:1 allocation. Pass
## allocation_ratio = 2 to set_design() to recover the original 2:1 split.

library(adabay)

des <- set_design(endpoint     = "tte",
                  d_total      = 250,
                  d_per_look   = c(50, 100, 150, 200, 250),
                  effect_scale = "log_hazard_ratio",
                  alternative  = "less")

pri <- set_prior(endpoint = "tte",
                 arms = list(c = list(family = "gamma", a = 1, b = 1),
                             t = list(family = "gamma", a = 1, b = 1)))

dec <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE)

## Recruitment rate back-calculated from CheckMate-141: 361 patients were
## randomised between 29 May 2014 and 31 July 2015 (about 14.1 months), i.e. an
## accrual rate of roughly 25.6 patients per month overall, pooled across both
## arms (Ferris et al. NEJM 2016; ClinicalTrials.gov NCT02105636). set_accrual()
## always takes the pooled rate; it is split into per-arm rates by
## set_design()'s allocation_ratio (1:1 here, so ~12.8 patients per month per
## arm, matching the manuscript's design-stage assumption).
acc <- set_accrual(model = "poisson", rate = 25.6)

## Control: median OS 5.1 months (CheckMate-141 chemotherapy arm).
## Treatment: median OS 7.5 months (CheckMate-141 nivolumab arm).
oc_h0 <- evaluate_design(des, pri, dec,
                         effect  = list(lambda_c = log(2)/5.1,
                                        lambda_t = log(2)/5.1),
                         accrual = acc,
                         n_trials = 1e6, cores = 8, seed = 1L)
oc_h1 <- evaluate_design(des, pri, dec,
                         effect  = list(lambda_c = log(2)/5.1,
                                        lambda_t = log(2)/7.5),
                         accrual = acc,
                         n_trials = 1e6, cores = 8, seed = 1L)

print(oc_h0)
print(oc_h1)
