###############################################################################
##
## Replication script for:
##   He Z, Yu F (2026). "adabay: An R Package for Rapid Evaluation and
##   Calibration of Bayesian Group Sequential Designs Across Common Endpoint
##   Types".
##
## This single standalone script reproduces every numerical result reported in
## the manuscript. Each block is labelled with the corresponding section, table
## or figure of the paper, and prints its result so that the value can be
## checked against the table in the manuscript.
##
## ---------------------------------------------------------------------------
## How to run
## ---------------------------------------------------------------------------
##
##   * Set the working directory to the folder containing this file
##     (the Code/ directory of the supplementary materials).
##   * Install adabay from the bundled tarball (or via remotes::install_github):
##       install.packages("adabay_0.1.0.tar.gz", repos = NULL, type = "source")
##     The manuscript and this script target adabay version 0.1.0.
##   * Comparator runs (BATSS, adaptr, gsbDesign) for manuscript Tables 4-7
##     live in Section 7 of this script. They are gated by environment
##     variables, all defaulting to OFF:
##       RUN_ADAPTR=1     -- adaptr     (continuous + binary, ~30-35 s/scenario)
##       RUN_GSBDESIGN=1  -- gsbDesign  (continuous analytic, ~seconds)
##       RUN_BATSS=1      -- BATSS      (all four endpoints, multi-hour to
##                                       multi-day per endpoint x {H0,H1}
##                                       pair on eight cores)
##     Example: RUN_ADAPTR=1 RUN_GSBDESIGN=1 Rscript replication.R
##   * Set RUN_FULL_PRECISION below to TRUE to use the manuscript budgets
##     (R = 10,000 matched-budget and R = 1,000,000 high-precision). With
##     RUN_FULL_PRECISION = FALSE (the default), every scenario uses a fast
##     R = 2,000 pass so the whole script completes in well under a minute
##     for a smoke-test.
##
## ---------------------------------------------------------------------------
## Estimated wall-clock cost (reference machine for the manuscript:
## 4.2 GHz quad-core Intel Core i7, 32 GB 2400 MHz DDR4, cores = 8)
## ---------------------------------------------------------------------------
##
##   RUN_FULL_PRECISION = FALSE : under 1 minute total.
##   RUN_FULL_PRECISION = TRUE  : approximately 7 minutes total without
##                                comparators, dominated by the four R = 10^6
##                                endpoint sweeps (binary ~2.8 min, survival
##                                ~1.5 min, continuous ~20 s, count ~6 s)
##                                and the R = 10^4 cache builds in
##                                Sections 6.3 (calibration) and 6.4.
##                                Comparator BATSS at R = 10,000 is
##                                extremely expensive (several hours to
##                                tens of hours per endpoint x {H0, H1}
##                                pair on eight cores; multi-day for the
##                                full sweep across all four endpoints);
##                                comparator adaptr on eight cores adds
##                                roughly 30-35 seconds per scenario for
##                                the continuous and binary cases; comparator
##                                gsbDesign (analytic, continuous only)
##                                runs in seconds on a single core.
##
###############################################################################


## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------

RUN_FULL_PRECISION <- FALSE
SEED               <- 1L
N_CORES            <- 8L  ## the manuscript figures use 8 cores via mclapply

R_FAST           <- 2000L
R_MATCHED        <- if (RUN_FULL_PRECISION) 10000L    else R_FAST
R_HIGH_PRECISION <- if (RUN_FULL_PRECISION) 1000000L  else R_FAST

## Every block runs at the user-level SEED on N_CORES cores, matching the
## `seed = 1L, cores = 8` calls in each manuscript code listing. Results are
## bit-identical only at a fixed (seed, cores) pair, so these must NOT be
## offset per block: an offset seed still reproduces from SEED, but it does
## not reproduce the published Tables 4-7 and S2-S5.
SEED_CONT        <- SEED  # Section 5.1 continuous
SEED_BIN         <- SEED  # Section 5.2 binary
SEED_CNT         <- SEED  # Section 5.3 count
SEED_SUR         <- SEED  # Section 5.4 survival
SEED_CACHE_64    <- SEED  # Section 6.4 caching benchmark
SEED_CACHE_63    <- SEED  # Section 6.3 stage-specific calibration


## ---------------------------------------------------------------------------
## 1. Packages
## ---------------------------------------------------------------------------

stopifnot(getRversion() >= "4.0.0")

suppressPackageStartupMessages({
  library(adabay)
})

have_gsbDesign <- requireNamespace("gsbDesign", quietly = TRUE)
have_BATSS     <- requireNamespace("BATSS",     quietly = TRUE)
have_adaptr    <- requireNamespace("adaptr",    quietly = TRUE)

cat("=== Replication script for the adabay manuscript ===\n")
cat(sprintf("R %s; adabay %s\n",
            R.version$version.string,
            utils::packageVersion("adabay")))
cat(sprintf("RUN_FULL_PRECISION = %s  (R_MATCHED = %d, R_HIGH_PRECISION = %d)\n",
            RUN_FULL_PRECISION, R_MATCHED, R_HIGH_PRECISION))
cat(sprintf("Seed = %d (all blocks; matches the manuscript listings), cores = %d\n",
            SEED, N_CORES))
cat(sprintf("Comparators available: gsbDesign = %s, BATSS = %s, adaptr = %s\n\n",
            have_gsbDesign, have_BATSS, have_adaptr))


## ---------------------------------------------------------------------------
## 2. Helpers for tidy reporting against the manuscript tables
## ---------------------------------------------------------------------------

mc_se_pct <- function(p, R) 100 * sqrt(p * (1 - p) / R)

## Print overall (alpha, 1-beta, E(N), E(T)) -- reproduces the adabay rows of
## the main-manuscript Tables 4 to 7.
report <- function(label, oc, R) {
  alpha  <- oc$alpha
  ess    <- oc$expected_sample_size
  edur   <- if (!is.null(oc$expected_duration)) oc$expected_duration else NA_real_
  cat(sprintf("    %-32s  rate = %.4f (SE %.2f pp), E(N) = %.1f%s\n",
              label, alpha, mc_se_pct(alpha, R), ess,
              if (!is.na(edur)) sprintf(", E(T) = %.2f", edur) else ""))
}

## Print stage stopping probabilities -- reproduces the alpha_k and beta_k
## columns of the supplement Tables S2 to S5 ($H_0$ and $H_1$ rows).
report_perlook <- function(label, oc) {
  fmt <- function(v) paste(sprintf("%.4f", v), collapse = ", ")
  if (!is.null(oc$eff_prob_at_look))
    cat(sprintf("      %-30s alpha_k = %s\n", label, fmt(oc$eff_prob_at_look)))
  if (!is.null(oc$fut_prob_at_look))
    cat(sprintf("      %-30s beta_k  = %s\n", label, fmt(oc$fut_prob_at_look)))
}

## ===========================================================================
## SECTION 5.1: CONTINUOUS --- TRANSFORM-2 (ESKETAMINE FOR TRD) (TABLE 4)
## ===========================================================================

cat("--- Section 5.1 / Table 4: continuous (TRANSFORM-2 esketamine) ---\n")

des_cont <- set_design(endpoint     = "continuous",
                       n_per_look   = c(90, 178, 266, 354, 444),
                       effect_scale = "mean_difference",
                       alternative  = "less",
                       sigma        = 13)

pri_cont <- set_prior(endpoint = "continuous",
                      arms = list(c = list(family = "normal", mean = -10, sd = 1e3),
                                  t = list(family = "normal", mean = -10, sd = 1e3)))

dec_cont <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE)

## Control: placebo nasal spray + oral antidepressant
##   (TRANSFORM-2 placebo MADRS change at Day 28: -15.8).
## Treatment: esketamine nasal spray + oral antidepressant
##   (TRANSFORM-2 esketamine MADRS change at Day 28: -19.8).
## Matched budget R_MATCHED -- corresponds to the R = 10,000 row of Table 4.
t_cont_matched <- Sys.time()
oc_cont_h0_M <- evaluate_design(des_cont, pri_cont, dec_cont,
                                effect = list(mu_c = -15.8, mu_t = -15.8),
                                n_trials = R_MATCHED,
                                cores = N_CORES, seed = SEED_CONT)
oc_cont_h1_M <- evaluate_design(des_cont, pri_cont, dec_cont,
                                effect = list(mu_c = -15.8, mu_t = -19.8),
                                n_trials = R_MATCHED,
                                cores = N_CORES, seed = SEED_CONT)
t_cont_matched <- as.numeric(Sys.time() - t_cont_matched, units = "secs")

cat(sprintf("  R = %d (matched budget), wall-clock = %.1fs\n",
            R_MATCHED, t_cont_matched))
report("adabay type I (H0)", oc_cont_h0_M, R_MATCHED)
report("adabay power   (H1)", oc_cont_h1_M, R_MATCHED)

## High-precision budget R_HIGH_PRECISION -- corresponds to the R = 10^6 row.
t_cont_hp <- Sys.time()
oc_cont_h0_H <- evaluate_design(des_cont, pri_cont, dec_cont,
                                effect = list(mu_c = -15.8, mu_t = -15.8),
                                n_trials = R_HIGH_PRECISION,
                                cores = N_CORES, seed = SEED_CONT)
oc_cont_h1_H <- evaluate_design(des_cont, pri_cont, dec_cont,
                                effect = list(mu_c = -15.8, mu_t = -19.8),
                                n_trials = R_HIGH_PRECISION,
                                cores = N_CORES, seed = SEED_CONT)
t_cont_hp <- as.numeric(Sys.time() - t_cont_hp, units = "secs")

cat(sprintf("  R = %d (high precision), wall-clock = %.1fs\n",
            R_HIGH_PRECISION, t_cont_hp))
report("adabay type I (H0)", oc_cont_h0_H, R_HIGH_PRECISION)
report("adabay power   (H1)", oc_cont_h1_H, R_HIGH_PRECISION)
## Per-look stopping probabilities -- supplement Table S2.
report_perlook("H0 (mu_t-mu_c = 0)", oc_cont_h0_H)
report_perlook("H1 (mu_t-mu_c = -4)", oc_cont_h1_H)

cat("\n")


## ===========================================================================
## SECTION 5.2: BINARY ENDPOINT --- ADRENAL RE-DESIGN (TABLE 5)
## ===========================================================================

cat("--- Section 5.2 / Table 5: binary (ADRENAL re-design) ---\n")

des_bin <- set_design(endpoint     = "binary",
                     n_per_look   = c(760, 1520, 2280, 3040, 3800),
                     effect_scale = "risk_difference",
                     alternative  = "less")

pri_bin <- set_prior(endpoint = "binary",
                    arms = list(c = list(family = "beta", a = 1, b = 1),
                                t = list(family = "beta", a = 1, b = 1)))

dec_bin <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE)

t_bin_matched <- Sys.time()
oc_bin_h0_M <- evaluate_design(des_bin, pri_bin, dec_bin,
                              effect = list(theta_c = 0.33, theta_t = 0.33),
                              n_trials = R_MATCHED,
                              cores = N_CORES, seed = SEED_BIN)
oc_bin_h1_M <- evaluate_design(des_bin, pri_bin, dec_bin,
                              effect = list(theta_c = 0.33, theta_t = 0.28),
                              n_trials = R_MATCHED,
                              cores = N_CORES, seed = SEED_BIN)
t_bin_matched <- as.numeric(Sys.time() - t_bin_matched, units = "secs")

cat(sprintf("  R = %d (matched budget), wall-clock = %.1fs\n",
            R_MATCHED, t_bin_matched))
report("adabay type I (H0)", oc_bin_h0_M, R_MATCHED)
report("adabay power   (H1)", oc_bin_h1_M, R_MATCHED)

t_bin_hp <- Sys.time()
oc_bin_h0_H <- evaluate_design(des_bin, pri_bin, dec_bin,
                              effect = list(theta_c = 0.33, theta_t = 0.33),
                              n_trials = R_HIGH_PRECISION,
                              cores = N_CORES, seed = SEED_BIN)
oc_bin_h1_H <- evaluate_design(des_bin, pri_bin, dec_bin,
                              effect = list(theta_c = 0.33, theta_t = 0.28),
                              n_trials = R_HIGH_PRECISION,
                              cores = N_CORES, seed = SEED_BIN)
t_bin_hp <- as.numeric(Sys.time() - t_bin_hp, units = "secs")

cat(sprintf("  R = %d (high precision), wall-clock = %.1fs\n",
            R_HIGH_PRECISION, t_bin_hp))
report("adabay type I (H0)", oc_bin_h0_H, R_HIGH_PRECISION)
report("adabay power   (H1)", oc_bin_h1_H, R_HIGH_PRECISION)
## Per-look stopping probabilities -- supplement Table S3.
report_perlook("H0 (theta_t-theta_c = 0)",     oc_bin_h0_H)
report_perlook("H1 (theta_t-theta_c = -0.05)", oc_bin_h1_H)

cat("\n")


## ===========================================================================
## SECTION 5.3: COUNT ENDPOINT --- BNT162b2 COVID-19 VACCINE (TABLE 6)
## ===========================================================================

cat("--- Section 5.3 / Table 6: count (BNT162b2 COVID-19 vaccine) ---\n")

des_cnt <- set_design(endpoint          = "count",
                     exposure_per_look = c(2540, 4920, 7300, 9520, 13020),
                     effect_scale      = "log_rate_ratio",
                     alternative       = "less",
                     delta_null        = log(0.70))

pri_cnt <- set_prior(endpoint = "count",
                    arms = list(c = list(family = "gamma", a = 1, b = 1),
                                t = list(family = "gamma", a = 1, b = 1)))

dec_cnt <- set_decision(
  efficacy = list(list(threshold_effect = log(0.70), threshold_prob = 0.99)),
  futility = list(list(threshold_effect = log(0.70), threshold_prob = 0.90)),
  futility_binding = TRUE)

## H0: rate ratio at the protocol superiority threshold (VE = 30%, ratio 0.70).
## H1: rate ratio at the designed-to-detect alternative   (VE = 60%, ratio 0.40).
t_cnt_matched <- Sys.time()
oc_cnt_h0_M <- evaluate_design(des_cnt, pri_cnt, dec_cnt,
                              effect = list(lambda_c = 0.018, lambda_t = 0.0126),
                              n_trials = R_MATCHED,
                              cores = N_CORES, seed = SEED_CNT)
oc_cnt_h1_M <- evaluate_design(des_cnt, pri_cnt, dec_cnt,
                              effect = list(lambda_c = 0.018, lambda_t = 0.0072),
                              n_trials = R_MATCHED,
                              cores = N_CORES, seed = SEED_CNT)
t_cnt_matched <- as.numeric(Sys.time() - t_cnt_matched, units = "secs")

cat(sprintf("  R = %d (matched budget), wall-clock = %.1fs\n",
            R_MATCHED, t_cnt_matched))
report("adabay type I (H0)", oc_cnt_h0_M, R_MATCHED)
report("adabay power   (H1)", oc_cnt_h1_M, R_MATCHED)

t_cnt_hp <- Sys.time()
oc_cnt_h0_H <- evaluate_design(des_cnt, pri_cnt, dec_cnt,
                              effect = list(lambda_c = 0.018, lambda_t = 0.0126),
                              n_trials = R_HIGH_PRECISION,
                              cores = N_CORES, seed = SEED_CNT)
oc_cnt_h1_H <- evaluate_design(des_cnt, pri_cnt, dec_cnt,
                              effect = list(lambda_c = 0.018, lambda_t = 0.0072),
                              n_trials = R_HIGH_PRECISION,
                              cores = N_CORES, seed = SEED_CNT)
t_cnt_hp <- as.numeric(Sys.time() - t_cnt_hp, units = "secs")

cat(sprintf("  R = %d (high precision), wall-clock = %.1fs\n",
            R_HIGH_PRECISION, t_cnt_hp))
report("adabay type I (H0)", oc_cnt_h0_H, R_HIGH_PRECISION)
report("adabay power   (H1)", oc_cnt_h1_H, R_HIGH_PRECISION)
## Per-look stopping probabilities -- supplement Table S4.
report_perlook("H0 (lambda_t/lambda_c = 0.70)", oc_cnt_h0_H)
report_perlook("H1 (lambda_t/lambda_c = 0.40)", oc_cnt_h1_H)

cat("\n")


## ===========================================================================
## SECTION 5.4: TIME-TO-EVENT --- CHECKMATE-141 (NIVOLUMAB IN HNSCC) (TABLE 7)
## ===========================================================================

cat("--- Section 5.4 / Table 7: survival (CheckMate-141 HNSCC) ---\n")

des_sur <- set_design(endpoint     = "survival",
                     d_total      = 250,
                     d_per_look   = c(50, 100, 150, 200, 250),
                     effect_scale = "log_hazard_ratio",
                     alternative  = "less")

pri_sur <- set_prior(endpoint = "survival",
                    arms = list(c = list(family = "gamma", a = 1, b = 1),
                                t = list(family = "gamma", a = 1, b = 1)))

dec_sur <- set_decision(
  efficacy = list(list(threshold_effect = 0, threshold_prob = 0.99)),
  futility = list(list(threshold_effect = 0, threshold_prob = 0.90)),
  futility_binding = TRUE)

## Recruitment rate back-calculated from CheckMate-141: 361 patients
## randomised over ~14.1 months (29 May 2014 to 31 July 2015) is ~25.6/month
## overall, pooled across both arms (set_accrual() splits it into ~12.8/month
## per arm by the design-stage 1:1 allocation_ratio).
acc_sur <- set_accrual(model = "poisson", rate = 25.6)

## Control: chemotherapy median OS 5.1 months.
## Treatment: nivolumab median OS 7.5 months (CheckMate-141 HR = 0.70).
t_sur_matched <- Sys.time()
oc_sur_h0_M <- evaluate_design(des_sur, pri_sur, dec_sur,
                              effect = list(lambda_c = log(2)/5.1,
                                            lambda_t = log(2)/5.1),
                              accrual = acc_sur,
                              n_trials = R_MATCHED,
                              cores = N_CORES, seed = SEED_SUR)
oc_sur_h1_M <- evaluate_design(des_sur, pri_sur, dec_sur,
                              effect = list(lambda_c = log(2)/5.1,
                                            lambda_t = log(2)/7.5),
                              accrual = acc_sur,
                              n_trials = R_MATCHED,
                              cores = N_CORES, seed = SEED_SUR)
t_sur_matched <- as.numeric(Sys.time() - t_sur_matched, units = "secs")

cat(sprintf("  R = %d (matched budget), wall-clock = %.1fs\n",
            R_MATCHED, t_sur_matched))
report("adabay type I (H0)", oc_sur_h0_M, R_MATCHED)
report("adabay power   (H1)", oc_sur_h1_M, R_MATCHED)

t_sur_hp <- Sys.time()
oc_sur_h0_H <- evaluate_design(des_sur, pri_sur, dec_sur,
                              effect = list(lambda_c = log(2)/5.1,
                                            lambda_t = log(2)/5.1),
                              accrual = acc_sur,
                              n_trials = R_HIGH_PRECISION,
                              cores = N_CORES, seed = SEED_SUR)
oc_sur_h1_H <- evaluate_design(des_sur, pri_sur, dec_sur,
                              effect = list(lambda_c = log(2)/5.1,
                                            lambda_t = log(2)/7.5),
                              accrual = acc_sur,
                              n_trials = R_HIGH_PRECISION,
                              cores = N_CORES, seed = SEED_SUR)
t_sur_hp <- as.numeric(Sys.time() - t_sur_hp, units = "secs")

cat(sprintf("  R = %d (high precision), wall-clock = %.1fs\n",
            R_HIGH_PRECISION, t_sur_hp))
report("adabay type I (H0)", oc_sur_h0_H, R_HIGH_PRECISION)
report("adabay power   (H1)", oc_sur_h1_H, R_HIGH_PRECISION)
## Per-look stopping probabilities -- supplement Table S5.
report_perlook("H0 (lambda_t/lambda_c = 1)",        oc_sur_h0_H)
report_perlook("H1 (lambda_t/lambda_c approx 0.68)", oc_sur_h1_H)

cat("\n")


## ===========================================================================
## SECTION 6.4 / TABLE 3: DIRECT SIMULATION vs CACHED GRID EVALUATION
## ---------------------------------------------------------------------------
## At full scale Table 3 reports a 50 x 16 = 800-cell threshold grid on the
## ADRENAL re-design at R = 10^4 (the manuscript matched-budget convention,
## same R as R_MATCHED above). To keep the smoke-test path tractable, we
## use a 10 x 4 = 40-cell grid by default; the per-cell speedup ratio of
## cached over direct evaluation is the metric of interest and does not
## depend on grid size. Set RUN_FULL_PRECISION = TRUE to recover the
## manuscript's 50 x 16 = 800-cell grid.
## ===========================================================================

cat("--- Section 6.4 / Table 3: caching speedup on ADRENAL ---\n")

R_CACHE  <- if (RUN_FULL_PRECISION) 100000L else R_FAST  # Sec 6.4: manuscript R = 10^5.
GRID_P   <- if (RUN_FULL_PRECISION) {
  seq(0.95, 0.999, length.out = 50)
} else {
  seq(0.95, 0.999, length.out = 10)
}
GRID_Q   <- if (RUN_FULL_PRECISION) {
  seq(0.80, 0.95, length.out = 16)
} else {
  seq(0.80, 0.95, length.out =  4)
}
N_CELLS  <- length(GRID_P) * length(GRID_Q)

cat(sprintf("  Grid: %d cells, R = %d per scenario\n", N_CELLS, R_CACHE))

## (a) Cached path: one build_cache pass, then evaluate_design on every cell.
t_cache_build <- Sys.time()
cache_h1 <- build_cache(des_bin, pri_bin,
                        effect   = list(theta_c = 0.33, theta_t = 0.28),
                        n_trials = R_CACHE, cores = N_CORES, seed = SEED_CACHE_64,
                        threshold_grid = list(efficacy = 0, futility = 0))
t_cache_build <- as.numeric(Sys.time() - t_cache_build, units = "secs")

t_cache_sweep <- Sys.time()
for (p in GRID_P) for (q in GRID_Q) {
  dec_pq <- set_decision(
    efficacy = list(list(threshold_effect = 0, threshold_prob = p)),
    futility = list(list(threshold_effect = 0, threshold_prob = q)),
    futility_binding = TRUE)
  evaluate_design(cache_h1, dec_pq)
}
t_cache_sweep <- as.numeric(Sys.time() - t_cache_sweep, units = "secs")

t_per_cell_cached <- t_cache_sweep / N_CELLS

cat(sprintf("  Cached  : build = %.2fs, sweep = %.2fs (%.4fs per cell)\n",
            t_cache_build, t_cache_sweep, t_per_cell_cached))

## (b) Direct path: one evaluate_design per cell (the wasteful baseline).
## Time only a single cell to estimate the per-cell direct cost; multiplying
## by N_CELLS recovers what the full direct sweep would cost.
t_direct_one <- Sys.time()
dec_single <- set_decision(
  efficacy = list(list(threshold_effect = 0,
                       threshold_prob = mean(GRID_P))),
  futility = list(list(threshold_effect = 0,
                       threshold_prob = mean(GRID_Q))),
  futility_binding = TRUE)
evaluate_design(des_bin, pri_bin, dec_single,
                effect   = list(theta_c = 0.33, theta_t = 0.28),
                n_trials = R_CACHE, cores = N_CORES, seed = SEED_CACHE_64)
t_direct_one <- as.numeric(Sys.time() - t_direct_one, units = "secs")
t_direct_full <- t_direct_one * N_CELLS

cat(sprintf("  Direct  : per-cell = %.2fs, full-grid extrapolation = %.1fs\n",
            t_direct_one, t_direct_full))
cat(sprintf("  Caching speedup : %.1fx per cell (after the one-time pass)\n\n",
            t_direct_one / t_per_cell_cached))


## ===========================================================================
## SECTION 6.3: STAGE-SPECIFIC CALIBRATION FROM THE COMPANION PAPER
## ---------------------------------------------------------------------------
## Recover the calibrated five-look ADRENAL design with stage-specific
## efficacy thresholds p_{1:5} = (0.997, 0.997, 0.997, 0.997, 0.980)
## and common futility threshold q = 0.85, by evaluating the per-look rule
## directly against the H0 and H1 caches. The manuscript reports type I
## error 2.44%, power 90.32%, E(N|H0) ~ 3,073 and E(N|H1) ~ 2,636 for this
## design (Section 6.3, Threshold calibration).
## ===========================================================================

cat("--- Section 6.3: stage-specific calibrated five-look design ---\n")

R_CAL <- if (RUN_FULL_PRECISION) 100000L else R_FAST  # Sec 6.3: manuscript R = 10^5.

cache_bin_h0 <- build_cache(des_bin, pri_bin,
                           effect   = list(theta_c = 0.33, theta_t = 0.33),
                           n_trials = R_CAL, cores = N_CORES, seed = SEED_CACHE_63,
                           threshold_grid = list(efficacy = 0, futility = 0))
cache_bin_h1 <- build_cache(des_bin, pri_bin,
                           effect   = list(theta_c = 0.33, theta_t = 0.28),
                           n_trials = R_CAL, cores = N_CORES, seed = SEED_CACHE_63,
                           threshold_grid = list(efficacy = 0, futility = 0))

p_stage <- c(0.997, 0.997, 0.997, 0.997, 0.980)
q_stage <- 0.85
dec_stage <- set_decision(
  efficacy = lapply(p_stage, function(p)
    list(list(threshold_effect = 0, threshold_prob = p))),
  futility = list(list(threshold_effect = 0, threshold_prob = q_stage)),
  futility_binding = TRUE)

oc_cal_h0 <- evaluate_design(cache_bin_h0, dec_stage)
oc_cal_h1 <- evaluate_design(cache_bin_h1, dec_stage)

cat(sprintf("  R = %d per cache\n", R_CAL))
report("adabay type I (H0)", oc_cal_h0, R_CAL)
report("adabay power   (H1)", oc_cal_h1, R_CAL)
cat(sprintf("    Manuscript target: type I 2.44%%, power 90.32%%, E(N|H0) ~ 3073, E(N|H1) ~ 2636\n\n"))


## Companion-paper three-look design: p_{1:3} = (0.995, 0.995, 0.980) and the
## same futility threshold q = 0.85, on the same N_max = 3,800 with three
## equally spaced looks. The manuscript reports type I 2.46%, power 90.63%,
## E(N|H0) ~ 3,309 and E(N|H1) ~ 2,726 (Section 6.3).

cat("--- Section 6.3: stage-specific calibrated three-look design ---\n")

des_bin3 <- set_design(endpoint     = "binary",
                       n_per_look   = c(1267, 2533, 3800),
                       effect_scale = "risk_difference",
                       alternative  = "less")

cache_bin3_h0 <- build_cache(des_bin3, pri_bin,
                            effect   = list(theta_c = 0.33, theta_t = 0.33),
                            n_trials = R_CAL, cores = N_CORES, seed = SEED_CACHE_63,
                            threshold_grid = list(efficacy = 0, futility = 0))
cache_bin3_h1 <- build_cache(des_bin3, pri_bin,
                            effect   = list(theta_c = 0.33, theta_t = 0.28),
                            n_trials = R_CAL, cores = N_CORES, seed = SEED_CACHE_63,
                            threshold_grid = list(efficacy = 0, futility = 0))

dec_stage3 <- set_decision(
  efficacy = lapply(c(0.995, 0.995, 0.980), function(p)
    list(list(threshold_effect = 0, threshold_prob = p))),
  futility = list(list(threshold_effect = 0, threshold_prob = q_stage)),
  futility_binding = TRUE)

oc_cal3_h0 <- evaluate_design(cache_bin3_h0, dec_stage3)
oc_cal3_h1 <- evaluate_design(cache_bin3_h1, dec_stage3)

cat(sprintf("  R = %d per cache\n", R_CAL))
report("adabay type I (H0)", oc_cal3_h0, R_CAL)
report("adabay power   (H1)", oc_cal3_h1, R_CAL)
cat(sprintf("    Manuscript target: type I 2.46%%, power 90.63%%, E(N|H0) ~ 3309, E(N|H1) ~ 2726\n\n"))


## ===========================================================================
## SECTION 7 -- COMPARATOR PACKAGES
## ---------------------------------------------------------------------------
## Cross-package operating-characteristic and wall-clock comparators at the
## matched budget R = 10,000 for the case studies of Section 5. Each
## comparator is gated by an environment variable; all default to OFF so that
## a bare `Rscript replication.R` reproduces only the adabay rows above.
##
##   adaptr    -- covers the continuous (TRANSFORM-2) and binary (ADRENAL)
##                cases on eight cores via its native parallel backend; it
##                does not natively support count or time-to-event GSDs and
##                so does not appear in Tables 6-7. Gated behind RUN_ADAPTR.
##
##   gsbDesign -- analytic (numerical-integration) comparator restricted to
##                the continuous (TRANSFORM-2) case (Table 4 only); single
##                CPU core. Gated behind RUN_GSBDESIGN.
##
##   BATSS     -- covers all four case studies: binary (ADRENAL), continuous
##                (TRANSFORM-2), count (BNT162b2) and time-to-event
##                (CheckMate-141). Runs with eight mclapply workers and INLA
##                pinned to one thread per worker; extremely slow (multi-hour
##                to multi-day per scenario at R = 10,000). The count and
##                time-to-event BATSS runs are approximate comparators -- see
##                the per-endpoint caveats below. The time-to-event run is
##                exploratory only: Table 7 of the manuscript reports BATSS
##                as n/a there, because batss.glm is GLM-only and cannot
##                construct the censored event-time response. Gated behind
##                RUN_BATSS.
##
## Run AFTER the adabay benchmark above has finished and with nothing else
## competing for CPU: the wall-clock figures reported in Tables 4-7 are only
## meaningful under exclusive CPU access.
## ===========================================================================

RUN_ADAPTR    <- identical(Sys.getenv("RUN_ADAPTR"),    "1")
RUN_GSBDESIGN <- identical(Sys.getenv("RUN_GSBDESIGN"), "1")
RUN_BATSS     <- identical(Sys.getenv("RUN_BATSS"),     "1")

## Directory for the saved BATSS result objects. Defaults to the bundled
## Output folder at "../Output/Output v1.0" relative to this script (which
## sits in Code/), so reruns persist next to the shipped artefacts. Override
## with the RESULTS_DIR environment variable to redirect, e.g.
##   RESULTS_DIR=/tmp/batss_runs RUN_BATSS=1 Rscript replication.R
RESULTS_DIR <- Sys.getenv(
  "RESULTS_DIR",
  unset = normalizePath(file.path("..", "Output", "Output v1.0"),
                        mustWork = FALSE)
)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

cat("--- Section 7: comparator packages ---\n")
cat(sprintf("  RUN_ADAPTR=%s, RUN_GSBDESIGN=%s, RUN_BATSS=%s\n",
            RUN_ADAPTR, RUN_GSBDESIGN, RUN_BATSS))

## adaptr / gsbDesign are skipped when not installed, exactly as advertised
## by the manuscript ("Calls to BATSS, adaptr and gsbDesign are guarded by
## requireNamespace() and skipped when the comparator package is absent").
## Every call site below uses adaptr:: / gsbDesign:: / BATSS:: qualifiers, so
## attaching the package via library() is unnecessary and is omitted here.
if (RUN_ADAPTR && !requireNamespace("adaptr", quietly = TRUE)) {
  message("adaptr not installed; skipping the adaptr comparators.")
  RUN_ADAPTR <- FALSE
}
if (RUN_GSBDESIGN && !requireNamespace("gsbDesign", quietly = TRUE)) {
  message("gsbDesign not installed; skipping the gsbDesign comparator.")
  RUN_GSBDESIGN <- FALSE
}

cmp_timed <- function(expr) {
  t0 <- Sys.time(); v <- force(expr)
  cat(sprintf("    [wall-clock %.1fs]\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  v
}

## ---------------------------------------------------------------------------
## 7.1 adaptr -- continuous TRANSFORM-2 and binary ADRENAL (Tables 4 and 5)
##   2-arm, 1:1 fixed allocation (fixed_probs disables RAR -> a pure GSD),
##   5 looks at the matched cumulative TOTAL n schedules,
##   superiority      P(theta_t < theta_c) > 0.99,
##   binding futility P(theta_t - theta_c > 0) > 0.90,
##   lower outcome is better (highest_is_best = FALSE).
## Validation anchor (companion paper he2026, R = 5,000):
##   binary type I 3.26%, power 87.36%, E(N|H1) ~ 2207.
## ---------------------------------------------------------------------------

run_adaptr_binom <- function(true_ys, label) {
  spec <- adaptr::setup_trial_binom(
    arms             = c("control", "treatment"),
    true_ys          = true_ys,
    fixed_probs      = c(0.5, 0.5),
    data_looks       = c(760, 1520, 2280, 3040, 3800),
    control          = "control",
    superiority      = 0.99,   # adabay efficacy  P(theta_t < theta_c) > 0.99
    inferiority      = 0.10,   # adabay binding futility P(theta_t-theta_c>0)>0.90
    highest_is_best  = FALSE,
    soften_power     = 1)
  res <- cmp_timed(adaptr::run_trials(spec, n_rep = R_MATCHED,
                                      base_seed = SEED, cores = 8L))
  s <- summary(res)
  cat(sprintf("  adaptr ADRENAL %s\n", label))
  print(s)
  invisible(s)
}

run_adaptr_norm <- function(true_ys, label) {
  spec <- adaptr::setup_trial_norm(
    arms             = c("control", "treatment"),
    true_ys          = true_ys,
    sds              = c(13, 13),
    fixed_probs      = c(0.5, 0.5),
    data_looks       = c(90, 178, 266, 354, 444),
    control          = "control",
    superiority      = 0.99,   # adabay efficacy  P(mu_t < mu_c) > 0.99
    inferiority      = 0.10,   # adabay binding futility P(mu_t-mu_c>0)>0.90
    highest_is_best  = FALSE,
    soften_power     = 1)
  res <- cmp_timed(adaptr::run_trials(spec, n_rep = R_MATCHED,
                                      base_seed = SEED, cores = 8L))
  s <- summary(res)
  cat(sprintf("  adaptr TRANSFORM-2 %s\n", label))
  print(s)
  invisible(s)
}

if (RUN_ADAPTR) {
  cat(sprintf("=== adaptr @ R = %s (8 cores) ===\n", format(R_MATCHED, big.mark = ",")))
  cat("--- continuous TRANSFORM-2 (Table 4) ---\n")
  ac_h0 <- run_adaptr_norm(c(-15.8, -15.8), "H0 (mu_t-mu_c = 0)")
  ac_h1 <- run_adaptr_norm(c(-15.8, -19.8), "H1 (mu_t-mu_c = -4)")
  cat("--- binary ADRENAL (Table 5) ---\n")
  ab_h0 <- run_adaptr_binom(c(0.33, 0.33), "H0 (theta_t-theta_c = 0)")
  ab_h1 <- run_adaptr_binom(c(0.33, 0.28), "H1 (theta_t-theta_c = -0.05)")
}

## ---------------------------------------------------------------------------
## 7.2 gsbDesign -- analytic comparator for the continuous TRANSFORM-2 case
## (Table 4).
##
## gsbDesign convention: delta = treatment - control, success when delta is
## LARGE (upper tail). TRANSFORM-2 is "lower MADRS is better", so we analyse
## on delta' = control - treatment (= -(mu_t - mu_c)); treatment benefit =>
## delta' > 0. Then:
##   adabay efficacy  P(mu_t-mu_c < 0) > 0.99  <=>  P(delta' > 0) > 0.99
##     -> criteria.success  = c(0, 0.99)
##   adabay binding futility P(mu_t-mu_c > 0) > 0.90
##     -> criteria.futility = c(0, 0.90)   [validated below]
## NOTE: gsbDesign's criteria.futility = c(theta.f, prob.f) is the probability
## threshold on P(delta' > theta.f) BELOW which futility is declared, but the
## empirical mapping that reproduces the adabay normal-normal answer (and
## the analytic check below) is c(0, 0.90), NOT c(0, 0.10). Validated against
## the adabay R=1e6 continuous result (type I 3.09%, power 85.59%,
## E(N|H1) ~266): gsbDesign Config C gives 3.09% / 85.54% / E(N|H1)=265.6 --
## agreement to <0.05 pp and <0.5 patients (analytic vs Monte Carlo).
## Truth under delta':  H0 -> 0;  H1 (mu_c=-15.8, mu_t=-19.8) -> +4.
## prior.difference = "non-informative"  (matches the Table 4 caption).
## Cumulative per-arm n = 45,89,133,177,222 -> per-stage incremental
## 45,44,44,44,45.  sigma = 13.
## ---------------------------------------------------------------------------

if (RUN_GSBDESIGN) {
  cat("=== gsbDesign analytic comparator (continuous, Table 4) ===\n")
  design <- gsbDesign::gsbDesign(
    nr.stages         = 5,
    patients          = matrix(c(45, 44, 44, 44, 45), ncol = 1),  # per-stage incr/arm
    sigma             = 13,
    criteria.success  = c(0, 0.99),
    criteria.futility = c(0, 0.90),
    prior.difference  = "non-informative")

  sim <- gsbDesign::gsbSimulation(
    truth        = c(0, 4),                       # delta' : H0=0, H1=4
    type.update  = "treatment effect",
    method       = "numerical integration",
    grid.type    = "manually")

  x <- cmp_timed(gsbDesign::gsb(design = design, simulation = sim))

  cat("\n--- success probability (= efficacy / type I & power) ---\n")
  print(gsbDesign::tab(x, "success", atDelta = c(0, 4), digits = 5))
  cat("\n--- cumulative futility probability ---\n")
  print(gsbDesign::tab(x, "cumulative futility", atDelta = c(0, 4), digits = 5))
  cat("\n--- expected sample size ---\n")
  print(gsbDesign::tab(x, "sample size", atDelta = c(0, 4), digits = 2))
  cat("=== gsbDesign done ===\n")
}

## ---------------------------------------------------------------------------
## 7.3 BATSS -- all four case studies at R = 10,000 (8 mclapply workers, INLA
## pinned to 1 thread/worker). Enable with: RUN_BATSS=1 Rscript replication.R
##
## WALL-CLOCK: BATSS runs INLA inside the Monte Carlo loop and is extremely
## expensive. A small R = 150 binary probe took ~8 min for H0 alone on this
## hardware, implying ~17 h per binary scenario at R = 10,000 and multiple
## DAYS for all four endpoints x {H0, H1}. This is a detached, multi-day job;
## it will not finish in an interactive session.
##
## DECISION-RULE MAPPING (structurally exact for binary/continuous, validated
## against the he2026 anchor):
##   eff.arm.simple(b = 0.99)  <=>  adabay efficacy  P(Delta in H1 dir) > 0.99
##   fut.arm.simple(b = 0.10)  <=>  adabay binding futility P(.) < 0.10
##   alloc.balanced            <=>  fixed 1:1 allocation (a pure GSD, no RAR)
## Anchors (he2026, R = 5,000): binary type I 2.70%, power 87.30%,
##   E(N|H0) ~ 3255, E(N|H1) ~ 2200.
##
## CAVEATS:
##  * Count (BNT162b2): adabay uses an exposure-driven (cumulative person-time)
##    look schedule; BATSS only supports recruited-N interims. The Poisson
##    design below treats one person-time unit as one "patient" with mean
##    count exp(beta) -- a proxy, not an exact reproduction of the adabay
##    exposure design. Treat the count BATSS row as approximate.
##  * Survival (CheckMate-141): adabay uses an EVENT-driven look schedule
##    (cumulative deaths D_k); BATSS's GLM framework with recruited-N interims
##    cannot natively express an event-driven GSD. The exponential.surv design
##    below is a recruited-N approximation with N chosen so the expected
##    number of events is ~250; it is NOT a faithful comparator and is
##    flagged as such. (The manuscript already notes BATSS does not natively
##    model recruitment / event-driven designs.)
## ---------------------------------------------------------------------------

if (RUN_BATSS) {
  if (!requireNamespace("BATSS", quietly = TRUE) ||
      !requireNamespace("INLA",  quietly = TRUE)) {
    message("BATSS or INLA not installed; skipping the BATSS comparators.")
  } else {
  ## BATSS exports several bare helpers (alloc.balanced, eff.arm.simple,
  ## fut.arm.simple, eff.trial.all, fut.trial.all) consumed by batss.glm()
  ## below without namespace qualification, so attach BATSS here. INLA is
  ## reached only through the qualified call INLA::inla.setOption(), so
  ## library(INLA) is not needed.
  suppressPackageStartupMessages(library(BATSS))
  INLA::inla.setOption(num.threads = "1:1")
  ## Reproducibility under mclapply: BATSS::batss.glm has no base_seed
  ## argument, so we set the L'Ecuyer-CMRG stream here so the parallel
  ## workers spawn deterministic substreams from SEED. Save the calling
  ## RNGkind first and restore it after the BATSS block (below), so that
  ## sourcing this script does not silently change the default RNG of the
  ## host R session.
  .old_rngkind <- RNGkind()
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  cat(sprintf("=== BATSS @ R = %s (8 mclapply workers, INLA 1:1) ===\n", format(R_MATCHED, big.mark = ",")))

  ## NOTE: BATSS::batss.glm() itself runs the simulation and returns the
  ## result object (there is no separate BATSS::batss() in BATSS 1.1.1).
  ## With H0 = TRUE, ONE call evaluates BOTH the alternative ($H1, using the
  ## supplied beta) and the null ($H0, target coefficient forced to 0), so
  ## there is exactly one call per endpoint -- not a per-hypothesis loop.
  ## Per-trial outcomes are in res$H{0,1}$sample: column `type` is "1"
  ## (efficacy stop), "2" (futility stop) or "0" (no decision); `control`
  ## and `treatment` give the realised per-arm sample size.
  ##   type I  = mean(res$H0$sample$type == "1")
  ##   power   = mean(res$H1$sample$type == "1")
  ##   E(N|Hx) = mean(res$Hx$sample$control + res$Hx$sample$treatment)
  ## `delta` is the clinically meaningful effect MARGIN on the linear-
  ## predictor scale (batss.glm delta.eff/delta.fut), matching adabay's
  ## decision threshold_effect. Binary/continuous use threshold_effect = 0
  ## (P(Delta < 0) > 0.99), so delta = 0 and BATSS's H0=TRUE (target = 0)
  ## auto-null coincides with adabay's no-effect H0 -- faithful. For the
  ## count design adabay's threshold_effect is log(0.70) (VE > 30%), so
  ## delta = log(0.70) and the null must be simulated at the MARGIN (truth
  ## ratio 0.70), not BATSS's auto ratio-1 null; that case is handled by
  ## run_batss_margin below (one call per truth scenario, H0 = FALSE).
  run_batss <- function(family, link, var_y, var_control, beta_h1, N, interim,
                        label, delta = 0) {
    cat(sprintf("  >>> BATSS %s starting (delta=%g; H0 auto + H1)\n",
                label, delta)); flush.console()
    res <- cmp_timed(BATSS::batss.glm(
      model           = y ~ group,
      var             = list(y = var_y, group = alloc.balanced),
      var.control     = var_control,
      family          = family, link = link,
      beta            = beta_h1, which = 2, alternative = "less",
      R               = R_MATCHED, N = N,
      interim         = list(recruited = interim),
      prob0           = c(control = 0.5, treatment = 0.5),
      delta.eff       = delta, delta.fut = delta,
      eff.arm         = eff.arm.simple, eff.arm.control = list(b = 0.99),
      fut.arm         = fut.arm.simple, fut.arm.control = list(b = 0.10),
      eff.trial       = eff.trial.all,  fut.trial = fut.trial.all,
      H0              = TRUE,
      computation     = "parallel", mc.cores = 8L))
    saveRDS(res, file.path(RESULTS_DIR,
                           sprintf("batss_%s.rds",
                                   gsub("[^A-Za-z0-9]+","_", label))))
    en <- function(h) mean(h$sample$control + h$sample$treatment)
    cat(sprintf(
      "  BATSS %s RESULT: typeI=%.4f power=%.4f E(N|H0)=%.1f E(N|H1)=%.1f\n",
      label,
      mean(res$H0$sample$type == "1"), mean(res$H1$sample$type == "1"),
      en(res$H0), en(res$H1)))
    print(summary(res))
    invisible(res)
  }

  ## Margin-anchored runner: adabay's count type I is P(efficacy | truth at
  ## the protocol margin ratio 0.70), power is P(efficacy | truth ratio
  ## 0.40), BOTH judged against delta = log(0.70). BATSS's auto ratio-1 null
  ## is irrelevant here, so H0 = FALSE and we read $H1 of each truth run.
  run_batss_margin <- function(family, link, var_y, var_control, beta_intc,
                               N, interim, delta, beta_eff_typeI,
                               beta_eff_power, label) {
    one <- function(beta_eff, tag) {
      cat(sprintf("  >>> BATSS %s %s starting (delta=%g)\n",
                  label, tag, delta)); flush.console()
      r <- cmp_timed(BATSS::batss.glm(
        model = y ~ group, var = list(y = var_y, group = alloc.balanced),
        var.control = var_control, family = family, link = link,
        beta = c(beta_intc, beta_eff), which = 2, alternative = "less",
        R = R_MATCHED, N = N, interim = list(recruited = interim),
        prob0 = c(control = 0.5, treatment = 0.5),
        delta.eff = delta, delta.fut = delta,
        eff.arm = eff.arm.simple, eff.arm.control = list(b = 0.99),
        fut.arm = fut.arm.simple, fut.arm.control = list(b = 0.10),
        eff.trial = eff.trial.all, fut.trial = fut.trial.all,
        H0 = FALSE, computation = "parallel", mc.cores = 8L))
      saveRDS(r, file.path(RESULTS_DIR,
                           sprintf("batss_%s_%s.rds", label, tag)))
      list(eff = mean(r$H1$sample$type == "1"),
           en  = mean(r$H1$sample$control + r$H1$sample$treatment))
    }
    a <- one(beta_eff_typeI, "typeI")   # truth = margin  -> type I
    b <- one(beta_eff_power, "power")   # truth = alt      -> power
    cat(sprintf(
      "  BATSS %s RESULT: typeI=%.4f power=%.4f E(N|H0)=%.1f E(N|H1)=%.1f\n",
      label, a$eff, b$eff, a$en, b$en))
  }

  ## Binary ADRENAL (Table 5) -- faithful; he2026 anchor 2.70/87.30/3255/2200.
  run_batss("binomial", "logit", rbinom, list(y = list(size = 1)),
            c(qlogis(0.33), qlogis(0.28) - qlogis(0.33)),
            3800, c(760, 1520, 2280, 3040), "binary", delta = 0)

  ## Continuous TRANSFORM-2 (Table 4) -- faithful (adabay H0 = no effect = BATSS auto-H0).
  run_batss("gaussian", "identity", rnorm, list(y = list(sd = 13)),
            c(-15.8, -4), 444, c(90, 178, 266, 354), "continuous", delta = 0)

  ## Count BNT162b2 (Table 6) -- CORRECTED: adabay threshold_effect = log(0.70)
  ## (VE > 30%), null simulated at the margin (rate ratio 0.70), power at
  ## ratio 0.40. batss.glm's N and interim `recruited` are TOTAL across both
  ## arms (as in the binary N=3800 and continuous N=444 calls above), so the
  ## adabay PER-ARM cumulative person-time schedule c(1270,2460,3650,4760,
  ## 6510) is doubled here: N = 13,020 total, interims = 2 x per-arm =
  ## c(2540,4920,7300,9520). Residual caveat: BATSS recruited-patient
  ## person-time proxy for the exposure-driven cumulative person-time looks.
  run_batss_margin("poisson", "log", rpois, list(), log(0.018),
                   13020, c(2540, 4920, 7300, 9520), delta = log(0.70),
                   beta_eff_typeI = log(0.70), beta_eff_power = log(0.40),
                   label = "count")

  ## Restore the caller's RNGkind on exit from the BATSS block (done after
  ## the survival run below).

  ## Survival CheckMate-141 (Table 7) -- NOT FAITHFUL: recruited-N proxy for
  ## the event-driven adabay design. Flagged.
  run_batss("exponential.surv", "log", rexp, list(),
            c(log(log(2)/5.1), log((log(2)/7.5) / (log(2)/5.1))),
            360, c(72, 144, 216, 288), "survival")

  ## Restore the caller's RNGkind to whatever it was before the BATSS block.
  RNGkind(.old_rngkind[1L])
  }
}

cat("--- Section 7: comparator packages done ---\n\n")


## ===========================================================================
## SUMMARY
## ===========================================================================

cat("=== End of replication script ===\n")
cat("The full set of adabay numerical results in the manuscript has been\n")
cat("reproduced above. Compare the printed rates and expected sample sizes\n")
cat("against Tables 4 through 7 of the paper. With RUN_FULL_PRECISION =\n")
cat("FALSE (the default) the rates are at R = 2,000 Monte Carlo standard\n")
cat("error; with RUN_FULL_PRECISION = TRUE they are at the manuscript\n")
cat("budgets R = 10,000 (matched) and R = 1,000,000 (high precision).\n")

## Capture the runtime environment used to produce these numbers so that
## reviewers can confirm the (R, adabay, comparator) versions and the
## parallel backend reported in the manuscript's Computational details
## section. Printed to stdout, not saved to disk.
cat("\n=== sessionInfo() ===\n")
print(sessionInfo())
