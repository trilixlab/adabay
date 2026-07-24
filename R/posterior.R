## Closed-form per-look posterior computation.
##
## Each function returns the right-tail posterior probability
##   P(Delta > e | data)
## for a vector of effect-size thresholds e under the user-specified prior
## (single conjugate component or finite mixture). The posterior of Delta is
## a mixture of L_c x L_t conjugate components, and the tail probability is
## the convex combination of per-component tail probabilities (Eq. S1 of the
## supplement).
##
## Slow inner integrals are dispatched to the C++ implementations declared in
## src/posterior.cpp via Rcpp::sourceCpp; the R-level code below is otherwise
## fully vectorised across the R virtual trials.
##
## All functions are internal.

## Continuous endpoint -------------------------------------------------------

## Update one normal-known-variance component, vectorised over R trials.
.update_normal_known <- function(comp, sigma2, n, ybar) {
  rho2      <- comp$sd^2
  prec_post <- 1 / rho2 + n / sigma2
  rho2_post <- 1 / prec_post
  list(mean = rho2_post * (comp$mean / rho2 + n * ybar / sigma2),
       sd   = sqrt(rho2_post))
}

## Update one normal--inverse-gamma component, vectorised over R trials.
.update_nig <- function(comp, n, ybar, S2) {
  k_post  <- comp$kappa + n
  list(nu    = (comp$kappa * comp$nu + n * ybar) / k_post,
       kappa = k_post,
       alpha = comp$alpha + n / 2,
       beta  = comp$beta + 0.5 * S2 +
         0.5 * comp$kappa * n / k_post * (ybar - comp$nu)^2)
}

## Mixture-component posterior reweighting -----------------------------------
##
## Under a finite conjugate-mixture prior the posterior mixture weights are the
## prior weights reweighted by each component's marginal likelihood of the data,
##   w*_l proportional to w_l * L_l(D_k)                 (supplement Eq. S11),
## and, the two arms being independent, the joint posterior weight factorises as
## w*_{ij} = w*_{c,i} w*_{t,j}                           (supplement Eq. S12).
## Each helper returns the per-component log marginal likelihood, vectorised over
## the R virtual trials, up to an additive constant common across an arm's
## components (such constants cancel in the within-arm normalisation below).

.logmarg_normal_known <- function(comp, sigma2, n, ybar) {
  stats::dnorm(ybar, mean = comp$mean,
               sd = sqrt(comp$sd^2 + sigma2 / n), log = TRUE)
}

.logmarg_nig <- function(comp, n, ybar, S2) {
  k_post <- comp$kappa + n
  a_post <- comp$alpha + n / 2
  b_post <- comp$beta + 0.5 * S2 +
    0.5 * comp$kappa * n / k_post * (ybar - comp$nu)^2
  0.5 * (log(comp$kappa) - log(k_post)) +
    comp$alpha * log(comp$beta) - a_post * log(b_post) +
    lgamma(a_post) - lgamma(comp$alpha)
}

.logmarg_beta <- function(comp, n, S) {
  lbeta(comp$a + S, comp$b + n - S) - lbeta(comp$a, comp$b)
}

.logmarg_gamma <- function(comp, E, S) {
  comp$a * log(comp$b) - (comp$a + S) * log(comp$b + E) +
    lgamma(comp$a + S) - lgamma(comp$a)
}

## Convert an R x L matrix of (log prior weight + log marginal likelihood) into
## an R x L matrix of posterior component weights summing to one across columns.
## A single component returns a column of ones, so single-component conjugate
## priors are bit-identical to the pre-reweighting behaviour of v0.1.0.
.posterior_component_weights <- function(logw) {
  if (ncol(logw) == 1L) return(matrix(1, nrow = nrow(logw), ncol = 1L))
  mx <- logw[cbind(seq_len(nrow(logw)), max.col(logw, ties.method = "first"))]
  w  <- exp(logw - mx)
  w / rowSums(w)
}

## Posterior component-weight matrix (R x L) for one arm, from its prior weights
## and the R x L matrix of per-component log marginal likelihoods.
.arm_weights <- function(weights, logmarg) {
  .posterior_component_weights(sweep(logmarg, 2L, log(weights), `+`))
}

## Posterior P(Delta > e) for the continuous endpoint.
##
## Returns a matrix nrow = R, ncol = length(thresholds).
.tail_prob_continuous <- function(prior, sigma, scale = "mean_difference",
                                  n_c, n_t, ybar_c, ybar_t,
                                  S2_c = NULL, S2_t = NULL,
                                  thresholds) {
  arms   <- prior$arms
  family <- arms$c$family
  R      <- length(ybar_c)
  out    <- matrix(0, nrow = R, ncol = length(thresholds))
  ## On the standardised-mean-difference scale the threshold e is in units of
  ## sigma, so P((mu_t - mu_c)/sigma > e) = P(mu_t - mu_c > e * sigma).
  if (identical(scale, "standardised_mean_difference")) {
    if (is.null(sigma))
      stop("Effect scale 'standardised_mean_difference' requires 'sigma' in set_design().",
           call. = FALSE)
    thresholds <- thresholds * sigma
  }
  ## Posterior mixture weights, reweighted per trial by the marginal likelihood.
  if (family == "normal") {
    lm_c <- vapply(arms$c$components,
                   function(cc) .logmarg_normal_known(cc, sigma^2, n_c, ybar_c),
                   numeric(R))
    lm_t <- vapply(arms$t$components,
                   function(ct) .logmarg_normal_known(ct, sigma^2, n_t, ybar_t),
                   numeric(R))
  } else if (family == "normal_inverse_gamma") {
    lm_c <- vapply(arms$c$components,
                   function(cc) .logmarg_nig(cc, n_c, ybar_c, S2_c), numeric(R))
    lm_t <- vapply(arms$t$components,
                   function(ct) .logmarg_nig(ct, n_t, ybar_t, S2_t), numeric(R))
  } else {
    stop("Unsupported continuous prior family.", call. = FALSE)
  }
  Wc <- .arm_weights(arms$c$weights, matrix(lm_c, nrow = R))
  Wt <- .arm_weights(arms$t$weights, matrix(lm_t, nrow = R))
  for (i in seq_along(arms$c$weights)) {
    for (j in seq_along(arms$t$weights)) {
      w  <- Wc[, i] * Wt[, j]
      cc <- arms$c$components[[i]]
      ct <- arms$t$components[[j]]
      if (family == "normal") {
        post_c <- .update_normal_known(cc, sigma^2, n_c, ybar_c)
        post_t <- .update_normal_known(ct, sigma^2, n_t, ybar_t)
        mu_d <- post_t$mean - post_c$mean
        sd_d <- sqrt(post_t$sd^2 + post_c$sd^2)
        ## Vectorise over thresholds: matrix of P(Delta > e_k) per (trial, k).
        for (k in seq_along(thresholds)) {
          out[, k] <- out[, k] +
            w * stats::pnorm(thresholds[k], mu_d, sd_d, lower.tail = FALSE)
        }
      } else {
        post_c <- .update_nig(cc, n_c, ybar_c, S2_c)
        post_t <- .update_nig(ct, n_t, ybar_t, S2_t)
        for (k in seq_along(thresholds)) {
          out[, k] <- out[, k] +
            w * .pdiff_t(thresholds[k],
                         loc1 = post_t$nu,
                         scale1 = sqrt(post_t$beta /
                                         (post_t$alpha * post_t$kappa)),
                         df1 = 2 * post_t$alpha,
                         loc2 = post_c$nu,
                         scale2 = sqrt(post_c$beta /
                                         (post_c$alpha * post_c$kappa)),
                         df2 = 2 * post_c$alpha)
        }
      }
    }
  }
  out
}

## P(T1 - T2 > e) for two independent non-standardised Student-t random
## variables, where T_a ~ t_{df_a}(loc_a, scale_a) means
##   T_a = loc_a + scale_a * W_a,  W_a ~ standard t_{df_a}.
##
## NOTE: this path uses stats::integrate at the R level (no C++ dispatch)
## and is the only inner integral that is not implemented in posterior.cpp.
## It is invoked only on the continuous endpoint with the
## normal-inverse-gamma prior; for the much more common normal known-
## variance prior the closed-form .update_normal_known() path is used.
##
## We compute
##   P(T1 - T2 > e)
##     = integral over s of P(T1 > e + s) * f_{T2}(s) ds
##     = integral over u of pt( (scale2 * u + loc2 - loc1 + e) / scale1; df1,
##                             lower.tail = FALSE )
##                       * dt(u; df2) du,
## where u = (s - loc2) / scale2 is the standardised T2. This integrates over
## the standard-t support [-10, 10], so the bounds do not depend on loc2 or
## scale2.
##
## Vectorised across R trials. Each of (loc1, scale1, df1, loc2, scale2,
## df2) may be a scalar or a length-R vector; scalars are recycled to length
## R, which is needed because the normal--inverse-gamma update produces a
## scalar \tilde\alpha and \tilde\kappa shared across all trials, while
## \tilde\nu and \tilde\beta vary per trial.
.pdiff_t <- function(e, loc1, scale1, df1, loc2, scale2, df2) {
  R <- max(length(loc1), length(scale1), length(df1),
           length(loc2), length(scale2), length(df2))
  loc1   <- rep_len(loc1,   R); scale1 <- rep_len(scale1, R)
  df1    <- rep_len(df1,    R)
  loc2   <- rep_len(loc2,   R); scale2 <- rep_len(scale2, R)
  df2    <- rep_len(df2,    R)
  vapply(seq_len(R), function(i) {
    fn <- function(u) {
      x <- (scale2[i] * u + loc2[i] - loc1[i] + e) / scale1[i]
      stats::pt(x, df = df1[i], lower.tail = FALSE) *
        stats::dt(u, df = df2[i])
    }
    tryCatch(stats::integrate(fn, lower = -10, upper = 10)$value,
             error = function(ex) NA_real_)
  }, numeric(1))
}

## Binary endpoint ----------------------------------------------------------

## Posterior P(Delta > e) for the binary endpoint on the requested scale.
.tail_prob_binary <- function(prior, scale,
                              n_c, n_t, S_c, S_t,
                              thresholds) {
  arms <- prior$arms
  R    <- length(S_c)
  out  <- matrix(0, nrow = R, ncol = length(thresholds))
  lm_c <- vapply(arms$c$components,
                 function(cc) .logmarg_beta(cc, n_c, S_c), numeric(R))
  lm_t <- vapply(arms$t$components,
                 function(ct) .logmarg_beta(ct, n_t, S_t), numeric(R))
  Wc <- .arm_weights(arms$c$weights, matrix(lm_c, nrow = R))
  Wt <- .arm_weights(arms$t$weights, matrix(lm_t, nrow = R))
  for (i in seq_along(arms$c$weights)) {
    for (j in seq_along(arms$t$weights)) {
      w   <- Wc[, i] * Wt[, j]
      cc  <- arms$c$components[[i]]
      ct  <- arms$t$components[[j]]
      a_c <- cc$a + S_c
      b_c <- cc$b + n_c - S_c
      a_t <- ct$a + S_t
      b_t <- ct$b + n_t - S_t
      for (k in seq_along(thresholds)) {
        out[, k] <- out[, k] + w * .pdelta_binary(thresholds[k], scale,
                                                  a_c, b_c, a_t, b_t)
      }
    }
  }
  out
}

## Vectorised dispatch on scale; calls C++ when available, otherwise R fallback.
.pdelta_binary <- function(e, scale, a_c, b_c, a_t, b_t) {
  if (scale == "risk_difference") {
    return(pdelta_binary_rd_cpp(e, a_c, b_c, a_t, b_t))
  }
  if (scale == "risk_ratio") {
    return(pdelta_binary_rr_cpp(e, a_c, b_c, a_t, b_t))
  }
  if (scale == "odds_ratio") {
    return(pdelta_binary_or_cpp(e, a_c, b_c, a_t, b_t))
  }
  stop("Unsupported binary scale.", call. = FALSE)
}

## Count and tte endpoints --------------------------------------------

## Posterior P(Delta > e) for count and tte endpoints.
.tail_prob_gamma <- function(prior, scale,
                             E_c, E_t, S_c, S_t,
                             thresholds) {
  arms <- prior$arms
  R    <- length(S_c)
  out  <- matrix(0, nrow = R, ncol = length(thresholds))
  lm_c <- vapply(arms$c$components,
                 function(cc) .logmarg_gamma(cc, E_c, S_c), numeric(R))
  lm_t <- vapply(arms$t$components,
                 function(ct) .logmarg_gamma(ct, E_t, S_t), numeric(R))
  Wc <- .arm_weights(arms$c$weights, matrix(lm_c, nrow = R))
  Wt <- .arm_weights(arms$t$weights, matrix(lm_t, nrow = R))
  for (i in seq_along(arms$c$weights)) {
    for (j in seq_along(arms$t$weights)) {
      w   <- Wc[, i] * Wt[, j]
      cc  <- arms$c$components[[i]]
      ct  <- arms$t$components[[j]]
      a_c <- cc$a + S_c
      ## Count exposure is deterministic (length 1); recycle to length R so the
      ## C++ tail routine never indexes past the end of a scalar (bug B2).
      b_c <- rep_len(cc$b + E_c, R)
      a_t <- ct$a + S_t
      b_t <- rep_len(ct$b + E_t, R)
      for (k in seq_along(thresholds)) {
        out[, k] <- out[, k] + w * .pdelta_gamma(thresholds[k], scale,
                                                 a_c, b_c, a_t, b_t)
      }
    }
  }
  out
}

## Vectorised gamma-tail dispatch.
.pdelta_gamma <- function(e, scale, a_c, b_c, a_t, b_t) {
  if (scale %in% c("rate_ratio", "hazard_ratio")) {
    if (e <= 0) return(rep(1, length(a_c)))
    z <- e / (e + b_c / b_t)
    return(stats::pbeta(z, a_t, a_c, lower.tail = FALSE))
  }
  if (scale %in% c("log_rate_ratio", "log_hazard_ratio")) {
    z <- exp(e) / (exp(e) + b_c / b_t)
    return(stats::pbeta(z, a_t, a_c, lower.tail = FALSE))
  }
  if (scale == "rate_difference") {
    return(pdelta_gamma_rd_cpp(e, a_c, b_c, a_t, b_t))
  }
  stop("Unsupported gamma scale.", call. = FALSE)
}
