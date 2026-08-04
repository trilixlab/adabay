#' Approximate a prior by a finite mixture of conjugate components
#'
#' Approximates a user-supplied prior on an arm-specific parameter by a
#' finite mixture of conjugate components (beta, normal or gamma kernels)
#' fitted by minimisation of the empirical forward Kullback--Leibler
#' divergence. The number of components is selected by
#' \code{RBesT::automixfit()}'s AIC rule when the suggested \pkg{RBesT}
#' package is installed. Otherwise the internal fallback locates an elbow in
#' the divergence decrement -- the smallest \eqn{L} past which adding a
#' component no longer improves the divergence by at least \code{tol_kl} --
#' and then continues past that elbow one component at a time until the
#' maximum absolute tail-probability error falls within \code{tol_tail},
#' returning the smallest-error fit if none does. Either way the count is
#' capped at \code{n_components_max}.
#'
#' Forward KL is mass-covering: it pressures the approximation to put mass
#' wherever the true prior puts mass, but it does not, on its own, guarantee
#' accurate tail probabilities. Because Bayesian GSDs stop on posterior tail
#' probabilities, \code{fit_mixture()} additionally reports the empirical
#' tail-probability error of the fitted mixture at a set of quantile
#' thresholds (defaults to the 0.025, 0.05, 0.10, 0.50, 0.90, 0.95 and 0.975
#' empirical quantiles), and emits a \code{warning()} when the maximum
#' absolute tail-probability error exceeds \code{tol_tail}. The diagnostics
#' are stored on the returned object and summarised by the \code{print()}
#' and \code{plot()} methods.
#'
#' If the suggested package \pkg{RBesT} is available, this function delegates
#' to \code{RBesT::automixfit}, which provides robust kernel-specific EM
#' updates for the beta, normal and gamma kernels. If \pkg{RBesT} is not
#' available, a built-in EM fallback is used (currently for the beta kernel
#' only); other kernels return an informative error.
#'
#' @param endpoint One of \code{"binary"}, \code{"count"}, \code{"tte"}
#'   or \code{"continuous"}. The endpoint determines the kernel family.
#' @param arm Either \code{"c"} (control) or \code{"t"} (treatment), used as
#'   metadata on the returned object.
#' @param prior Either a function \code{prior(x)} returning the density at
#'   \code{x}, or a numeric vector of independent samples drawn from the
#'   prior. A density function is accepted only for the beta kernel (binary
#'   endpoints), where \code{n_samples} draws are produced by
#'   inverse-transform sampling on a fine grid over \eqn{[0, 1]}. For the
#'   gamma and normal kernels the support is unbounded and no fixed grid is
#'   reliable, so \code{prior} must be a numeric vector of samples; supplying
#'   a function raises an error.
#' @param n_components_max Maximum number of mixture components to try.
#' @param tol_kl Mixture-fit tolerance, supplied on the natural (untransformed)
#'   scale, e.g. \code{1e-3}. When \pkg{RBesT} is installed -- the default
#'   path -- it is forwarded to \code{RBesT::automixfit()} as the
#'   per-parameter EM convergence tolerance \code{eps}, which \pkg{RBesT}
#'   expresses on a \eqn{-\log_{10}} scale, so \code{tol_kl = 1e-3} is passed
#'   as \code{eps = 3}; the number of components is then chosen by the AIC
#'   rule inside \code{RBesT::automixfit()}. When \pkg{RBesT} is absent, the
#'   internal expectation-maximisation fallback uses \code{tol_kl} directly as
#'   the Kullback--Leibler decrement defining the elbow: the elbow is the
#'   smallest \eqn{L} past which adding a component no longer improves the
#'   divergence by at least \code{tol_kl}. That elbow is the starting point,
#'   not the answer, because the \code{tol_tail} stage may select a larger
#'   \eqn{L}.
#' @param tol_tail Tolerance for the maximum absolute tail-probability error
#'   across the default tail thresholds; \code{warning()} fires when the
#'   error exceeds this value. Default \code{5e-3} (0.5 percentage points).
#' @param tail_quantiles Numeric vector of quantile levels at which the
#'   tail-probability error is computed. Defaults to
#'   \code{c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)}.
#' @param n_samples Number of prior samples to use when \code{prior} is given
#'   as a function. Defaults to 10\eqn{^{4}}.
#' @param seed Optional integer seed. When \code{prior} is a density
#'   function, sample generation (\code{stats::runif()}-based inverse-
#'   transform sampling) draws from the ambient RNG stream by default, so
#'   the fitted mixture depends on whatever ran before this call; passing
#'   \code{seed} calls \code{set.seed(seed)} immediately beforehand so the
#'   result is reproducible regardless of prior RNG state or execution
#'   order. Has no effect when \code{prior} is already a numeric vector of
#'   samples (no randomness is drawn in that case). Default \code{NULL}
#'   leaves the RNG untouched (current behaviour).
#' @return An object of class \code{"adabay_prior_arm"} suitable for
#'   [set_prior()]; carries \code{family}, \code{components}, \code{weights}
#'   and a \code{diagnostics} attribute that includes the empirical samples
#'   used, the fitted forward KL, the per-quantile tail-probability error
#'   table, the maximum absolute tail error, a density overlay grid, the
#'   \code{seed} used (if any), and a logical \code{warning_fired} flag.
#' @examples
#' \donttest{
#'   ## A logit-normal prior on the control rate, approximated by a beta
#'   ## mixture. 'seed' makes the fit reproducible without relying on
#'   ## ambient RNG state (prior is a density function here, so samples are
#'   ## drawn internally).
#'   logit_normal_density <- function(theta) {
#'     z <- log(theta / (1 - theta))
#'     dnorm(z, mean = -0.7, sd = 0.4) / (theta * (1 - theta))
#'   }
#'   pri_c <- fit_mixture(endpoint = "binary", arm = "c",
#'                        prior = logit_normal_density,
#'                        n_components_max = 5, tol_kl = 1e-3, seed = 1L)
#'   print(pri_c)
#'   plot(pri_c)
#'   set_prior(endpoint = "binary",
#'             arms = list(c = pri_c,
#'                         t = list(family = "beta", a = 1, b = 1)))
#' }
#' @export
fit_mixture <- function(endpoint, arm,
                        prior, n_components_max = 5L, tol_kl = 1e-3,
                        tol_tail = 5e-3,
                        tail_quantiles = c(0.025, 0.05, 0.10, 0.50,
                                           0.90, 0.95, 0.975),
                        n_samples = 10000L,
                        seed = NULL) {
  endpoint <- .match_choice(endpoint, .adabay_endpoints)
  arm      <- .match_choice(arm, c("c", "t"))
  if (length(n_components_max) != 1L || !is.finite(n_components_max) ||
      n_components_max < 1L)
    stop("'n_components_max' must be a positive integer.", call. = FALSE)
  n_components_max <- as.integer(n_components_max)
  ## Validate before sort(): sort() drops NAs (na.last = NA), so an is.finite
  ## check afterwards can never fire, and an all-NA input would silently become
  ## a zero-length vector that passes every any() test.
  tail_quantiles <- as.numeric(tail_quantiles)
  if (length(tail_quantiles) == 0L || any(!is.finite(tail_quantiles)) ||
      any(tail_quantiles <= 0 | tail_quantiles >= 1))
    stop("'tail_quantiles' must be a non-empty numeric vector with every ",
         "value strictly in (0, 1).", call. = FALSE)
  tail_quantiles <- sort(unique(tail_quantiles))
  if (length(tol_kl) != 1L || !is.finite(tol_kl) || tol_kl <= 0)
    stop("'tol_kl' must be a positive scalar.", call. = FALSE)
  if (length(tol_tail) != 1L || !is.finite(tol_tail) || tol_tail <= 0)
    stop("'tol_tail' must be a positive scalar.", call. = FALSE)
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.finite(seed) || seed != round(seed)))
    stop("'seed' must be a single integer or NULL.", call. = FALSE)

  ## Determine kernel family.
  family <- switch(endpoint,
                   continuous = "normal",
                   binary     = "beta",
                   count      = "gamma",
                   tte        = "gamma")
  rbest_kernel <- switch(family,
                         beta   = "beta",
                         normal = "norm",
                         gamma  = "gamma")

  ## Resolve prior into samples. When 'prior' is a density function, samples
  ## are drawn by inverse-transform sampling, so a supplied seed is applied
  ## here to make that draw reproducible. Seeding here (rather than saving/
  ## restoring .Random.seed) matches the convention already used by
  ## evaluate_design()/build_cache(): a supplied seed deterministically sets
  ## the RNG stream from this point on, at the cost of also affecting
  ## subsequent unrelated draws in the caller's session. When 'prior' is a
  ## numeric vector no sampling occurs, so 'seed' has no effect.
  if (!is.null(seed) && is.function(prior)) set.seed(seed)
  samples <- .resolve_prior_samples(prior, n_samples = n_samples,
                                    family = family)

  ## Use RBesT if available; otherwise fall back to internal EM (beta only).
  if (requireNamespace("RBesT", quietly = TRUE)) {
    ## RBesT expresses its per-parameter EM convergence tolerance on a
    ## -log10 scale: EM_bmm_ab()/EM_gmm() rescale a length-1 'eps' as
    ## rep(10^(-eps), 3). Passing tol_kl directly would request a tolerance
    ## of 10^(-1e-3) ~ 0.998 rather than 1e-3, so convert here.
    fit <- RBesT::automixfit(sample = samples,
                             Nc = seq_len(n_components_max),
                             type = rbest_kernel,
                             eps = -log10(tol_kl))
    trace <- attr(fit, "traceMix")
    final <- trace[[length(trace)]]
    parsed <- .parse_rbest_mix(final, family)
    method <- "RBesT::automixfit"
    kl     <- .empirical_kl(samples, family, parsed$weights, parsed$components)
  } else {
    if (family != "beta")
      stop("RBesT is required for non-beta kernels. Install RBesT or supply a conjugate prior directly.",
           call. = FALSE)
    parsed <- .em_beta_mixture(samples,
                               n_components_max = n_components_max,
                               tol_kl = tol_kl,
                               tail_quantiles = tail_quantiles,
                               tol_tail = tol_tail)
    method <- "internal-EM"
    kl     <- parsed$kl
  }

  ## Tail-probability error diagnostics.
  tail_table <- .tail_error_table(samples, family,
                                  parsed$weights, parsed$components,
                                  tail_quantiles)
  max_tail   <- max(abs(tail_table$error))
  warning_fired <- max_tail > tol_tail
  if (warning_fired) {
    warning(sprintf(
      paste0("fit_mixture(): maximum absolute tail-probability error %.4f ",
             "exceeds tol_tail = %.4f. Consider increasing 'n_components_max', ",
             "providing more prior samples, or supplying a hand-chosen ",
             "conjugate-mixture prior directly. See attr(<fit>, \"diagnostics\")",
             " for details."),
      max_tail, tol_tail),
      call. = FALSE)
  }

  ## Density overlay grid (used by plot.adabay_prior_arm).
  overlay <- .density_overlay(samples, family,
                              parsed$weights, parsed$components)

  diagnostics <- list(
    L_hat        = length(parsed$weights),
    method       = method,
    kernel       = family,
    samples      = samples,
    kl           = kl,
    tail_quantiles = tail_quantiles,
    tail_table   = tail_table,
    max_tail_error = max_tail,
    tol_tail     = tol_tail,
    warning_fired = warning_fired,
    overlay      = overlay,
    seed         = seed)

  out <- list(family = family,
              weights = parsed$weights,
              components = parsed$components)
  attr(out, "diagnostics") <- diagnostics
  attr(out, "arm") <- arm
  class(out) <- c("adabay_prior_arm", "list")
  out
}

## Helpers -------------------------------------------------------------------

.resolve_prior_samples <- function(prior, n_samples, family) {
  if (is.numeric(prior)) {
    if (any(!is.finite(prior)))
      stop("'prior' samples must be finite.", call. = FALSE)
    return(as.numeric(prior))
  }
  if (!is.function(prior))
    stop("'prior' must be a numeric vector or a density function.",
         call. = FALSE)
  ## Inverse-transform sampling on a fine grid. Only the beta kernel has a
  ## bounded support on which a fixed grid is safe: on [0, 1] the 5000-node
  ## grid has spacing 2e-4. For the gamma and normal kernels any fixed range
  ## wide enough to be generally useful is far too coarse for the small rates
  ## typical of count and time-to-event designs (a 0-1000 range gives spacing
  ## 0.2, so a prior centred near 0.02 would collapse onto a single node and
  ## be returned, silently, as a degenerate sample). Rather than return a
  ## corrupted prior we require pre-drawn samples for those kernels.
  if (family != "beta")
    stop("For the '", family, "' kernel (endpoint ",
         if (family == "gamma") "\"count\" or \"tte\"" else "\"continuous\"",
         "), 'prior' must be supplied as a numeric vector of samples rather ",
         "than a density function: inverse-transform sampling on a fixed grid ",
         "is not reliable on an unbounded support. Draw from the prior ",
         "directly, e.g. prior = ",
         if (family == "gamma") "rgamma(1e4, shape, rate)." else "rnorm(1e4, mean, sd).",
         call. = FALSE)
  grid <- seq(1e-6, 1 - 1e-6, length.out = 5000)
  d    <- pmax(0, prior(grid))
  if (!any(d > 0))
    stop("Density evaluated to zero on the sampling grid.", call. = FALSE)
  cdf  <- cumsum(d) / sum(d)
  u    <- stats::runif(n_samples)
  stats::approx(x = cdf, y = grid, xout = u, ties = mean, rule = 2)$y
}

.parse_rbest_mix <- function(mix, family) {
  ## RBesT mixture matrices have rows: weight, then kernel-specific parameters.
  weights <- as.numeric(mix[1, ])
  if (family == "beta") {
    components <- lapply(seq_along(weights), function(j)
      list(a = as.numeric(mix[2, j]), b = as.numeric(mix[3, j])))
  } else if (family == "gamma") {
    components <- lapply(seq_along(weights), function(j)
      list(a = as.numeric(mix[2, j]), b = as.numeric(mix[3, j])))
  } else if (family == "normal") {
    components <- lapply(seq_along(weights), function(j)
      list(mean = as.numeric(mix[2, j]), sd = as.numeric(mix[3, j])))
  }
  list(weights = weights, components = components)
}

## Mixture-CDF and mixture-density helpers.
.mix_pdf <- function(x, family, weights, components) {
  vapply(x, function(z)
    sum(weights * vapply(components, function(c) {
      if (family == "beta")   stats::dbeta(z, c$a, c$b)
      else if (family == "gamma") stats::dgamma(z, shape = c$a, rate = c$b)
      else                    stats::dnorm(z, mean = c$mean, sd = c$sd)
    }, numeric(1))),
    numeric(1))
}

.mix_cdf <- function(x, family, weights, components) {
  vapply(x, function(z)
    sum(weights * vapply(components, function(c) {
      if (family == "beta")   stats::pbeta(z, c$a, c$b)
      else if (family == "gamma") stats::pgamma(z, shape = c$a, rate = c$b)
      else                    stats::pnorm(z, mean = c$mean, sd = c$sd)
    }, numeric(1))),
    numeric(1))
}

## Tail-error table: per-quantile P_true(<= q) (always equals the quantile
## level by construction of the empirical CDF) vs P_mixture(<= q).
.tail_error_table <- function(samples, family, weights, components,
                              tail_quantiles) {
  q_emp <- stats::quantile(samples, probs = tail_quantiles, names = FALSE,
                           type = 7L)
  p_emp <- tail_quantiles  # by definition of the empirical quantile
  p_mix <- .mix_cdf(q_emp, family, weights, components)
  data.frame(quantile_level = tail_quantiles,
             threshold      = q_emp,
             p_empirical    = p_emp,
             p_mixture      = p_mix,
             error          = p_mix - p_emp)
}

## Forward KL approximation (empirical-density vs mixture-density) on a
## kernel-appropriate grid.
.empirical_kl <- function(samples, family, weights, components) {
  rng <- switch(family,
                beta   = c(1e-4, 1 - 1e-4),
                normal = stats::quantile(samples, c(0.001, 0.999), names = FALSE),
                gamma  = c(1e-4, stats::quantile(samples, 0.999, names = FALSE)))
  grid <- seq(rng[1], rng[2], length.out = 1000L)
  pdf_mix <- .mix_pdf(grid, family, weights, components)
  kde     <- stats::density(samples,
                            from = rng[1], to = rng[2], n = 1000L)
  pdf_emp <- stats::approx(x = kde$x, y = kde$y, xout = grid, rule = 2)$y
  pdf_emp <- pmax(pdf_emp, 1e-12)
  pdf_mix <- pmax(pdf_mix, 1e-12)
  mean(pdf_emp * (log(pdf_emp) - log(pdf_mix))) * (max(grid) - min(grid))
}

## Density overlay grid for the plot method.
.density_overlay <- function(samples, family, weights, components) {
  rng <- switch(family,
                beta   = c(1e-4, 1 - 1e-4),
                normal = stats::quantile(samples, c(0.001, 0.999), names = FALSE),
                gamma  = c(1e-4, stats::quantile(samples, 0.999, names = FALSE)))
  grid <- seq(rng[1], rng[2], length.out = 512L)
  pdf_mix <- .mix_pdf(grid, family, weights, components)
  data.frame(x = grid, density_mixture = pdf_mix)
}

#' @noRd
#' @export
print.adabay_prior_arm <- function(x, ...) {
  d <- attr(x, "diagnostics")
  cat(sprintf("Fitted conjugate-mixture prior arm '%s'\n",
              if (!is.null(attr(x, "arm"))) attr(x, "arm") else "?"))
  cat(sprintf("  Kernel family:        %s\n", x$family))
  cat(sprintf("  Mixture components:   %d\n", length(x$weights)))
  for (i in seq_along(x$components)) {
    cat(sprintf("    [%d] weight = %.3f, %s\n", i, x$weights[i],
                paste(names(x$components[[i]]),
                      formatC(unlist(x$components[[i]]),
                              format = "g", digits = 4),
                      sep = " = ", collapse = ", ")))
  }
  cat(sprintf("  Method:               %s\n", d$method))
  cat(sprintf("  Forward KL (samples): %.5f\n", d$kl))
  cat(sprintf("  Max |tail error|:     %.4f (tol = %.4f)%s\n",
              d$max_tail_error, d$tol_tail,
              if (d$warning_fired) "  <-- exceeds tolerance" else ""))
  if (!is.null(d$seed))
    cat(sprintf("  Seed:                 %s\n", format(d$seed)))
  cat("  Per-quantile tail-probability table:\n")
  tt <- d$tail_table
  for (i in seq_len(nrow(tt))) {
    cat(sprintf("    q=%.3f  threshold=%.4g  P_empirical=%.4f  P_mixture=%.4f  error=%+.4f\n",
                tt$quantile_level[i], tt$threshold[i],
                tt$p_empirical[i], tt$p_mixture[i], tt$error[i]))
  }
  invisible(x)
}

#' @noRd
#' @export
plot.adabay_prior_arm <- function(x, file = NULL,
                               width = 7, height = 4, dpi = 300, ...) {
  d <- attr(x, "diagnostics")
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  if (!has_ggplot) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::hist(d$samples, freq = FALSE, breaks = 50,
                   col = "grey80", border = "white",
                   main = "Prior overlay: empirical vs mixture",
                   xlab = sprintf("Arm parameter (%s kernel)", x$family))
    graphics::lines(d$overlay$x, d$overlay$density_mixture,
                    col = "tomato", lwd = 2)
    graphics::legend("topright",
                     legend = c("Empirical samples", "Fitted mixture"),
                     fill = c("grey80", NA), border = c("white", NA),
                     col = c(NA, "tomato"), lty = c(NA, 1), lwd = c(NA, 2),
                     bty = "n")
    graphics::abline(v = d$tail_table$threshold, col = "grey60", lty = 2)
    return(invisible(NULL))
  }
  ## ggplot2 path: density overlay with tail-threshold rug.
  df_emp <- data.frame(x = d$samples)
  df_mix <- d$overlay
  df_tail <- d$tail_table
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data = df_emp, ggplot2::aes(x = .data$x),
                          fill = "grey80", colour = NA, alpha = 0.8) +
    ggplot2::geom_line(data = df_mix,
                       ggplot2::aes(x = .data$x, y = .data$density_mixture),
                       colour = "#C44E52", linewidth = 0.8) +
    ggplot2::geom_rug(data = df_tail,
                      ggplot2::aes(x = .data$threshold),
                      sides = "b", length = ggplot2::unit(0.03, "npc"),
                      colour = "#4C72B0") +
    ggplot2::labs(x = sprintf("Arm parameter (%s kernel)", x$family),
                  y = "Density",
                  title = sprintf("Prior overlay (max |tail error| = %.4f)",
                                  d$max_tail_error)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (!is.null(file)) {
    ggplot2::ggsave(filename = file, plot = p,
                    width = width, height = height, dpi = dpi, ...)
  }
  p
}

## Beta-mixture EM fallback (used only when RBesT is unavailable).
##
## Model selection has two stages. First, the elbow rule on the forward-KL
## decrement (as before) picks a candidate L: the smallest L past which
## adding a component no longer improves KL by at least tol_kl. But forward
## KL is mass-covering, not tail-accurate (see fit_mixture() docs) -- the
## elbow pick can still miss tol_tail, exactly when a component that barely
## moves the *bulk* of the density materially corrects a *tail*. So second,
## if the elbow pick's tail error exceeds tol_tail (and tol_tail/
## tail_quantiles were supplied), search does not stop there: it continues
## past the elbow, one L at a time up to n_components_max, and returns the
## smallest L that satisfies tol_tail. If none does, the fit with the
## smallest maximum tail error among those searched is returned (adding a
## component is a fresh EM fit that may settle at a different local optimum,
## so more components do not monotonically reduce the tail error -- returning
## fits[[n_components_max]] blindly could hand back a fit strictly worse than
## one already computed).
.em_beta_mixture <- function(x, n_components_max = 5L, tol_kl = 1e-3,
                             tail_quantiles = NULL, tol_tail = NULL) {
  fits <- vector("list", n_components_max)
  prev_kl <- Inf
  elbow_L <- NA_integer_
  for (L in seq_len(n_components_max)) {
    fits[[L]] <- .em_beta(x, L = L)
    if (is.na(elbow_L) && L > 1L &&
        (prev_kl - fits[[L]]$kl) < tol_kl) {
      elbow_L <- L - 1L
    }
    prev_kl <- fits[[L]]$kl
  }
  if (is.na(elbow_L)) elbow_L <- n_components_max

  if (!is.null(tol_tail) && !is.null(tail_quantiles)) {
    best_L <- elbow_L
    best_err <- Inf
    for (L in elbow_L:n_components_max) {
      tt <- .tail_error_table(x, "beta", fits[[L]]$weights,
                              fits[[L]]$components, tail_quantiles)
      err <- max(abs(tt$error))
      if (err <= tol_tail) return(fits[[L]])
      if (err < best_err) {
        best_err <- err
        best_L <- L
      }
    }
    return(fits[[best_L]])
  }
  fits[[elbow_L]]
}

.em_beta <- function(x, L, max_iter = 200L, tol = 1e-6) {
  n <- length(x)
  if (L == 1L) {
    ## Method-of-moments fit.
    mu <- mean(x); v <- stats::var(x)
    common <- mu * (1 - mu) / max(v, 1e-8) - 1
    a <- mu * common; b <- (1 - mu) * common
    a <- max(a, 0.5); b <- max(b, 0.5)
    components <- list(list(a = a, b = b))
    weights <- 1
  } else {
    ## Initialise via deterministic quantile-based partition. This avoids
    ## consuming the caller's RNG state and gives a reproducible warm start
    ## without overriding the user-supplied seed.
    breaks <- stats::quantile(x, probs = seq(0, 1, length.out = L + 1L),
                              names = FALSE)
    breaks[1L]      <- -Inf
    breaks[L + 1L]  <-  Inf
    z <- as.integer(cut(x, breaks = breaks, include.lowest = TRUE))
    weights <- as.numeric(table(factor(z, levels = seq_len(L)))) / n
    components <- lapply(seq_len(L), function(l) {
      sub <- x[z == l]
      mu <- if (length(sub) > 1) mean(sub) else mean(x)
      v  <- if (length(sub) > 1) stats::var(sub) else stats::var(x)
      common <- mu * (1 - mu) / max(v, 1e-8) - 1
      list(a = max(mu * common, 0.5), b = max((1 - mu) * common, 0.5))
    })
    for (it in seq_len(max_iter)) {
      ## E-step. vapply over components keeps log_dens an n x L numeric
      ## matrix regardless of L (sapply would collapse to a vector when L==1).
      log_dens <- vapply(seq_len(L), function(l)
        log(weights[l]) +
          stats::dbeta(x, components[[l]]$a, components[[l]]$b, log = TRUE),
        numeric(n))
      m <- apply(log_dens, 1, max)
      lse <- m + log(rowSums(exp(log_dens - m)))
      gamma <- exp(log_dens - lse)
      ## M-step.
      Nl <- colSums(gamma)
      weights_new <- Nl / n
      components_new <- lapply(seq_len(L), function(l) {
        if (Nl[l] < 1e-8) return(components[[l]])
        w <- gamma[, l] / Nl[l]
        mu <- sum(w * x); v <- sum(w * (x - mu)^2)
        common <- mu * (1 - mu) / max(v, 1e-8) - 1
        list(a = max(mu * common, 0.5), b = max((1 - mu) * common, 0.5))
      })
      comp_deltas <- vapply(seq_len(L), function(l)
        max(abs(unlist(components_new[[l]]) - unlist(components[[l]]))),
        numeric(1))
      delta <- max(abs(weights_new - weights), comp_deltas)
      weights <- weights_new; components <- components_new
      if (delta < tol) break
    }
  }
  ## Compute approximate empirical KL via a high-resolution density grid.
  grid <- seq(1e-4, 1 - 1e-4, length.out = 1000)
  approx_pdf <- vapply(grid, function(z)
    sum(weights * vapply(components, function(c)
      stats::dbeta(z, c$a, c$b), numeric(1))), numeric(1))
  kde <- stats::density(x, from = 0, to = 1, n = 1000)
  pdf_emp <- stats::approx(x = kde$x, y = kde$y, xout = grid, rule = 2)$y
  pdf_emp <- pmax(pdf_emp, 1e-12)
  kl <- mean(pdf_emp * (log(pdf_emp) - log(pmax(approx_pdf, 1e-12)))) *
    (max(grid) - min(grid))
  list(weights = weights, components = components, kl = kl)
}
