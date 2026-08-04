// adabay: C++ implementations of per-look posterior tail probabilities.
//
// Each function takes a single threshold e and vectors of posterior
// hyperparameters across R virtual trials, and returns a length-R numeric
// vector of P(Delta > e | data). Integrals are evaluated by Gauss-Legendre
// quadrature with a fixed node count, which avoids the per-trial integrate()
// call and makes the inner loop trivially vectorised in C++.
//
// The number of quadrature nodes (NODES) is chosen large enough that the
// numerical error is well below the Monte Carlo error at typical R = 1e5.
//
// Compiled by R CMD INSTALL and exposed through the Rcpp::compileAttributes()
// stubs in src/RcppExports.cpp, registered via useDynLib(adabay, .registration
// = TRUE). The R wrappers in posterior.R call these entry points
// unconditionally -- there is no R-level fallback for the quadrature scales.

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

static const int NODES = 64;

// Pre-computed Gauss-Legendre nodes and weights on [-1, 1] for N = 64.
// These are loaded lazily on first use; we use the simple iterative
// recurrence to compute them rather than embedding 64 floating-point
// constants inline.
static std::vector<double> gl_nodes;
static std::vector<double> gl_weights;

static void compute_gl(int n) {
  // Standard Newton iteration on the Legendre polynomial.
  gl_nodes.assign(n, 0.0);
  gl_weights.assign(n, 0.0);
  const double pi = 3.14159265358979323846;
  for (int i = 0; i < (n + 1) / 2; ++i) {
    double z = std::cos(pi * (i + 0.75) / (n + 0.5));
    double z1, pp;
    do {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 1; j <= n; ++j) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (std::fabs(z - z1) > 1e-14);
    gl_nodes[i]         = -z;
    gl_nodes[n - 1 - i] =  z;
    double w = 2.0 / ((1.0 - z * z) * pp * pp);
    gl_weights[i]         = w;
    gl_weights[n - 1 - i] = w;
  }
}

static void ensure_gl() {
  if ((int)gl_nodes.size() != NODES) compute_gl(NODES);
}

// Internal helpers that map [-1, 1] Gauss-Legendre to [a, b].
static inline double gl_x(int k, double a, double b) {
  return 0.5 * (b - a) * gl_nodes[k] + 0.5 * (b + a);
}
static inline double gl_w(int k, double a, double b) {
  return 0.5 * (b - a) * gl_weights[k];
}

// [[Rcpp::export]]
NumericVector pdelta_binary_rd_cpp(double e,
                                   NumericVector a_c, NumericVector b_c,
                                   NumericVector a_t, NumericVector b_t) {
  ensure_gl();
  int R = a_c.size();
  NumericVector out(R);
  if (e <= -1.0) { std::fill(out.begin(), out.end(), 1.0); return out; }
  if (e >=  1.0) { std::fill(out.begin(), out.end(), 0.0); return out; }
  double lo = std::max(0.0, -e);
  double hi = std::min(1.0, 1.0 - e);
  if (lo >= hi) {
    double v = (e < 0.0) ? 1.0 : 0.0;
    std::fill(out.begin(), out.end(), v);
    return out;
  }
  for (int r = 0; r < R; ++r) {
    double s = 0.0;
    for (int k = 0; k < NODES; ++k) {
      double th = gl_x(k, lo, hi);
      double w  = gl_w(k, lo, hi);
      double dc = R::dbeta(th, a_c[r], b_c[r], 0);
      double pt = R::pbeta(th + e, a_t[r], b_t[r], 0, 0);
      s += w * dc * pt;
    }
    if (e < 0.0) s += R::pbeta(-e, a_c[r], b_c[r], 1, 0);
    out[r] = s;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pdelta_binary_rr_cpp(double e,
                                   NumericVector a_c, NumericVector b_c,
                                   NumericVector a_t, NumericVector b_t) {
  ensure_gl();
  int R = a_c.size();
  NumericVector out(R);
  if (e <= 0.0) { std::fill(out.begin(), out.end(), 1.0); return out; }
  double lo = 0.0;
  double hi = std::min(1.0, 1.0 / e);
  for (int r = 0; r < R; ++r) {
    double s = 0.0;
    for (int k = 0; k < NODES; ++k) {
      double th = gl_x(k, lo, hi);
      double w  = gl_w(k, lo, hi);
      double dc = R::dbeta(th, a_c[r], b_c[r], 0);
      double pt = R::pbeta(std::min(1.0, th * e), a_t[r], b_t[r], 0, 0);
      s += w * dc * pt;
    }
    out[r] = s;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pdelta_binary_or_cpp(double e,
                                   NumericVector a_c, NumericVector b_c,
                                   NumericVector a_t, NumericVector b_t) {
  ensure_gl();
  int R = a_c.size();
  NumericVector out(R);
  if (e <= 0.0) { std::fill(out.begin(), out.end(), 1.0); return out; }
  double lo = 0.0;
  double hi = 1.0 - 1e-10;
  for (int r = 0; r < R; ++r) {
    double s = 0.0;
    for (int k = 0; k < NODES; ++k) {
      double th = gl_x(k, lo, hi);
      double w  = gl_w(k, lo, hi);
      double odds_t = th / (1.0 - th) * e;
      double g = odds_t / (1.0 + odds_t);
      double dc = R::dbeta(th, a_c[r], b_c[r], 0);
      double pt = R::pbeta(g, a_t[r], b_t[r], 0, 0);
      s += w * dc * pt;
    }
    out[r] = s;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pdelta_gamma_rd_cpp(double e,
                                  NumericVector a_c, NumericVector b_c,
                                  NumericVector a_t, NumericVector b_t) {
  ensure_gl();
  int R = a_c.size();
  NumericVector out(R);
  for (int r = 0; r < R; ++r) {
    double lo = std::max(0.0, -e);
    // Upper integration limit: q-th percentile of Gamma(a_c, b_c) for q very
    // close to 1, ensuring negligible truncation error.
    double hi = R::qgamma(1.0 - 1e-12, a_c[r], 1.0 / b_c[r], 1, 0);
    if (hi <= lo) {
      out[r] = (e < 0.0) ? R::pgamma(-e, a_c[r], 1.0 / b_c[r], 1, 0) : 0.0;
      continue;
    }
    double s = 0.0;
    for (int k = 0; k < NODES; ++k) {
      double lc = gl_x(k, lo, hi);
      double w  = gl_w(k, lo, hi);
      double dc = R::dgamma(lc, a_c[r], 1.0 / b_c[r], 0);
      double pt = R::pgamma(lc + e, a_t[r], 1.0 / b_t[r], 0, 0);
      s += w * dc * pt;
    }
    if (e < 0.0) s += R::pgamma(-e, a_c[r], 1.0 / b_c[r], 1, 0);
    out[r] = s;
  }
  return out;
}

// Vectorised time-to-event (tte) sufficient-statistic simulation under staggered
// Poisson enrolment. For each of R virtual trials and each arm a in
// {control (_c), treatment (_t)}:
//   * draw enrolment (arrival) times A_{a,i} as a homogeneous Poisson
//     process with per-arm rate rate_a (cumulative Exp(rate_a) gaps);
//   * draw time-to-event T_{a,i} ~ Exp(lambda_a);
//   * form the calendar event time X_{a,i} = A_{a,i} + T_{a,i}.
// The k-th event-driven look fires at the calendar time c_k equal to the
// D_per_look[k]-th smallest X pooled across both arms. At that look the
// per-arm sufficient statistic is
//   D_{a,k} = #{ i : X_{a,i} <= c_k }                               (events)
//   E_{a,k} = sum_{ i : A_{a,i} <= c_k } min(T_{a,i}, c_k - A_{a,i}) (exposure)
// so subjects not yet enrolled by c_k contribute zero exposure and subjects
// enrolled but event-free are censored at the calendar cutoff.
//
// Returns a list with elements D_c, D_t, E_c, E_t, C, each an R x K matrix;
// C(r, k) is the calendar time c_k actually realised in trial r at look k,
// i.e. the exact (simulated, not analytic) duration of trial r if it had
// run to look k -- consumed by .expected_duration() in R/sim-engine.R to
// compute the expected trial duration directly from simulation rather than
// from a deterministic mean-field approximation.
//
// By exponential memorylessness the law of (D_a, E_a) at an event-driven
// look is invariant to the enrolment schedule (lambda_a E_a | D_a ~
// Gamma(D_a, 1)); this routine realises that law by literally drawing the
// arrival times rather than assuming simultaneous entry. The per-arm pools
// n_pool_c, n_pool_t are sized on the R side (mean enrolment by the final
// look plus an 8-SD margin) so that, with overwhelming probability, the
// last-enrolled subject in each arm arrives after the final look; the
// runtime backstop below makes that requirement explicit.
// [[Rcpp::export]]
List sim_tte_cpp(int R, int n_pool_c, int n_pool_t,
                      double lambda_c, double lambda_t,
                      double rate_c, double rate_t,
                      IntegerVector D_per_look) {
  int K = D_per_look.size();
  NumericMatrix D_c(R, K), D_t(R, K);
  NumericMatrix E_c(R, K), E_t(R, K);
  NumericMatrix C(R, K);
  std::vector<double> A_c(n_pool_c), T_c(n_pool_c), X_c(n_pool_c);
  std::vector<double> A_t(n_pool_t), T_t(n_pool_t), X_t(n_pool_t);
  std::vector<double> pooled(n_pool_c + n_pool_t);
  for (int r = 0; r < R; ++r) {
    double a = 0.0;
    for (int i = 0; i < n_pool_c; ++i) {
      a += R::rexp(1.0 / rate_c);
      A_c[i] = a;
      T_c[i] = R::rexp(1.0 / lambda_c);
      X_c[i] = a + T_c[i];
      pooled[i] = X_c[i];
    }
    a = 0.0;
    for (int i = 0; i < n_pool_t; ++i) {
      a += R::rexp(1.0 / rate_t);
      A_t[i] = a;
      T_t[i] = R::rexp(1.0 / lambda_t);
      X_t[i] = a + T_t[i];
      pooled[n_pool_c + i] = X_t[i];
    }
    std::sort(pooled.begin(), pooled.end());
    for (int k = 0; k < K; ++k) {
      int target = D_per_look[k];
      if (target > n_pool_c + n_pool_t)
        Rcpp::stop("Pool size insufficient for target event count");
      double c_k = pooled[target - 1];
      C(r, k) = c_k;
      int dc = 0, dt = 0;
      double ec = 0.0, et = 0.0;
      for (int i = 0; i < n_pool_c; ++i) {
        if (A_c[i] <= c_k) {
          ec += std::min(T_c[i], c_k - A_c[i]);
          if (X_c[i] <= c_k) ++dc;
        }
      }
      for (int i = 0; i < n_pool_t; ++i) {
        if (A_t[i] <= c_k) {
          et += std::min(T_t[i], c_k - A_t[i]);
          if (X_t[i] <= c_k) ++dt;
        }
      }
      D_c(r, k) = dc; D_t(r, k) = dt;
      E_c(r, k) = ec; E_t(r, k) = et;
    }
    // Backstop: the last-enrolled subject in each arm must arrive after the
    // final look's calendar cutoff, so every subject who could contribute an
    // event before the last look is represented in the pool. With the R-side
    // sizing this never triggers; it guards a caller that shrinks the pool.
    double c_last = pooled[D_per_look[K - 1] - 1];
    if (A_c[n_pool_c - 1] <= c_last || A_t[n_pool_t - 1] <= c_last)
      Rcpp::stop("Enrolment pool exhausted before the final look; "
                 "increase the recruitment pool size.");
  }
  return List::create(_["D_c"] = D_c, _["D_t"] = D_t,
                      _["E_c"] = E_c, _["E_t"] = E_t, _["C"] = C);
}
