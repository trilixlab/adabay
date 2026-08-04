\name{adabay_calibration_best-methods}
\alias{adabay_calibration_best-methods}
\alias{print.adabay_calibration_best}
\title{Print method for objects of class \code{adabay_calibration_best}}
\usage{
\method{print}{adabay_calibration_best}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_calibration_best} object, the \code{$best} element
returned by \code{\link{calibrate_design}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the selected threshold pair \code{(p, q)} together with the
type I error rate under \eqn{H_0}, the power under \eqn{H_1} and the
expected sample size under \eqn{H_1} at those thresholds. The futility
threshold is reported as \code{"disabled"} when \code{q >= 1}. A note
records that the reported type I error rate and power follow the
calibration's futility convention: under the default non-binding rule the
efficacy-crossing probability ignores futility crossings, so the full
per-hypothesis operating characteristics in \code{$oc_h0} and \code{$oc_h1}
should be consulted for the realised behaviour.
}
\details{
The object is the single grid point that meets the \code{alpha_target} and
\code{power_target} constraints at minimum expected sample size under the
alternative, selected from the user-supplied \code{efficacy_grid} by
\code{futility_grid} sweep. It is deliberately not of class
\code{adabay_oc}: its \code{power} field is the efficacy-crossing
probability under \eqn{H_1}, whereas \code{oc_h1$alpha} carries the same
quantity under the \code{adabay_oc} naming, and printing it through the
\code{adabay_oc} method would mislabel it.
}
\seealso{\code{\link{calibrate_design}}, \code{\link{build_cache}},
\code{\link{evaluate_design}}.}
