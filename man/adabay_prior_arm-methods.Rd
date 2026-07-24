\name{adabay_prior_arm-methods}
\alias{adabay_prior_arm-methods}
\alias{print.adabay_prior_arm}
\alias{plot.adabay_prior_arm}
\title{Methods for objects of class \code{adabay_prior_arm}}
\usage{
\method{print}{adabay_prior_arm}(x, \dots)

\method{plot}{adabay_prior_arm}(x, file = NULL, width = 7, height = 4, dpi = 300, \dots)
}
\arguments{
\item{x}{A \code{adabay_prior_arm} object returned by \code{\link{fit_mixture}}.}
\item{file}{Optional output path. When supplied, the plot is written via
  \code{ggplot2::ggsave}; the file extension determines the device.}
\item{width, height}{Plot dimensions in inches when \code{file} is supplied.}
\item{dpi}{Resolution for raster outputs.}
\item{\dots}{Additional arguments forwarded to \code{ggsave} (for
  \code{plot}) or ignored (for \code{print}).}
}
\value{
\code{print} returns its input invisibly. \code{plot} returns a
\code{ggplot} object (or \code{NULL}, invisibly, when \code{ggplot2} is
not installed).
}
\description{
\code{print()} reports the kernel family, the number of mixture
components with each component's weight and conjugate parameters (in the
same format as \code{print.adabay_prior}), the fit method, the empirical
forward Kullback--Leibler divergence, the seed passed to
\code{\link{fit_mixture}} (if any), and the per-quantile tail-probability
error table. \code{plot()} overlays the empirical-sample density on the
fitted mixture density and marks the tail-quantile thresholds as a rug.
}
\seealso{\code{\link{fit_mixture}}, \code{\link{set_prior}}.}
