\name{adabay_cache-methods}
\alias{adabay_cache-methods}
\alias{print.adabay_cache}
\title{Print method for objects of class \code{adabay_cache}}
\usage{
\method{print}{adabay_cache}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_cache} object from \code{\link{build_cache}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the endpoint, look count, virtual-trial count, and the
sorted vector of effect-size thresholds cached for grid evaluation by
\code{\link{evaluate_design}} or \code{\link{calibrate_design}}.
}
\seealso{\code{\link{build_cache}}, \code{\link{evaluate_design}},
\code{\link{calibrate_design}}.}
