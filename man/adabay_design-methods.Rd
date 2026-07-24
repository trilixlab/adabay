\name{adabay_design-methods}
\alias{adabay_design-methods}
\alias{print.adabay_design}
\title{Print method for objects of class \code{adabay_design}}
\usage{
\method{print}{adabay_design}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_design} object from \code{\link{set_design}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the endpoint, look schedule, effect scale, alternative,
null-effect value and the per-arm cumulative sample-size (or exposure, or
event-count) vectors of a Bayesian group sequential design skeleton.
}
\seealso{\code{\link{set_design}}.}
