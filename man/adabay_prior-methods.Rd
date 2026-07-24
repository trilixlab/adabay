\name{adabay_prior-methods}
\alias{adabay_prior-methods}
\alias{print.adabay_prior}
\title{Print method for objects of class \code{adabay_prior}}
\usage{
\method{print}{adabay_prior}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_prior} object from \code{\link{set_prior}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the endpoint and the per-arm conjugate (possibly mixture)
prior specification: family, component count, and per-component weights
and hyperparameters.
}
\seealso{\code{\link{set_prior}}, \code{\link{fit_mixture}}.}
