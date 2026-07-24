\name{adabay_accrual-methods}
\alias{adabay_accrual-methods}
\alias{print.adabay_accrual}
\title{Print method for objects of class \code{adabay_accrual}}
\usage{
\method{print}{adabay_accrual}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_accrual} object from \code{\link{set_accrual}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the recruitment model (Poisson, piecewise, or callback),
its rate (where applicable), the per-patient follow-up and the dropout
specification.
}
\seealso{\code{\link{set_accrual}}.}
