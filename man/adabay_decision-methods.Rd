\name{adabay_decision-methods}
\alias{adabay_decision-methods}
\alias{print.adabay_decision}
\title{Print method for objects of class \code{adabay_decision}}
\usage{
\method{print}{adabay_decision}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_decision} object from \code{\link{set_decision}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the efficacy and futility rules, the binding-futility flag,
and indicate whether the rules are common across looks or per-look.
}
\seealso{\code{\link{set_decision}}.}
