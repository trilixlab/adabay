\name{adabay_calibration-methods}
\alias{adabay_calibration-methods}
\alias{print.adabay_calibration}
\title{Print method for objects of class \code{adabay_calibration}}
\usage{
\method{print}{adabay_calibration}(x, \dots)
}
\arguments{
\item{x}{A \code{adabay_calibration} object from
\code{\link{calibrate_design}}.}
\item{\dots}{Currently unused.}
}
\value{The input, invisibly.}
\description{
Pretty-print the size of the evaluated threshold grid, followed by the
selected design in \code{$best} rendered through
\code{print.adabay_calibration_best}. When no grid point met the
\code{alpha_target} and \code{power_target} constraints, \code{$best} is
\code{NULL} and this is reported instead. The grid itself is not printed:
it is a data frame of one row per threshold combination and is available
as \code{$grid}.
}
\details{
The object returned by \code{\link{calibrate_design}} is a list of two
elements: \code{$best}, of class \code{adabay_calibration_best}, and
\code{$grid}, a data frame with columns \code{p}, \code{q},
\code{type_I}, \code{power} and \code{E_N}. Printing the object without a
method would emit the whole grid, which at a typical sweep is several
hundred rows.
}
\seealso{\code{\link{calibrate_design}}, \code{\link{build_cache}},
\code{\link{evaluate_design}}.}
