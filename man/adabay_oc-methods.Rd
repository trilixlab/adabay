\name{adabay_oc-methods}
\alias{adabay_oc-methods}
\alias{print.adabay_oc}
\alias{plot.adabay_oc}
\alias{summary.adabay_oc}
\title{Methods for objects of class \code{adabay_oc}}
\usage{
\method{print}{adabay_oc}(x, \dots)

\method{plot}{adabay_oc}(x, file = NULL, width = 7, height = 4, dpi = 300, \dots)

\method{summary}{adabay_oc}(object, \dots)
}
\arguments{
\item{x, object}{A \code{adabay_oc} object returned by \code{\link{evaluate_design}}.}
\item{file}{Optional output path. When supplied, the plot is written via
  \code{ggplot2::ggsave}; the file extension (\code{.jpeg}, \code{.jpg},
  \code{.png}, \code{.pdf}) determines the device.}
\item{width, height}{Plot dimensions in inches when \code{file} is supplied.}
\item{dpi}{Resolution for raster outputs.}
\item{\dots}{Additional arguments forwarded to \code{ggsave} (for
  \code{plot}) or ignored (for \code{print} and \code{summary}).}
}
\value{
\code{print} returns its input invisibly. \code{summary} returns a named
list with elements \code{alpha} (overall type I error rate, under the
regulatory convention for non-binding futility), \code{eff_prob} (overall
efficacy probability under binding conduct), \code{fut_prob} (overall
futility probability under binding conduct), \code{ess}, \code{edur},
\code{n_trials}, \code{eff_prob_at_look} (stage efficacy stopping
probabilities) and \code{fut_prob_at_look} (stage futility stopping
probabilities). \code{plot} returns a \code{ggplot} object (or \code{NULL},
invisibly, when \code{ggplot2} is not installed).
}
\description{
Print, summary, and per-look stopping-probability bar chart for an
operating-characteristics object produced by
\code{\link{evaluate_design}}. The plot
uses \code{ggplot2::theme_bw()} and the BMJ palette via
\code{ggsci::scale_fill_bmj()} when both packages are available, falling
back to base graphics otherwise.
}
\seealso{
\code{\link{evaluate_design}}, \code{\link{build_cache}}.
}
