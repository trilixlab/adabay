#' Set accrual: recruitment, follow-up and dropout
#'
#' Configures the recruitment process, the per-patient follow-up duration and
#' the dropout process used to compute the expected study duration alongside
#' the operating characteristics.
#'
#' @param model Recruitment model. One of:
#'   \describe{
#'     \item{\code{"poisson"}}{homogeneous Poisson process with constant rate
#'           \code{rate}.}
#'     \item{\code{"piecewise"}}{piecewise-constant Poisson process with
#'           breakpoints in \code{breakpoints} and rates in \code{rate}.}
#'     \item{\code{"callback"}}{user-supplied callback that returns
#'           recruitment times for a given target sample size.}
#'   }
#' @param rate Recruitment rate, pooled across both arms (matching
#'   \code{n_per_look}/\code{exposure_per_look}/\code{d_per_look} in
#'   \code{\link{set_design}}, which are likewise pooled). Split into
#'   per-arm rates by the design's \code{allocation_ratio}:
#'   \code{rate_c = rate / (1 + allocation_ratio)}, \code{rate_t = rate_c *
#'   allocation_ratio}, so \code{rate_c + rate_t == rate} exactly. With the
#'   default \code{allocation_ratio = 1} this is a 1:1 split.
#' @param breakpoints Numeric vector of calendar-time breakpoints (only used
#'   with \code{model = "piecewise"}).
#' @param callback Function with signature \code{function(n, ...)} returning
#'   recruitment times (only used with \code{model = "callback"}).
#' @param follow_up Per-patient follow-up duration. Either a single
#'   non-negative numeric or a function with signature
#'   \code{function(n, ...)} returning a vector of follow-up times.
#' @param dropout Dropout model. One of \code{"none"} (the default,
#'   administrative censoring at the look times) or \code{"exponential"}
#'   with rate \code{dropout_rate}.
#' @param dropout_rate Per-arm dropout rate when \code{dropout =
#'   "exponential"}. Length-1 (pooled) or length-2 numeric.
#' @return An object of class \code{"adabay_accrual"}.
#' @examples
#' set_accrual(model = "poisson", rate = 40)
#' @export
set_accrual <- function(model = c("poisson", "piecewise", "callback"),
                        rate = NULL,
                        breakpoints = NULL,
                        callback = NULL,
                        follow_up = 0,
                        dropout = c("none", "exponential"),
                        dropout_rate = NULL) {
  model   <- .match_choice(model, c("poisson", "piecewise", "callback"))
  dropout <- .match_choice(dropout, c("none", "exponential"))
  if (model == "poisson") {
    if (is.null(rate))
      stop("Provide 'rate' for the Poisson recruitment model.",
           call. = FALSE)
    .assert_numeric(rate, len = 1L, positive = TRUE)
  } else if (model == "piecewise") {
    if (is.null(rate) || is.null(breakpoints))
      stop("Provide 'rate' and 'breakpoints' for the piecewise recruitment model.",
           call. = FALSE)
    .assert_numeric(breakpoints, positive = TRUE)
    if (length(rate) != length(breakpoints) + 1L)
      stop("'rate' must have length(breakpoints) + 1 entries.",
           call. = FALSE)
    .assert_numeric(rate, positive = TRUE)
  } else if (model == "callback") {
    if (!is.function(callback))
      stop("Provide a function 'callback' for the callback recruitment model.",
           call. = FALSE)
  }
  if (is.numeric(follow_up)) {
    .assert_numeric(follow_up, len = 1L)
    if (follow_up < 0) stop("'follow_up' must be non-negative.", call. = FALSE)
  } else if (!is.function(follow_up)) {
    stop("'follow_up' must be a non-negative numeric scalar or a function.",
         call. = FALSE)
  }
  if (dropout == "exponential") {
    if (is.null(dropout_rate))
      stop("Provide 'dropout_rate' for the exponential dropout model.",
           call. = FALSE)
    .assert_numeric(dropout_rate, positive = TRUE)
    if (!length(dropout_rate) %in% c(1L, 2L))
      stop("'dropout_rate' must have length 1 or 2.", call. = FALSE)
  }

  structure(
    list(model        = model,
         rate         = rate,
         breakpoints  = breakpoints,
         callback     = callback,
         follow_up    = follow_up,
         dropout      = dropout,
         dropout_rate = dropout_rate),
    class = "adabay_accrual")
}
