#' Bayesian credible interval when n=1
#'
#' @param x a double
#' @param level The level of the credible interval
#'
#' @author David Gerard
#'
#' @export
bci1 <- function(x, level = 0.95) {
  stopifnot(level >= 0, level <= 1)
  stopifnot(length(x) == 1, length(level) == 1)
  alpha <- 1 - level
  retlist <- list()
  retlist$lower <- qinvnorm(p = alpha / 2, imean = 1 / x, isd = 1 / abs(x))
  retlist$upper <- qinvnorm(p = 1 - alpha / 2, imean = 1 / x, isd = 1 / abs(x))
  retlist$med <- qinvnorm(p = 0.5, imean = 1 / x, isd = 1 / abs(x))
  return(retlist)
}

#' True level of the Bayesian credible interval at a given
#' coefficient of variation
#'
#' @param level The credible interval level.
#' @param lambda One over the coefficient of variation. mu/sigma
#'
#' @return The true level of the (1-alpha) credible interval.
#'
#' @author David Gerard
#'
#' @export
blevel <- function(level, lambda) {
  alpha <- 1 - level
  c1 <- qinvnorm(p = alpha / 2, imean = 1, isd = 1)
  c2 <- qinvnorm(p = 1 - alpha / 2, imean = 1, isd = 1)

  if (c1 <= 0 && c2 >= 0) {
    tl <- stats::pnorm(lambda * (1 - 1 / c2)) + 1 - stats::pnorm(lambda * (1 - 1 / c1))
  } else if (c1 > 0 && c2 > 0) {
    tl <- stats::pnorm(lambda * (1 - 1 / c2)) - stats::pnorm(lambda * (1 - 1 / c1))
  } else {
    stop("invalid c1 or c2")
  }
  return(tl)
}
