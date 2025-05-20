#' Bayesian credible interval when n=1
#'
#' @param x a double
#' @param A The prior center.
#' @param level The level of the credible interval
#'
#' @author David Gerard
#'
#' @export
bci1 <- function(x, A = 0, level = 0.95) {
  stopifnot(level >= 0, level <= 1)
  stopifnot(length(level) == 1, length(A) == 1)
  retmat <- matrix(NA_real_, ncol = 4, nrow = length(x))
  colnames(retmat) <- c("x", "med", "lower", "upper")
  retmat[, "x"] <- x
  x <- x - A
  alpha <- 1 - level
  for (i in seq_along(x)) {
    retmat[i, "lower"] <- qinvnorm(p = alpha / 2, imean = 1 / x[[i]], isd = 1 / abs(x[[i]])) + A
    retmat[i, "upper"] <- qinvnorm(p = 1 - alpha / 2, imean = 1 / x[[i]], isd = 1 / abs(x[[i]])) + A
    retmat[i, "med"] <- qinvnorm(p = 0.5, imean = 1 / x[[i]], isd = 1 / abs(x[[i]])) + A
  }
  return(retmat)
}

#' True level of the Bayesian credible interval at a given
#' coefficient of variation
#'
#' @param level The credible interval level.
#' @param nu (mu - A)/sigma
#'
#' @return The true level of the (1-alpha) credible interval.
#'
#' @author David Gerard
#'
#' @export
blevel <- function(level, nu) {
  alpha <- 1 - level
  c1 <- qinvnorm(p = alpha / 2, imean = 1, isd = 1)
  c2 <- qinvnorm(p = 1 - alpha / 2, imean = 1, isd = 1)

  if (c1 <= 0 && c2 >= 0) {
    tl <- stats::pnorm(nu * (1 - 1 / c2)) + 1 - stats::pnorm(nu * (1 - 1 / c1))
  } else if (c1 > 0 && c2 > 0) {
    tl <- stats::pnorm(nu * (1 - 1 / c2)) - stats::pnorm(nu * (1 - 1 / c1))
  } else {
    stop("invalid c1 or c2")
  }
  return(tl)
}

bcin <- function(x, A = 0, level = 0.95) {
  y <- x - A
  alpha <- length(x) + 1
  mu <- sum(y) / sum(y^2)
  tau <- 1 / sqrt(sum(y^2))

  al <- (1 - level) / 2
  mginvnorm(alpha = alpha, mu = mu, tau = tau)
  qginvnorm(p = 0.5, alpha = alpha, mu = mu, tau = tau)
}
