#' The generalized inverse normal distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the \emph{generalized} inverse normal distribution when parameterized
#' by the shape, mean, and standard deviation of the inverse (reciprocal).
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param mu vector of means of inverse.
#' @param tau vector of standard deviations of inverse.
#' @param alpha vector of shape parameters
#' @param log logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P(X<=x)
#'     otherwise, P(X>x).
#'
#' @references
#' \itemize{
#'  \item{Robert, C. (1991). Generalized inverse normal distributions. Statistics & Probability Letters, 11(1), 37-41. \doi{10.1016/0167-7152(91)90174-P}}
#' }
#'
#' @return Either a random sample (\code{rginvnorm}),
#'     the density (\code{dginvnorm}), the tail
#'     probability (\code{pginvnorm}), or the quantile
#'     (\code{qginvnorm}) of the inverse normal distribution.
#'
#' @name ginvnorm
#'
#' @author David Gerard
NULL

#' @describeIn ginvnorm Density function.
#'
#' @export
dginvnorm <- function(x, alpha, mu = 0, tau = 1, log = FALSE) {
  stopifnot(alpha > 1, tau > 0)

  ## normalizing constant
  lK <- (alpha - 1) * log(tau) +
    -mu^2 / (2 * tau^2) +
    (alpha - 1) * log(2) / 2 +
    lgamma((alpha - 1) / 2) +
    hg1f1_special(a = alpha - 1, z = mu^2 / (2 * tau^2), log = TRUE) # same as log(chgm(a = (alpha - 1)/2, b = 0.5, x = mu^2 / (2 * tau^2)))

  dval <- -1/(2 * tau^2) * (1 / x - mu)^2 - alpha * log(abs(x)) - lK

  if (!log) {
    dval <- exp(dval)
  }
  return(dval)
}

#' @describeIn ginvnorm Probability function.
#'
#' @export
pginvnorm <- function(q, alpha, mu = 0, tau = 1, lower.tail = TRUE) {
  ## some infinite values make integrate() complain, so give a slight nudge
  nudge <- stats::runif(n = 1, min = -1, max = 1) / 10^6
  f <- function(x) dginvnorm(x = x, alpha = alpha, mu = mu, tau = tau, log = FALSE)
  if (lower.tail) {
    ret <- stats::integrate(f = f, lower = -Inf, upper = q + nudge)$value
  } else {
    ret <- stats::integrate(f = f, lower = q + nudge, upper = Inf)$value
  }
  return(ret)
}

#' Find bounds of quantile in qginvnorm, starting at mode
#'
#' @noRd
find_bounds_q <- function(p, alpha, mu = 0, tau = 1) {
  mvec <- xginvnorm(alpha = alpha, mu = mu, tau = tau)
  mcent <- ifelse(mu < 0, mvec[[1]], mvec[[2]])
}

#' @describeIn ginvnorm Quantile function.
#'
#' @export
qginvnorm <- function(p, alpha, mu = 0, tau = 1) {
  f <- function(q) pginvnorm(q = q, alpha = alpha, mu = mu, tau = tau) - p
  stats::uniroot(f = f, interval = c(-100, 100))$root
}

#' @describeIn ginvnorm Mean.
#'
#' @export
mginvnorm <- function(alpha, mu = 0, tau = 1) {
  mu / tau^2 * chgm((alpha - 1) / 2, 3 / 2, mu^2 / (2 * tau^2)) / chgm((alpha - 1) / 2, 1 / 2, mu^2 / (2 * tau^2))
}

#' @describeIn ginvnorm Modes.
#'
#' @export
xginvnorm <- function(alpha, mu = 0, tau = 1) {
  c(
    -(mu + sqrt(mu^2 + 4 * alpha * tau^2)) / (2 * alpha * tau^2),
    -(mu - sqrt(mu^2 + 4 * alpha * tau^2)) / (2 * alpha * tau^2)
  )
}
