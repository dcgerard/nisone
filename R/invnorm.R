#' The inverse normal distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the inverse normal distribution when parameterized by the mean and
#' standard deviation of the inverse (reciprocal).
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n sample size
#' @param imean vector of means of inverse.
#' @param isd vector of standard deviations of inverse.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P(X<=x)
#'     otherwise, P(X>x).
#'
#' @returns Either a random sample (\code{rinvnorm}),
#'     the density (\code{dinvnorm}), the tail
#'     probability (\code{pinvnorm}), or the quantile
#'     (\code{qinvnorm}) of the inverse normal distribution.
#'
#' @examples
#' x <- seq(-4, 4, length.out = 300)
#' y <- dinvnorm(x = x, imean = 0.1, isd = 1)
#' graphics::plot(x, y, type = "l")
#'
#' p <- pinvnorm(q = x, imean = 0.1, isd = 1)
#' graphics::plot(x, p, type = "l", ylim = c(0, 1))
#' graphics::abline(h = c(0, 1), lty = 2)
#'
#' qinvnorm(p = c(0.025, 0.5, 0.975), imean = 0.1, isd = 1)
#'
#' s <- rinvnorm(n = 10000, imean = 0.1, isd = 1)
#' stats::quantile(s, c(0.025, 0.5, 0.975))
#'
#' @name invnorm
#'
#' @references
#' \itemize{
#'  \item{Robert, C. (1991). Generalized inverse normal distributions. Statistics & Probability Letters, 11(1), 37-41. \doi{10.1016/0167-7152(91)90174-P}}
#' }
#'
#' @author David Gerard
NULL

#' @describeIn invnorm Density function.
#'
#' @export
dinvnorm <- function(x, imean = 0, isd = 1, log = FALSE) {
  dval <- stats::dnorm(x = 1 / x, mean = imean, sd = isd, log = TRUE) - 2 * log(abs(x))
  if (!log) {
    dval <- exp(dval)
  }
  return(dval)
}

#' @describeIn invnorm Probability function.
#'
#' @export
pinvnorm <- function(q, imean = 0, isd = 1, lower.tail = TRUE, log.p = FALSE) {
  pvec <- stats::pnorm((imean * q - 1) / (isd * q)) - stats::pnorm(imean / isd)
  pvec[q > 0] <- pvec[q > 0] + 1
  pvec[q == 0] <- 1 - stats::pnorm(imean / isd)

  if (!lower.tail) {
    pvec <- 1 - pvec
  }
  if (log.p) {
    pvec <- log(pvec)
  }

  return(pvec)
}

#' @describeIn invnorm Quantile function.
#'
#' @export
qinvnorm <- function(p, imean = 0, isd = 1) {
  cutval <- stats::pnorm(imean / isd)
  pvec <- rep(NA_real_, length.out = length(p))
  pvec[p < 1 - cutval] <- 1 / (imean - isd * stats::qnorm(p[p < 1 - cutval] + cutval))
  pvec[p >= 1 - cutval] <- 1 / (imean - isd * stats::qnorm(p[p >= 1 - cutval] + cutval - 1))
  return(pvec)
}

#' @describeIn invnorm Random generation.
#'
#' @export
rinvnorm <- function(n, imean = 0, isd = 1) {
  1 / stats::rnorm(n = n, mean = imean, sd = isd)
}
