#' The generalized inverse normal distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the \emph{generalized} inverse normal distribution when parameterized
#' by the shape, mean, and standard deviation of the inverse (reciprocal).
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n sample size
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
#'
#' @examples
#' set.seed(50)
#' samp <- rginvnorm(n = 1000, alpha = 4, mu = 0.5, tau = 1)
#' x <- seq(min(samp), max(samp), length.out = 500)
#' y <- dginvnorm(x = x, alpha = 4, mu = 0.5, tau = 1)
#' modes <- xginvnorm(alpha = 4, mu = 0.5, tau = 1)
#' graphics::hist(samp, freq = FALSE, breaks = 100, xlab = "x")
#' graphics::lines(x, y, col = "red")
#' graphics::abline(v = modes, col = "blue", lty = 2)
#'
NULL

#' @describeIn ginvnorm Density function.
#'
#' @export
dginvnorm <- function(x, alpha, mu = 0, tau = 1, log = FALSE) {
  stopifnot(alpha > 1, tau > 0)

  lK <- .gin_log_K(alpha = alpha, mu = mu, tau = tau)

  dval <- lK - (1 / x - mu)^2 / (2 * tau^2) - alpha * log(abs(x))

  if (!log) {
    dval <- exp(dval)
  }
  return(dval)
}

#' @describeIn ginvnorm Probability function.
#'
#' @export
pginvnorm <- function(q, alpha, mu = 0, tau = 1, lower.tail = TRUE) {
  stopifnot(alpha > 1, tau > 0)
  # Calculate P(X >= 1 / q) where X = 1/Z and Z is ginvnorm

  ## Precompute since can be slow for large mu/tau
  lK <- .gin_log_K(alpha = alpha, mu = mu, tau = tau)

  # density of x = 1/z
  f <- function(x) {
    exp(lK - (x - mu)^2 / (2 * tau^2) + (alpha - 2) * log(abs(x)))
  }

  ## Precompute this one because it is used whenever q >= 0
  if (any(q >= 0)) {
    negval <- stats::integrate(f = f, lower = -Inf, upper = 0)$value
  }

  lpvec <- rep(NA_real_, length.out = length(q))
  for (i in seq_along(q)) {
    if (is.na(q[[i]])) {
      lpvec[[i]] <- NA_real_
    } else if (q[[i]] == -Inf) {
      lpvec[[i]] <- 0
    } else if (q[[i]] == Inf) {
      lpvec[[i]] <- 1
    } else if (q[[i]] == 0) {
      # If q = 0, then integral is from -Inf to 0
      lpvec[[i]] <- negval
    } else if (q[[i]] > 0) {
      # If q > 0, then integral is from -Inf to 0 and from 1/q to Inf
      lpvec[[i]] <- negval + stats::integrate(f = f, lower = 1 / q[[i]], upper = Inf)$value
    } else {
      # If q < 0, then integral is from 1/q to 0
      lpvec[[i]] <- stats::integrate(f = f, lower = 1 / q[[i]], upper = 0)$value
    }
  }

  if (isTRUE(lower.tail)) {
    ret <- lpvec
  } else {
    ret <- 1 - lpvec
  }

  return(ret)
}

#' @describeIn ginvnorm Quantile function.
#'
#' @export
qginvnorm <- function(p, alpha, mu = 0, tau = 1) {
  stopifnot(alpha > 1, tau > 0)

  qvec <- rep(NA_real_, length.out = length(p))
  qvec[p == 0] <- -Inf
  qvec[p == 1] <- Inf

  if (all(p == 0 | p == 1)) {
    return(qvec)
  }

  ## Bracket
  minp <- min(p[p > 0 & p < 1])
  maxp <- max(p[p > 0 & p < 1])

  qgrid <- c(-exp(1:10), 0, exp(1: 10))
  pgrid <- pginvnorm(q = qgrid, alpha = alpha, mu = mu, tau = tau)

  while(minp < pgrid[[1]]) {
    qgrid <- c(qgrid[[1]] * exp(1), qgrid)
    pgrid <- c(pginvnorm(q = qgrid[[1]], alpha = alpha, mu = mu, tau = tau), pgrid)
  }
  while(maxp > pgrid[[length(pgrid)]]) {
    qgrid <- c(qgrid, qgrid[[length(qgrid)]] * exp(1))
    pgrid <- c(pgrid, pginvnorm(q = qgrid[[length(qgrid)]], alpha = alpha, mu = mu, tau = tau))
  }

  ## Get quantiles
  for (i in seq_along(p)) {
    if (p[[i]] > 0 && p[[i]] < 1) {
      lower <- max(qgrid[pgrid < p[[i]]])
      upper <- min(qgrid[pgrid > p[[i]]])
      f <- function(q) pginvnorm(q = q, alpha = alpha, mu = mu, tau = tau) - p[[i]]
      qvec[[i]] <- stats::uniroot(f = f, interval = c(lower, upper))$root
    }
  }

  return(qvec)
}

#' @describeIn ginvnorm Random generation.
#'
#' @export
rginvnorm <- function(n, alpha, mu = 0, tau = 1) {
  qginvnorm(p = stats::runif(n = n), alpha = alpha, mu = mu, tau = tau)
}

#' @describeIn ginvnorm Mean.
#'
#' @export
mginvnorm <- function(alpha, mu = 0, tau = 1) {
  # mu / tau^2 * chgm((alpha - 1) / 2, 3 / 2, mu^2 / (2 * tau^2)) / chgm((alpha - 1) / 2, 1 / 2, mu^2 / (2 * tau^2))
  mu / tau^2 * .hyp1f1((alpha - 1) / 2, 3 / 2, mu^2 / (2 * tau^2)) / .hyp1f1((alpha - 1) / 2, 1 / 2, mu^2 / (2 * tau^2))
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
