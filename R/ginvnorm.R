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
#' @param subdivisions The maximum number of subintervals used in \code{\link[stats]{integrate}()}.
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
#' graphics::hist(
#'   samp,
#'   freq = FALSE,
#'   breaks = 100,
#'   xlab = "x",
#'   main = "Generalized Inverse Normal Density")
#' graphics::lines(x, y, col = "#E69F00")
#' graphics::abline(v = modes, col = "#56B4E9", lty = 2)
#'
NULL

#' @describeIn ginvnorm Density function.
#'
#' @export
dginvnorm <- function(x, alpha, mu = 0, tau = 1, log = FALSE) {
  stopifnot(alpha > 1, tau > 0)
  if (alpha == 2) {
    return(dinvnorm(x = x, imean = mu, isd = tau, log = log))
  }

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
pginvnorm <- function(q, alpha, mu = 0, tau = 1, lower.tail = TRUE, subdivisions = 500L) {
  stopifnot(alpha > 1, tau > 0)
  if (alpha < 2) {
    warning("pginvnorm is not stable for alpha less than 2.")
    return(pginvnorm_old(q = q, alpha = alpha, mu = mu, tau = tau, lower.tail = lower.tail))
  } else if (alpha == 2) {
    return(pinvnorm(q = q, imean = mu, isd = tau, lower.tail = lower.tail))
  }

  # Calculate P(X >= 1 / (tau * q)) where X = 1/(tau * Z) and Z is ginvnorm
  # f(z|\alpha,\mu,\tau) = K(\alpha,\mu,\tau)|z|^{-\alpha}\exp(-(1/z-\mu)^2/(2\tau))
  # So, after change of variables
  # f(x|\alpha,\mu,\tau) = K(\alpha,\mu,\tau)|\tau|^{\alpha-1}|x|^{\alpha-2}\exp(-(x-\mu/\tau)^2/2)

  ## Precompute since can be slow for large mu/tau
  lK <- .gin_log_K(alpha = alpha, mu = mu, tau = tau)

  # density of x = 1/(tau * z)
  f <- function(x) {
    exp(lK + (alpha - 1) * log(abs(tau)) + (alpha - 2) * log(abs(x)) - (x - mu / tau)^2 / 2)
  }

  ## Get modes for distribution of x
  modes <- c(
    (mu/tau - sqrt((mu/tau)^2 + 4 * (alpha - 2))) / 2,
    (mu/tau + sqrt((mu/tau)^2 + 4 * (alpha - 2))) / 2
  )
  modes <- sort(modes)

  ## Start at modes, expand out until capture all mass
  ## The idea is that since we normalized by \tau, \pm 5 should work
  ##   most of the time, but maybe we need \pm 10 (hand-wave Chebyshev).
  ##   We expand around modes until we get all of the mass.
  totmass <- 0
  mid <- mean(modes)
  low1 <- modes[[1]] - 3
  up1 <- min(modes[[1]] + 3, mid)
  low2 <- max(modes[[2]] - 3, mid)
  up2 <- modes[[2]] + 3
  step <- 2
  while(abs(totmass - 1) > 1e-3) {
    low1 <- low1 - step
    up1 <- min(up1 + step, mid)
    low2 <- max(low2 - step, mid)
    up2 <- up2 + step
    mode1_mass <- stats::integrate(f = f, lower = low1, upper = up1, subdivisions = subdivisions)[[1]]
    mode2_mass <- stats::integrate(f = f, lower = low2, upper = up2, subdivisions = subdivisions)[[1]]
    totmass <- mode1_mass + mode2_mass

    if (totmass < 1e-6) {
      stop("Couldn't locate mass")
    }
  }

  ## Just to be safe, do it two more times
  low1 <- low1 - step * 2
  up1 <- min(up1 + step * 2, mid)
  low2 <- max(low2 - step * 2, mid)
  up2 <- up2 + step * 2

  ## Precompute this one because it is used whenever q >= 0
  if (any(q >= 0)) {
    ## Scenario: bounds
    ## (l,u,0) : (l,u)
    ## (l,0,u) : (l,0)
    ## (0,l,u) : (0,0)
    ## so bounds of integration are min(l,0) to min(0,u)
    negval <- stats::integrate(f, lower = min(low1, 0), upper = min(up1, 0), subdivisions = subdivisions)$value +
      stats::integrate(f, lower = min(low2, 0), upper = min(up2, 0), subdivisions = subdivisions)$value
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
      # If q > 0, then integral is from -Inf to 0 and from 1/(tau * q) to Inf
      ## Scenario: bounds (a = 1/(tau * q))
      ## (l,u,a) : (a,a)
      ## (l,a,u) : (a,u)
      ## (a,l,u) : (l,u)
      ## so bounds of integration are max(l,a) to max(u,a)
      tq_bound <- 1 / (tau * q[[i]])
      lpvec[[i]] <- negval +
        stats::integrate(f = f, lower = max(low1, tq_bound), upper = max(up1, tq_bound), subdivisions = subdivisions)$value +
        stats::integrate(f = f, lower = max(low2, tq_bound), upper = max(up2, tq_bound), subdivisions = subdivisions)$value
    } else {
      # If q < 0, then integral is from 1/(tau * q) to 0
      ## Scenario  : bounds (a = 1/(tau * q))
      ## (l,u,a,0) : skipped
      ## (l,a,u,0) : (a,u)
      ## (l,a,0,u) : (a,0)
      ## (a,l,u,0) : (l,u)
      ## (a,l,0,u) : (l,0)
      ## (a,0,l,u) : skipped
      ## so bounds of integration are max(l,a) to min(u,0)
      tq_bound <- 1 / (tau * q[[i]])
      lpvec[[i]] <- 0
      if (!(up1 < tq_bound || 0 < low1)) {
        lpvec[[i]] <- lpvec[[i]] +
          stats::integrate(f = f, lower = max(low1, tq_bound), upper = min(up1, 0), subdivisions = subdivisions)$value
      }
      if (!(up2 < tq_bound || 0 < low2)) {
        lpvec[[i]] <- lpvec[[i]] +
          stats::integrate(f = f, lower = max(low2, tq_bound), upper = min(up2, 0), subdivisions = subdivisions)$value
      }
    }
  }

  if (isTRUE(lower.tail)) {
    ret <- lpvec
  } else {
    ret <- 1 - lpvec
  }

  return(ret)
}

pginvnorm_old <- function(q, alpha, mu = 0, tau = 1, lower.tail = TRUE) {
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

  if (alpha == 2) {
    return(qinvnorm(p = p, imean = mu, isd = tau))
  }

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
      f <- function(x) (pginvnorm(q = x, alpha = alpha, mu = mu, tau = tau) - p[[i]]) * 1e3
      qvec[[i]] <- stats::uniroot(f = f, interval = c(lower, upper), check.conv = TRUE, tol = .Machine$double.eps^(0.75))$root
    }
  }

  return(qvec)
}

#' @describeIn ginvnorm Random generation.
#'
#' @export
rginvnorm <- function(n, alpha, mu = 0, tau = 1) {
  if (alpha == 2) {
    return(rinvnorm(n = n, imean = mu, isd = tau))
  }
  return(qginvnorm(p = stats::runif(n = n), alpha = alpha, mu = mu, tau = tau))
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
