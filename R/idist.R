#' General inverse distribution
#'
#' Let \eqn{X} come from some location-scale family
#' \deqn{f(x) = \frac{1}{\sigma}\rho\left(\frac{X - \mu}{\sigma}\right),}
#' for some standard distribution \eqn{\rho()}. You provide \eqn{\rho()},
#' \eqn{\mu} and \eqn{\sigma} and these functions will provide the density,
#' distribution, quantile, and random generation for \eqn{Z = 1/X}. Convenience
#' arguments for the normal, Cauchy, and uniform are provided.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n sample size
#' @param center The center parameter \eqn{\mu}.
#' @param scale The scale parameter \eqn{\sigma}.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P(X<=x)
#'     otherwise, P(X>x).
#' @param fam One of \code{"normal"}, \code{"cauchy"}, or \code{"uniform"}. If
#'     \code{ddist}, \code{pdist}, \code{qdist}, or \code{rdist} are specified
#'     then this argument is ignored.
#' @param ddist Density function of standard distribution \eqn{\rho()}. Should have a `log` argument.
#' @param pdist Cumulative distribution function of standard distribution.
#' @param qdist Quantile function of standard distribution.
#' @param rdist Random generation of standard distribution.
#' @param ... Additional arguments for \code{ddist}, \code{pdist},
#'     \code{qdist}, and \code{rdist}.
#'
#' @return Either a random sample (\code{ridist}),
#'     the density (\code{didist}), the tail
#'     probability (\code{pidist}), or the quantile
#'     (\code{qidist}) of the inverse distribution.
#'
#' @author David Gerard
#'
#' @name idist
#'
#' @examples
#' didist(1, fam = "normal")
#' didist(1, center = 3, scale = 2, ddist = dt, df = 1)
#'
#' pidist(1, fam = "normal")
#' pidist(1, center = 3, scale = 2, pdist = pt, df = 1)
#'
#' qidist(0.1, fam = "normal")
#' qidist(0.1, center = 3, scale = 2, qdist = qt, pdist = pt, df = 1)
#'
NULL

#' @describeIn idist Density function.
#'
#' @export
didist <- function(x, center = 0, scale = 1, fam = c("normal", "cauchy", "uniform"), ddist = NULL, log = FALSE, ...) {
  if (is.null(ddist)) {
    fam <- match.arg(fam)
    if (fam == "normal") {
      dfun <- function(x) {
        stats::dnorm(x = x, mean = center, sd = scale, log = TRUE)
      }
    } else if (fam == "cauchy") {
      dfun <- function(x) {
        stats::dt(x = (x - center) / scale, df = 1, log = TRUE) - log(scale)
      }
    } else if (fam == "uniform") {
      dfun <- function(x) {
        stats::dunif(x = (x - center) / scale, min = -0.5, max = 0.5, log = TRUE) - log(scale)
      }
    } else {
      stop("Not a supported distribution")
    }
  } else {
    dfun <- function(x) {
      ddist(x = (x - center) / scale, log = TRUE, ...) - log(scale)
    }
  }

  ret <- dfun(1 / x) - 2 * log(abs(x))
  if (!log) {
    ret <- exp(ret)
  }
  return(ret)
}

#' @describeIn idist Distribution function.
#'
#' @export
pidist <- function(q, center = 0, scale = 1, fam = c("normal", "cauchy", "uniform"), pdist = NULL, lower.tail = TRUE, log.p = FALSE, ...) {
  if (is.null(pdist)) {
    fam <- match.arg(fam)
    if (fam == "normal") {
      pfun <- stats::pnorm
    } else if (fam == "cauchy") {
      pfun <- function(q) {
        stats::pt(q = q, df = 1)
      }
    } else if (fam == "uniform") {
      pfun <- function(q) {
        stats::punif(q = q, min = -0.5, max = 0.5)
      }
    }
  } else {
    pfun <- function(q) {
      pdist(q = q, ...)
    }
  }
  pvec <- pfun((center * q - 1) / (scale * q)) - pfun(center / scale)
  pvec[q > 0] <- pvec[q > 0] + 1
  pvec[q == 0] <- 1 - pfun(center / scale)

  if (!lower.tail) {
    pvec <- 1 - pvec
  }
  if (log.p) {
    pvec <- log(pvec)
  }

  return(pvec)
}

#' @describeIn idist Quantile function.
#'
#' @export
qidist <- function(p, center = 0, scale = 1, fam = c("normal", "cauchy", "uniform"), qdist = NULL, pdist = NULL, ...) {

  if (xor(is.null(pdist), is.null(qdist))) {
    stop("Either both pdist and qdist need to be specified, or both need to be NULL")
  }

  if (is.null(qdist) || is.null(pdist)) {
    if (length(fam) > 1) {
      message("Using fam = 'normal'")
    }
    fam <- match.arg(fam)
    if (fam == "normal") {
      qfun <- stats::qnorm
      pfun <- stats::pnorm
    } else if (fam == "cauchy") {
      qfun <- function(p) {
        stats::qt(p = p, df = 1)
      }
      pfun <- function(q) {
        stats::pt(q = q, df = 1)
      }
    } else if (fam == "uniform") {
      qfun <- function(p) {
        stats::qunif(p = p, min = -0.5, max = 0.5)
      }
      pfun <- function(q) {
        stats::punif(q = q, min = -0.5, max = 0.5)
      }
    } else {
      stop("Not a supported family")
    }
  } else {
    qfun <- function(p) {
      qdist(p = p, ...)
    }
    pfun <- function(q) {
      pdist(q = q, ...)
    }
  }

  cutval <- pfun(center / scale)
  pvec <- rep(NA_real_, length.out = length(p))
  pvec[p < 1 - cutval] <- 1 / (center - scale * qfun(p[p < 1 - cutval] + cutval))
  pvec[p >= 1 - cutval] <- 1 / (center - scale * qfun(p[p >= 1 - cutval] + cutval - 1))
  return(pvec)
}

#' @describeIn idist Random generation.
#'
#' @export
ridist <- function(n, center = 0, scale = 1, fam = c("normal", "cauchy", "uniform"), rdist = NULL, ...) {
  if (is.null(rdist)) {
    fam <- match.arg(fam)
    if (fam == "normal") {
      rfun <- stats::rnorm
    } else if (fam == "cauchy") {
      rfun <- function(n) {
        stats::rt(n = n, df = 1)
      }
    } else if (fam == "uniform") {
      rfun <- function(n) {
        stats::runif(n = n, min = -0.5, max = 0.5)
      }
    } else {
      stop("Not a supported family")
    }
  } else {
    rfun <- function(n) {
      rdist(n = n, ...)
    }
  }

  return(1 / (rfun(n = n) * scale + center))
}
