## Normal distribution n=1 confidence intervals

arcoth <- function(x) {
  0.5 * log((x + 1) / (x - 1))
}

#' Worst-case coefficient of variation
#'
#' This assumes the random variable X is normal. Intervals are of the form
#' center +/- width * |X|.
#'
#' We use A = 0 for calculations, but the worst
#' case COV doesn't change if A is non-zero. So you can think of
#' intervals of the form center +/- width * |X - A| where
#' center is either X or (X + A)/2.
#'
#' @param width The half-width of the interval, in units of |X|.
#'      This needs to be at least 1 if \code{center = "X"}, or
#'      at least 0.5 if \code{center = "ave"}.
#' @param center What is the center of the interval? Either \code{"X"} or
#'     \code{"ave"}.
#'
#' @author David Gerard
#'
#' @returns Worst case coefficient of variation \eqn{\nu = |X - A|/\sigma}.
#'
#' @examples
#' wc_cov(width = 2, center = "X")
#' wc_cov(width = 2, center = "ave")
#'
#' @export
wc_cov <- function(width, center = c("X", "ave")) {
  center <- match.arg(center)
  if (center == "X") {
    stopifnot(width >= 1)
    if (width == 1) {
      return(0)
    }
    nu <- (1 - 1 / width^2) * sqrt(width * arcoth(width))
  } else if (center == "ave") {
    stopifnot(width >= 0.5)
    if (width == 0.5) {
      return(0)
    }
    nu <- (4 * width^2 - 1) * sqrt(arcoth(2 * width) / (2 * (4 * width^3 + width)))
  } else {
    stop("here")
  }
  return(nu)
}

#' Worst case exclusion probability of a n=1 confidence interval
#'
#' This assumes the random variable X follows a normal distribution.
#'
#' @inheritParams wc_cov
#'
#' @returns The worst case probability (a priori) that a confidence interval
#'      of half-width \code{width} will not capture the mean.
#'
#' @examples
#' wc_alpha(width = 2, center = "X")
#' wc_alpha(width = 2, center = "ave")
#'
#' @author David Gerard
#'
#' @export
wc_alpha <- function(width, center = c("X", "ave")) {
  center <- match.arg(center)
  nu <- wc_cov(width = width, center = center)
  if (center == "X") {
    if (width == 1) {
      return(0.5)
    }
    alpha <- stats::pnorm(nu * width / (width - 1)) -
      stats::pnorm(nu * width / (width + 1))
  } else if (center == "ave") {
    if (width == 0.5) {
      return(0.5)
    }
    alpha <- stats::pnorm(nu * (2 * width + 1) / (2 * width - 1)) -
      stats::pnorm(nu * (2 * width - 1) / (2 * width + 1))
  } else {
    stop("here")
  }
  return(alpha)
}

#' Given confidence level, provide half-width of CI
#'
#' \eqn{(1-\alpha)100\%}
#' confidence intervals are of the form
#' \eqn{X \pm \eta|X-A|} or \eqn{(X + A)/2 \pm \eta|X-A|}. You give \eqn{\alpha}
#' and this function will provide \eqn{\eta} if X follows a Cauchy,
#' normal, or uniform distribution.
#'
#' @param alpha The confidence error probability. We produce a (1-alpha)100 percent
#'     confidence interval.
#' @param center Either X, ave, or the asymptotic width (approx).
#' @param family Either the normal, Cauchy, or uniform distribution.
#'
#' @returns Returns half-width of worst-case confidence interval in units of |X - A|.
#'
#' @examples
#' wc_width(alpha = 0.2, center = "X")
#' wc_width(alpha = 0.2, center = "ave")
#' wc_width(alpha = 0.2, center = "approx")
#'
#' wc_width(alpha = 0.05, center = "X")
#' wc_width(alpha = 0.05, center = "ave")
#' wc_width(alpha = 0.05, center = "approx")
#'
#' @author David Gerard
#'
#' @export
wc_width <- function(alpha, center = c("X", "ave", "approx"), family = c("normal", "cauchy", "uniform")) {
  stopifnot(length(alpha) == 1,
            alpha <= 0.5,
            alpha >= 0)
  if (alpha == 0) {
    return(Inf)
  }
  center <- match.arg(center)
  family <- match.arg(family)

  if (family == "cauchy") {
    if (center == "X") {
      width <- 1 / sin(pi * alpha)
    } else if (center == "ave") {
      width <- 1 / (2 * tan(pi * alpha / 2))
    } else if (center == "approx") {
      width <- 1 / (pi * alpha)
    }
  } else if (family == "uniform") {
    if (center == "X" | center == "approx") {
      width <- 1 / alpha - 1
    } else if (center == "ave") {
      width <- 1 / (2 * alpha) - 0.5 + sqrt(1 / (2 * alpha)^2 - 1 / (2 * alpha))
    }
  } else if (family == "normal") {
    if (center == "approx" || alpha < 0.01) {
      width <- sqrt(2 / (pi * exp(1))) / alpha
    } else if (center == "X") {
      if (alpha == 0.5) {
        return(1)
      }
      f <- function(x) {
        alpha - wc_alpha(width = x, center = "X")
      }
      rout <- stats::uniroot(f = f, c(1, 50))
      width <- rout$root
    } else if (center == "ave") {
      if (alpha == 0.5) {
        return(0.5)
      }
      f <- function(x) {
        alpha - wc_alpha(width = x, center = "ave")
      }
      rout <- stats::uniroot(f = f, c(0.5, 50))
      width <- rout$root
    }
  } else {
    stop("Not a valid family")
  }
  return(width)
}

#' Confidence distribution CDF based on n=1 confidence interval
#'
#' Let \eqn{X \sim N(\mu, \sigma^2)}. Consider intervals of the form
#' \eqn{X \pm \eta |X|} or \eqn{X / 2 \pm \eta |X|}. These intervals
#' define a confidence distribution of \eqn{\mu}. This is the
#' CDF of the confidence distribution of \eqn{(\mu - X)/|X|}
#' or \eqn{(\mu - X/2)/|X|}. Please note that \eqn{\mu} is random
#' here not \eqn{X}. Also note that this confidence distribution does
#' not exist between the 25th and 75th percentiles.
#'
#' @param q The quantile. Only defined for \code{abs(q) >= 1} when
#'     \code{center = "X"}, and for \code{abs(q) >= 0.5} when
#'     \code{center = "ave"}.
#' @param center What is the center of the interval? Either \code{"X"} or
#'     \code{X/2}, (\code{"ave"}).
#'
#' @author David Gerard
#'
#' @returns n=1 confidence distribution CDF.
#'
#' @seealso [d_wc()]
#'
#' @examples
#' qseq <- seq(-10, 10, length.out = 100)
#' pseq <- p_wc(qseq)
#' plot(qseq, pseq, type = "l", xlab = "q", ylab = "CDF")
#'
#' @export
p_wc <- function(q, center = c("X", "ave")) {
  center <- match.arg(center)
  p <- rep(NA_real_, length.out = length(q))
  if (center == "X") {
    p[abs(q) < 1] <- NA_real_
    p[q == -1] <- 0.25
    p[q == 1] <- 0.75
    neg <- q < -1
    pos <- q > 1
    p[neg | pos] <- 0.5 * stats::pnorm((q[neg | pos] + 1) * sqrt(arcoth(q[neg | pos]) / q[neg | pos])) -
      0.5 * stats::pnorm((q[neg | pos] - 1) * sqrt(arcoth(q[neg | pos]) / q[neg | pos]))
    p[pos] <- 1 - p[pos]
  } else if (center == "ave") {
    p[abs(q) < 0.5] <- NA_real_
    p[q == -0.5] <- 0.25
    p[q == 0.5] <- 0.75
    neg <- q < -0.5
    pos <- q > 0.5
    p[neg | pos] <- 0.5 * stats::pnorm((2 * q[neg | pos] - 1)^2 * sqrt(arcoth(2 * q[neg | pos]) / (2 * q[neg | pos] * (4 * q[neg | pos]^2 + 1)))) -
      0.5 * stats::pnorm((2 * q[neg | pos] + 1)^2 * sqrt(arcoth(2 * q[neg | pos]) / (2 * q[neg | pos] * (4 * q[neg | pos]^2 + 1))))
    p[pos] <- 1 + p[pos]
  } else {
    stop("here")
  }
  return(p)
}

#' Confidence density based on n=1 confidence interval
#'
#' Let \eqn{X \sim N(\mu, \sigma^2)}. Consider intervals of the form
#' \eqn{X \pm \eta |X|} or \eqn{X / 2 \pm \eta |X|}. These intervals
#' define a confidence distribution of \eqn{\mu}. This is the
#' density of the confidence distribution of \eqn{(\mu - X)/|X|}
#' or \eqn{(\mu - X/2)/|X|}. Please note that \eqn{\mu} is random
#' here not \eqn{X}. Also note that this confidence distribution does
#' not exist between the 25th and 75th percentiles.
#'
#' @param x The value at which to calculate the density.
#'     Only defined for \code{abs(x) >= 1} when
#'     \code{center = "X"}, and for \code{abs(x) >= 0.5} when
#'     \code{center = "ave"}.
#' @param center What is the center of the interval? Either \code{"X"} or
#'     \code{"X/2"} (\code{"ave"}).
#'
#' @author David Gerard
#'
#' @returns n=1 confidence density.
#'
#' @seealso [p_wc()]
#'
#' @examples
#' xseq <- seq(-10, 10, length.out = 500)
#' dseq <- d_wc(xseq)
#' plot(xseq, dseq, type = "l")
#'
#' dseq <- d_wc(xseq, center = "ave")
#' plot(xseq, dseq, type = "l")
#'
#' @export
d_wc <- function(x, center = c("X", "ave")) {
  center <- match.arg(center)
  if (center == "X") {
    x[abs(x) <= 1] <- NA_real_
    dvec <- ((x + 1) / (x - 1))^(-(x - 1)^2 / (4 * x)) * sqrt(arcoth(x) / x) / (sqrt(2 * pi) * abs(x + 1))
  } else if (center == "ave") {
    x <- abs(x)
    x[abs(x) <= 0.5] <- NA_real_
    dvec <-
      -((-(1 + 2 * x)^2 * ((1 + 2 * x) / (-1 + 2 * x))^(-((1 + 2 * x)^4 / (
        8 * (x + 4 * x^3)))) * (-4 * (x + 4 * x^3) + (1 - 2 * x)^4 * log((
        1 + 2 * x) / (-1 + 2 * x))) + (1 - 2 * x)^2 * ((
        1 + 2 * x) / (-1 + 2 * x))^(-((1 - 2 * x)^4 / (
        8 * (x + 4 * x^3)))) * (-4 * (x + 4 * x^3) + (1 + 2 * x)^4 * log((
        1 + 2 * x) / (-1 + 2 * x)))) / (8 * sqrt(2 * pi) * sqrt(
        x + 4 * x^3
      ) * (x - 16 * x^5) * sqrt(log((1 + 2 * x) / (-1 + 2 * x)))))
  } else {
    stop("here")
  }
  return(dvec)
}

#' n=1 confidence interval
#'
#' When you have one observation, \eqn{X} from some symmetric distribution
#' with center \eqn{\mu}, for a pre-specified \eqn{A}, produces
#' confidence intervals for the center form  \eqn{X \pm \eta |X - A|}
#' (\code{type = "x"}) or \eqn{(X + A)/2 \pm \eta |X - A|}
#' (\code{type = "ave"}). The average intervals have smaller width on
#' average, and so are the default. We allow \eqn{X} to either come from
#' a normal distribution, a Cauchy, or a uniform. Note that if you use
#' \code{family = "uniform"} then these intervals are valid confidence intervals
#' for the median/mode for \emph{any} symmetric unimodal density, which
#' I think is really amazing.
#'
#' @param x A vector of single observations.
#' @param A Your "prior knowledge" or "augmented data value".
#' @param type Either centered at x or the average of x and A
#' @param family Either the normal distribution, Cauchy, or uniform.
#' @param level The level of the confidence interval.
#'
#' @examples
#' ci1(c(1, 2, 10), type = "x")
#' ci1(c(1, 2, 10), type = "ave")
#'
#' @returns
#' A matrix of with 4 columns: the data, the center, the lower bound,
#' the upper bound.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'  \item{Blachman, N., & Machol, R. (1987). Confidence intervals based on one or more observations. \emph{IEEE Transactions on Information Theory}, 33(3), 373-382. \doi{10.1109/TIT.1987.1057306}}
#' }
#'
#' @export
ci1 <- function(x, A = 0, type = c("ave", "x"), family = c("normal", "cauchy", "uniform"), level = 0.95) {
  type <- match.arg(type)
  family <- match.arg(family)
  alpha <- 1 - level

  if (type == "x") {
    width <- wc_width(alpha = alpha, center = "X", family = family)
    diff <- abs(x - A)
    lower <- x - width * diff
    upper <- x + width * diff
    center <- x
  } else if (type == "ave") {
    width <- wc_width(alpha = alpha, center = "ave", family = family)
    diff <- abs(x - A)
    lower <- (x + A) / 2 - width * diff
    upper <- (x + A) / 2 + width * diff
    center <- (x + A) / 2
  }
  return(cbind(x = x, center = center, lower = lower, upper = upper))
}
