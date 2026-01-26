## Normal distribution

arcoth <- function(x) {
  0.5 * log((x + 1) / (x - 1))
}

#' Worst-case coefficient of variation
#'
#' This assumes the random variable X is normal. Intervals are of the form
#' center +/- width * |X|.
#'
#' @param width The half-width of the interval, in units of |X|.
#'      This needs to be at least 1 if \code{center = "X"}, or
#'      at least 0.5 if \code{center = "X/2"}.
#' @param center What is the center of the interval? Either \code{"X"} or
#'     \code{"X/2"}.
#'
#' @author David Gerard
#'
#' @examples
#' wc_cov(width = 2, center = "X")
#' wc_cov(width = 2, center = "X/2")
#'
#' @export
wc_cov <- function(width, center = c("X", "X/2")) {
  center <- match.arg(center)
  if (center == "X") {
    stopifnot(width >= 1)
    if (width == 1) {
      return(0)
    }
    nu <- (1 - 1 / width^2) * sqrt(width * arcoth(width))
  } else if (center == "X/2") {
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
#' @return The worst case probability (a priori) that a confidence interval
#'      of half-width \code{width} will not capture the mean.
#'
#' @examples
#' wc_alpha(width = 2, center = "X")
#' wc_alpha(width = 2, center = "X/2")
#'
#' @author David Gerard
#'
#' @export
wc_alpha <- function(width, center = c("X", "X/2")) {
  center <- match.arg(center)
  nu <- wc_cov(width = width, center = center)
  if (center == "X") {
    if (width == 1) {
      return(0.5)
    }
    alpha <- stats::pnorm(nu * width / (width - 1)) -
      stats::pnorm(nu * width / (width + 1))
  } else if (center == "X/2") {
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
#' This assumes X follows a normal distribution. Intervals are of the form
#' X +/- t|X| or X/2 +/- t|X|.
#'
#' @param alpha The exclusion probability. We produce a (1-alpha)100 percent
#'     confidence interval.
#' @param center Either X, X/2, or the asymptotic width (approx).
#'
#' @return Returns half-width of worst-case confidence interval in units of |X|.
#'
#' @examples
#' wc_width(alpha = 0.2, center = "X")
#' wc_width(alpha = 0.2, center = "X/2")
#' wc_width(alpha = 0.2, center = "approx")
#'
#' wc_width(alpha = 0.05, center = "X")
#' wc_width(alpha = 0.05, center = "X/2")
#' wc_width(alpha = 0.05, center = "approx")
#'
#' @author David Gerard
#'
#' @export
wc_width <- function(alpha, center = c("X", "X/2", "approx")) {
  stopifnot(length(alpha) == 1,
            alpha <= 0.5,
            alpha >= 0)
  if (alpha == 0) {
    return(Inf)
  }
  center <- match.arg(center)
  if (center == "approx" || alpha < 0.01) {
    width <- sqrt(2 / (pi * exp(1))) / alpha
  }
  else if (center == "X") {
    if (alpha == 0.5) {
      return(1)
    }
    f <- function(width) {
      alpha - wc_alpha(width = width, center = "X")
    }
    rout <- stats::uniroot(f = f, c(1, 50))
    width <- rout$root
  } else if (center == "X/2") {
    if (alpha == 0.5) {
      return(0.5)
    }
    f <- function(width) {
      alpha - wc_alpha(width = width, center = "X/2")
    }
    rout <- stats::uniroot(f = f, c(0.5, 50))
    width <- rout$root
  }
  return(width)
}

#' Implied posterior CDF based on n=1 confidence interval
#'
#' This is the CDF of (mu-X) / |X|, which is a pivotal quantity.
#'
#' @param q The quantile. Only defined for abs(q) >= 1
#' @param center What is the center of the interval? Either \code{"X"} or
#'     \code{"X/2"}.
#'
#' @author David Gerard
#'
#' @examples
#' qseq <- seq(-10, 10, length.out = 100)
#' pseq <- p_wc(qseq)
#' plot(qseq, pseq, type = "l", xlab = "q", ylab = "CDF")
#'
#' @export
p_wc <- function(q, center = c("X", "X/2")) {
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
  } else if (center == "X/2") {
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

#' Implied posterior density based on n=1 confidence interval
#'
#' This is the density of (mu-X) / |X|, which is a pivotal quantity.
#'
#' @param x The value to evaluate the density. Only defined for abs(q) > 1
#' @param center What is the center of the interval? Either \code{"X"} or
#'     \code{"X/2"}.
#'
#' @author David Gerard
#'
#' @examples
#' xseq <- seq(-10, 10, length.out = 500)
#' dseq <- d_wc(xseq)
#' plot(xseq, dseq, type = "l")
#'
#' dseq <- d_wc(xseq, center = "X/2")
#' plot(xseq, dseq, type = "l")
#'
#' @export
d_wc <- function(x, center = c("X", "X/2")) {
  center <- match.arg(center)
  if (center == "X") {
    x[abs(x) <= 1] <- NA_real_
    dvec <- ((x + 1) / (x - 1))^(-(x - 1)^2 / (4 * x)) * sqrt(arcoth(x) / x) / (sqrt(2 * pi) * abs(x + 1))
  } else if (center == "X/2") {
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
#' @param x A vector of single observations.
#' @param A Where the CI should be centered
#' @param type Either centered at x or the average of x and A
#' @param level The level of the confidence interval.
#'
#' @examples
#' ci1(c(1, 2, 10), type = "x")
#' ci1(c(1, 2, 10), type = "ave")
#'
#' @references
#' \itemize{
#'   \item{Blachman, N., & Machol, R. (1987). Confidence intervals based on one or more observations. IEEE transactions on information theory, 33(3), 373-382.}
#' }
#'
#' @author David Gerard
#'
#' @export
ci1 <- function(x, A = 0, type = c("ave", "x"), level = 0.95) {
  type <- match.arg(type)
  alpha <- 1 - level
  if (type == "x") {
    width <- wc_width(alpha = alpha, center = "X")
    diff <- abs(x - A)
    lower <- x - width * diff
    upper <- x + width * diff
    center <- x
  } else if (type == "ave") {
    width <- wc_width(alpha = alpha, center = "X/2")
    diff <- abs(x - A)
    lower <- (x + A) / 2 - width * diff
    upper <- (x + A) / 2 + width * diff
    center <- (x + A) / 2
  }
  return(cbind(x = x, center = center, lower = lower, upper = upper))
}

