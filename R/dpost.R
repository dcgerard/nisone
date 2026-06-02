#' Posterior Distribution for n=1 Data
#'
#' Density, probability, quantile, and random generation from the posterior
#' distribution of the location parameter when one has \eqn{n = 1} observation
#' and uses a particular prior that results in valid confidence intervals
#' (asymptotically in the confidence level).
#'
#' Let \eqn{X} have PDF \eqn{\frac{1}{\sigma}\rho((x - \mu)/\sigma)} for a
#' known and symmetric \eqn{\rho(\cdot)}. Let \eqn{\nu = |\mu - A| / \sigma}
#' for a pre-specified \eqn{A}. We place prior \eqn{\pi(\mu) = |\mu - A|^{-1}}
#' and \eqn{\nu = \text{argmax}_{a > 0}a\rho(a)} with probability 1.
#'
#' These functions will calculate the posterior distribution and allow you to
#' interact with it through the distribution, density, quantile, and random
#' generation functions.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n sample size
#' @param A The prior value.
#' @param obs The observed value. This is \eqn{X} in the description.
#' @param nu The prior value of nu. If not known, use \code{nun1post()}
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P(X<=x)
#'     otherwise, P(X>x).
#' @param fam One of \code{"normal"} or \code{"cauchy"}. If
#'     \code{ddist}, \code{pdist}, \code{qdist}, or \code{rdist} are specified
#'     then this argument is ignored.
#' @param ddist Density function of standard distribution \eqn{\rho()}. Should have a `log` argument.
#' @param pdist Cumulative distribution function of standard distribution.
#' @param qdist Quantile function of standard distribution.
#' @param rdist Random generation of standard distribution.
#' @param ... Additional arguments for \code{ddist}, \code{pdist},
#'     \code{qdist}, and \code{rdist}.
#'
#' @name n1post
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' # Observe x = 2, assume t with 2 df, prior value is A = 1
#' nun1post(ddist = dt, df = 2) ## nu should be 1
#' x <- seq(-10, 10, length.out = 500)
#' y <- dn1post(x = x, A = 1, obs = 2, nu = 1, ddist = dt, df = 2)
#' z <- rn1post(n = 10000, A = 1, obs = 2, nu = 1, rdist = rt, df = 2)
#' z <- z[z >= -10 & z <= 10]
#' graphics::hist(z, freq = FALSE, breaks = 200, border = "grey", col = "grey")
#' graphics::lines(x, y, type = "l", col = "#E69F00")
#' graphics::abline(v = c(1, 2), col = c("#56B4E9", "#009E73"), lty = c(2, 3))
#' pn1post(q = 2, A = 1, obs = 2, nu = 1, pdist = pt, df = 2)
#' qn1post(p = 0.7113, A = 1, obs = 2, nu = 1, qdist = qt, pdist = pt, df = 2)
#'
NULL

#' @describeIn n1post Obtains prior value of \eqn{\nu}.
#'
#' @export
nun1post <- function(fam = c("normal", "cauchy"), ddist = NULL, ...) {
  fam <- match.arg(fam)
  if (!is.null(ddist)) {
    f <- function(a) {
      a * ddist(a, ...)
    }
    oout <- stats::optim(
      par = 1,
      fn = f,
      lower = 0,
      upper = Inf,
      control = list(fnscale = -1),
      method = "L-BFGS-B"
      )
    nu <- oout$par
  } else if (fam == "normal" || fam == "cauchy") {
    nu <- 1
  }
  return(nu)
}

#' @describeIn n1post Density function.
#'
#' @export
dn1post <- function(x, A, obs, nu = NULL, fam = c("normal", "cauchy"), ddist = NULL, log = FALSE, ...) {
  fam <- match.arg(fam)
  if (is.null(nu)) {
    nu <- nun1post(fam = fam, ddist = ddist, ...)
  }

  center <- 1 / (obs - A)
  scale <- 1 / abs(nu * (obs - A))

  didist(x = x - A, center = center, scale = scale, fam = fam, ddist = ddist, log = log, ...)
}

#' @describeIn n1post Distribution function.
#'
#' @export
pn1post <- function(q, A, obs, nu = NULL, fam = c("normal", "cauchy"), pdist = NULL, lower.tail = TRUE, log.p = FALSE, ...) {
  fam <- match.arg(fam)
  if (is.null(nu) && !is.null(pdist)) {
    stop("If you specify pdist, you need to specify nu. Use nun1post() to get nu")
  }

  center <- 1 / (obs - A)
  scale <- 1 / abs(nu * (obs - A))

  pidist(q = q - A, center = center, scale = scale, fam = fam, pdist = pdist, lower.tail = lower.tail, log.p = log.p, ...)
}

#' @describeIn n1post Quantile function.
#'
#' @export
qn1post <- function(p, A, obs, nu = NULL, fam = c("normal", "cauchy"), qdist = NULL, pdist = NULL, ...) {
  fam <- match.arg(fam)
  if (is.null(nu) && (!is.null(qdist) || !is.null(pdist))) {
    stop("If you specify qdist or pdist, you need to specify nu. Use nun1post() to get nu")
  }

  center <- 1 / (obs - A)
  scale <- 1 / abs(nu * (obs - A))

  qidist(p = p, center = center, scale = scale, fam = fam, qdist = qdist, pdist = pdist, ...) + A
}

#' @describeIn n1post Random generation.
#'
#' @export
rn1post <- function(n, A, obs, nu = NULL, fam = c("normal", "cauchy"), rdist = NULL, ...) {
  fam <- match.arg(fam)
  if (is.null(nu) && !is.null(rdist)) {
    stop("If you specify rdist, you need to specify nu. Use nun1post() to get nu")
  }

  center <- 1 / (obs - A)
  scale <- 1 / abs(nu * (obs - A))

  ridist(n = n, center = center, scale = scale, fam = fam, rdist = rdist, ...) + A
}
