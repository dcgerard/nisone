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

#' Bayesian credible interval when n >= 1
#'
#' Let \eqn{X_i \sim N(\mu, \sigma^2)}. Let \eqn{Y_i = X_i - A},
#' \eqn{\beta = \mu - A}, and \eqn{\nu = |\mu - A|/\sigma}. Then,
#' for a fixed value of \eqn{\nu}, this comes up with a \eqn{(1-\alpha)100%}
#' posterior credible interval of \eqn{\mu} where the prior is
#' \eqn{\pi(\beta) = |\beta|^{-1}} or \eqn{\pi(\mu) = |\mu - A|^{-1}}.
#' The full posterior distribution is generalized inverse normal (implemented
#' by \code{\link{ginvnorm}}) with shape parameter \eqn{n + 1}, inverse mean of
#' \code{sum(y) / sum(y^2)}, and inverse variance \code{1 / (nu^2 * sum(y^2))}.
#'
#' @param x a double
#' @param A The prior center.
#' @param nu The fixed value of nu.
#' @param level The level of the credible interval
#'
#' @return The credible interval of the provided level.
#'
#' @author David Gerard
#'
#' @export
bcin <- function(x, A = 0, nu = 1, level = 0.95) {
  y <- x - A

  if (is.null(nu)) {
    stop("Not ready for prime time yet")
    ## Exact values for n <= 10 from optimal_nu.R
    avec <- c(0.66154724673845, 0.866025840852927, 1.09736102743933, 1.34058742733214,
              1.58809677440661, 1.83703551422352, 2.08646555136187, 2.3360888370321,
              2.58580619053097)
    n <- length(x)
    if (n <= 10) {
      a <- avec[[n - 1]]
    } else {
      ## Asymptotic approximation
      a <- n / 4 + 1/12 + 2 / (81 * n)
    }
    nu <- sqrt(a * 2 * sum(y^2) / sum(y)^2)
  }

  alpha <- length(x) + 1
  mu <- sum(y) / sum(y^2)
  tau <- 1 / (nu * sqrt(sum(y^2)))
  al <- (1 - level) / 2
  retvec <- qginvnorm(p = c(al, 1 - al), alpha = alpha, mu = mu, tau = tau) + A
  names(retvec) <- c("lower", "upper")
  return(retvec)
}

#' Generalized inverse tail approximation.
#'
#' This function is the tail approximation via a Taylor series.
#' This is an interval function just for sanity checking my math.
#' It should approximate the tail for the generalized inverse normal.
#'
#' @param q The quantile.
#' @param alpha shape
#' @param mu inverse mean
#' @param tau inverse sd
#'
#' @return The approximation for Pr(X >= q) when X ~ Ginvnorm(alpha, mu, tau).
#' Should be accurate for large q.
#'
#' @author David Gerard
#'
#' @noRd
ginv_tail <- function(q, alpha, mu, tau, log = FALSE) {
  ret <- .gin_log_K(alpha = alpha, mu = mu, tau = tau) - mu^2 / (2 * tau^2) - log(alpha - 1) - (alpha - 1) * log(q)
  if (!log) {
    ret <- exp(ret)
  }
  return(ret)
}



