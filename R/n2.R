#' w as a function of z
#'
#' @param z standard normal value
#' @param n sample size
#' @param eta2 multiplier squared
#' @param nu (mu - A) / sigma
#'
#' @author David Gerard
#'
#' @noRd
w_z <- function(z, n, eta2, nu) {
  ## ((n * z - sqrt(n) * nu)^2 - eta2 * (z + sqrt(n) * nu)^2) / ((n + 1) * eta2)
  ((n^2 - eta2) * z^2 - 2 * (n^(1.5) + sqrt(n) * eta2) * nu * z + (1 - eta2) * n * nu^2) / ((n + 1) * eta2)
}

#' Integrand.
#'
#' This is the integral from 0 to w_z wrt chi-squared pdf times standard normal pdf
#'
#' @param z normal data
#' @param n sample size
#' @param eta2 squared multiplier
#' @param nu (mu - A) / sigma
#'
#' @author David Gerard
#'
#' @noRd
fn <- function(z, n, eta2, nu) {
  stats::pchisq(q = w_z(z = z, n = n, eta2 = eta2, nu = nu), df = n - 1) * stats::dnorm(x = z)
}

#' Alpha given eta2 and nu
#'
#' Error probability given eta2 and nu
#'
#' @param n sample size
#' @param eta2 squared multiplier
#' @param nu (mu - A) / 2
#'
#' @author David Gerard
#'
#' @noRd
obj_fn <- function(n, eta2, nu) {
  stats::integrate(f = fn, lower = -Inf, upper = Inf, n = n, eta2 = eta2, nu = nu)[[1]]
}

#' Largest alpha given eta2
#'
#' Optimizes \code{obj_fn} over nu to get worst error probability given eta2.
#'
#' @param n sample size
#' @param eta2 squared multiplier
#'
#' @author David Gerard
#'
#' @noRd
worst_alpha <- function(n, eta2) {
  oout <- stats::optim(
    par = 1,
    fn = obj_fn,
    method = "L-BFGS-B",
    lower = 0,
    upper = Inf,
    control = list(fnscale = -1),
    n = n,
    eta2 = eta2)
  c(alpha = oout$value, nu = oout$par)
}

#' Find bounds on eta2 given alpha
#'
#' Just a simple boundary finding algorithm. Given alpha, we want to find
#' eta2 such that worst_alpha(n, eta2) is alpha. This finds two bounds of
#' opposite signs for worst_alpha(n, eta2) - alpha.
#'
#' @param alpha error probability
#' @param n sample size
#'
#' @author David Gerard
#'
#' @noRd
find_bounds <- function(alpha, n) {
  low <- 0.01
  high <- exp(1)
  alpha_now <- 1
  while (alpha_now > alpha) {
    aout <- worst_alpha(n = n, eta2 = high)
    alpha_now <- aout[["alpha"]]
    if (alpha_now > alpha) {
      low <- high
      high <- high * exp(1)
    }
  }
  return(c(low, high))
}

#' Returns the eta such that worst_alpha(n, eta^2) = alpha
#'
#' @param alpha error probability
#' @param n sample size
#'
#' @author David Gerard
#'
#' @noRd
eta_alpha <- function(alpha, n) {
  bounds <- find_bounds(alpha = alpha, n = n)
  eout <- stats::uniroot(f = function(eta2) alpha - worst_alpha(n = n, eta2 = eta2)[["alpha"]], interval = bounds)
  aout <- worst_alpha(n = n, eta2 = eout$root)
  c(eta = sqrt(eout$root), nu = aout[["nu"]], alpha = aout[["alpha"]])
}

#' Augmented t-interval
#'
#' Caculates a t-interval using augmented data \code{c(x,A)}. The multiplier
#' of this interval bounds the level above \code{level}, so these intervals
#' are typically conservative. This method is very fast for \code{n <= 100}
#' and \code{level %in% c(0.8, 0.9, 0.95, 0.99)} because I saved those
#' multipliers in an internal dataset. But it can be slow (on the order of
#' a second) for other confidence levels or larger sample sizes.
#'
#' @param x The vector of data
#' @param A The prior mean
#' @param level The level of the interval
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' ## A is exactly correct
#' x <- rnorm(10)
#' aug_t(x)
#' stats::t.test(x)$conf.int
#'
#' ## A is very wrong
#' x <- rnorm(10)
#' aug_t(x, A = 100)
#' stats::t.test(x)$conf.int
#'
aug_t <- function(x, A = 0, level = 0.95) {
  n <- length(x)
  alpha <- 1 - level
  TOL <- sqrt(.Machine$double.eps)
  if (n %in% augtbounds$n && any(abs(alpha - augtbounds$alpha) < TOL)) {
    eta <- augtbounds$eta[augtbounds$n == n & abs(alpha - augtbounds$alpha) < TOL]
    nu <- augtbounds$nu[augtbounds$n == n & abs(alpha - augtbounds$alpha) < TOL]
  } else {
    eout <- eta_alpha(alpha = alpha, n = n)
    eta <- eout[["eta"]]
    nu <- eout[["nu"]]
  }

  mu_hat <- mean(c(x, A))
  sig2_hat <- stats::var(c(x, A))

  c(
    mu_hat - eta * sqrt(sig2_hat) / sqrt(n + 1),
    mu_hat + eta * sqrt(sig2_hat) / sqrt(n + 1)
  )
}
