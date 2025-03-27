## W as a function of Z
w_z <- function(z, n, eta, nu) {
  ## (n * z + nu)^2 / ((n + 1) * eta) - eta * (z - nu)^2 / ((n + 1) * eta)
  ((n^2 - eta) * z^2 + 2 * (n + eta) * nu * z + (1 - eta) * nu^2) / ((n + 1) * eta)
}

## integrand
fn <- function(z, n, eta, nu) {
  stats::pchisq(q = w_z(z = z, n = n, eta = eta, nu = nu), df = n - 1) * stats::dnorm(x = z)
}

## Probability, objective
obj_fn <- function(n, eta, nu) {
  stats::integrate(f = fn, lower = -Inf, upper = Inf, n = n, eta = eta, nu = nu)[[1]]
}

## worst case alpha, optimized over nu
worst_alpha <- function(n, eta) {
  oout <- stats::optim(
    par = 1,
    fn = obj_fn,
    method = "L-BFGS-B",
    lower = 0,
    upper = Inf,
    control = list(fnscale = -1),
    n = n,
    eta = eta)
  c(alpha = oout$value, nu = oout$par)
}

## find bounds of eta
find_bounds <- function(alpha, n) {
  low <- 0
  high <- exp(1)
  alpha_now <- 1
  while (alpha_now > alpha) {
    aout <- worst_alpha(n = n, eta = high)
    alpha_now <- aout[["alpha"]]
    if (alpha_now > alpha) {
      low <- high
      high <- high * exp(1)
    }
  }
  return(c(low, high))
}

eta_alpha <- function(alpha, n) {
  bounds <- find_bounds(alpha = alpha, n = n)
  eout <- stats::uniroot(f = \(eta) alpha - worst_alpha(n = n, eta = eta)[["alpha"]], interval = bounds)
  aout <- worst_alpha(n = n, eta = eout$root)
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
