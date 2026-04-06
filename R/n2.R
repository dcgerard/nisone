#' w as a function of z
#'
#' @param z standard normal value
#' @param n sample size
#' @param eta2 multiplier squared
#' @param nu (mu - A) / sigma
#' @param wt Prior weight for A. This is interpreted as how many replicates of
#'     A the t-statistic has.
#'
#' @author David Gerard
#'
#' @noRd
w_z <- function(z, n, eta2, nu, wt = 1) {
  ## ((n * z - sqrt(n) * nu)^2 - eta2 * (z + sqrt(n) * nu)^2) / ((n + 1) * eta2)
  ## ((n^2 - eta2) * z^2 - 2 * (n^(1.5) + sqrt(n) * eta2) * nu * z + (1 - eta2) * n * nu^2) / ((n + 1) * eta2) ## old way
  wt^2 * (n + wt - 1) * (n * z / wt - sqrt(n) * nu)^2 / ((n + wt) * n * eta2) - wt * (z + sqrt(n) * nu)^2 / (n + wt)
}

#' Integrand.
#'
#' This is the integral from 0 to w_z wrt chi-squared cdf times standard normal pdf
#'
#' @param z normal data
#' @param n sample size
#' @param eta2 squared multiplier
#' @param nu (mu - A) / sigma
#' @param wt Prior weight for A. This is interpreted as how many replicates of
#'     A the t-statistic has.
#'
#' @author David Gerard
#'
#' @noRd
fn <- function(z, n, eta2, nu, wt = 1) {
  stats::pchisq(q = w_z(z = z, n = n, eta2 = eta2, nu = nu, wt = wt), df = n - 1) * stats::dnorm(x = z)
}

#' Alpha given eta2 and nu
#'
#' Error probability given eta2 and nu
#'
#' @param n sample size
#' @param eta2 squared multiplier
#' @param nu (mu - A) / 2
#' @param wt Prior weight for A. This is interpreted as how many replicates of
#'     A the t-statistic has.
#'
#' @author David Gerard
#'
#' @noRd
obj_fn <- function(n, eta2, nu, wt = 1) {
  stats::integrate(f = fn, lower = -Inf, upper = Inf, n = n, eta2 = eta2, nu = nu, wt = wt)[[1]]
}

#' Largest alpha given eta2
#'
#' Optimizes \code{obj_fn} over nu to get worst error probability given eta2.
#'
#' @param n sample size
#' @param eta2 squared multiplier
#' @param wt Prior weight for A. This is interpreted as how many replicates of
#'     A the t-statistic has.
#'
#' @author David Gerard
#'
#' @noRd
worst_alpha <- function(n, eta2, wt = 1) {
  oout <- stats::optim(
    par = 1,
    fn = obj_fn,
    method = "L-BFGS-B",
    lower = 0,
    upper = Inf,
    control = list(fnscale = -1),
    n = n,
    eta2 = eta2,
    wt = wt)
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
#' @param wt Prior weight for A. This is interpreted as how many replicates of
#'     A the t-statistic has.
#'
#' @author David Gerard
#'
#' @noRd
find_bounds <- function(alpha, n, wt = 1) {
  low <- 0.01
  high <- exp(1)
  alpha_now <- 1
  while (alpha_now > alpha) {
    aout <- worst_alpha(n = n, eta2 = high, wt = wt)
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
eta_alpha <- function(alpha, n, wt = 1) {
  bounds <- find_bounds(alpha = alpha, n = n, wt = wt)
  eout <- stats::uniroot(f = function(eta2) alpha - worst_alpha(n = n, eta2 = eta2, wt = wt)[["alpha"]], interval = bounds)
  aout <- worst_alpha(n = n, eta2 = eout$root, wt = wt)
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
#' @param wt Prior weight for A. This is interpreted as how many replicates of
#'     A the t-statistic has.
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
aug_t <- function(x, A = 0, level = 0.95, wt = 1) {
  n <- length(x)
  alpha <- 1 - level
  TOL <- sqrt(.Machine$double.eps)
  if (n %in% augtbounds$n && any(abs(alpha - augtbounds$alpha) < TOL) && wt == 1) {
    eta <- augtbounds$eta[augtbounds$n == n & abs(alpha - augtbounds$alpha) < TOL]
    nu <- augtbounds$nu[augtbounds$n == n & abs(alpha - augtbounds$alpha) < TOL]
  } else {
    eout <- eta_alpha(alpha = alpha, n = n, wt = wt)
    eta <- eout[["eta"]]
    nu <- eout[["nu"]]
  }

  mu_hat <- mean_aug(x = x, A = A, wt = wt)
  sig2_hat <- var_aug(x = x, A = A, wt = wt)

  c(
    mu_hat - eta * sqrt(sig2_hat) / sqrt(n + wt),
    mu_hat + eta * sqrt(sig2_hat) / sqrt(n + wt)
  )
}

#' Augmented mean
#'
#' If wt is an integer, this is mean(c(x, rep(A, wt))). But we generalize it to
#' non-integer wt as well.
#'
#' @param x Data
#' @param A Prior value
#' @param wt Weight on prior value.
#'
#' @author David Gerard
#'
#' @export
#'
#' @return The augmented mean
#'
#' @examples
#' mean_aug(c(1, 2, 3), 4, wt = 10)
#'
mean_aug <- function(x, A, wt = 1) {
  n <- length(x)
  n / (n + wt) * mean(x) + wt / (n + wt) * A
}

#' Augmented variance
#'
#' If wt is an integer, this is var(c(x, rep(A, wt))). But we generalize it to
#' non-integer wt as well.
#'
#' @param x Data
#' @param A Prior value
#' @param wt Weight on prior value.
#'
#' @author David Gerard
#'
#' @export
#'
#' @return The augmented variance
#'
#' @examples
#' var_aug(c(1, 2, 3), 4, wt = 10)
#'
var_aug <- function(x, A, wt = 1) {
  n <- length(x)
  (n - 1) / (n + wt - 1) * stats::var(x) + n * wt * (mean(x) - A)^2 / ((n + wt - 1) * (n + wt))
}


