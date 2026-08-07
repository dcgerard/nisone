#' Posterior sample using inverse mean prior
#'
#' Let \eqn{X_1,\ldots,X_n \sim N(\mu,\mu^2/\nu^2)}. This function performs
#' Bayesian inference assuming prior \eqn{\pi(\mu) = |\mu - A|^{-1}} and
#' \eqn{\nu \sim \text{Gamma}(a,b)}. You can choose \eqn{a} and \eqn{b}.
#'
#' If \eqn{a = b \rightarrow \infty}, this gives you the same posterior as
#' \code{\link{bcin}()}. If \eqn{a = b = 1/2}, this gives you
#' \eqn{\frac{\mu - \hat{\mu}}{\hat{\sigma}/\sqrt{n+1}} | X \sim t_{n}}, where
#' \eqn{\hat{\mu}} and \eqn{\hat{\sigma}} are the mean and standard deviation,
#' respectively, of \code{c(x, A)}.
#'
#' @param x The data. A numeric vector.
#' @param A The augmented data.
#' @param iter Number of MCMC iterations.
#' @param warmup The burnin number of iterations.
#' @param nu2_shape Shape hyperparameter of Gamma distribution of nu2.
#' @param nu2_rate Rate hyperparameter of Gamma distribution of nu2.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- stats::rnorm(n = 2, mean = 5, sd = 3)
#' A <- 0
#'
#' ## For large nu2_rate and nu2_shape, get bcin ----
#' pd <- post_samp(x = x, A = A, nu2_shape = 100, nu2_rate = 100)
#' quantile(pd$mu, c(0.025, 0.975))
#' bcin(x = x, A = A)
#'
#' alpha <- length(x) + 1
#' mu <- sum(x) / sum(x^2)
#' tau <- 1 / sqrt(sum(x^2))
#' qvec <- qginvnorm(p = ppoints(length(pd$mu)), alpha = alpha, mu = mu, tau = tau)
#' stats::qqplot(x = pd$mu, y = qvec, xlim = c(-100, 100))
#' graphics::abline(0, 1, lty = 2, col = 2)
#'
#' ## For nu2_shape = nu2_rate = 0.5, get t_n ----
#' pd <- post_samp(x = x, A = A, nu2_shape = 0.5, nu2_rate = 0.5)
#' quantile(pd$mu, c(0.025, 0.975))
#' x_aug <- c(x, A)
#' c(
#'   mean(x_aug) - qt(0.975, df = length(x)) * sd(x_aug) / sqrt(length(x) + 1),
#'   mean(x_aug) + qt(0.975, df = length(x)) * sd(x_aug) / sqrt(length(x) + 1)
#' )
#'
#' tvec <- (pd$mu - mean(x_aug)) / (sd(x_aug) / sqrt(length(x_aug)))
#' qvec <- qt(ppoints(length(tvec)), df = length(x))
#' stats::qqplot(x = tvec, y = qvec, xlim = c(-50, 50))
#' graphics::abline(0, 1, lty = 2, col = 2)
#'
#' }
#'
#' @author David Gerard
#'
#' @noRd
post_samp <- function(x, A = 0, iter = 10000, warmup = floor(iter/2), nu2_shape = 0.5, nu2_rate = 0.5) {
  y <- x - A
  n <- length(y)

  ## statistics
  yrat <- sum(y) / sum(y^2)
  y2 <- sqrt(sum(y^2))

  ## parameter samples
  beta_vec <- rep(NA_real_, length.out = iter)
  nu2_vec <- rep(NA_real_, length.out = iter)

  ## initial values
  beta <- mean(y)
  nu2 <- sum(y^2) / sum(y)^2

  for (i in seq_len(iter)) {

    if (i %in% round(iter * seq(0, 1, by = 0.1))) {
      message(paste0(round(i / iter * 100), " percent done"))
    }

    ## nu given beta
    nu2 <- stats::rgamma(n = 1, shape = n / 2 + nu2_shape, rate = sum((y / beta - 1)^2) / 2 + nu2_rate)

    ## beta given nu
    beta <- rginvnorm(n = 1, alpha = n + 1, mu = yrat, tau = 1 / (sqrt(nu2) * y2))

    nu2_vec[[i]] <- nu2
    beta_vec[[i]] <- beta
  }

  mu_vec <- beta_vec + A

  ret <- data.frame(mu = mu_vec[(warmup+1):iter], nu2 = nu2_vec[(warmup+1):iter])

  return(ret)
}
