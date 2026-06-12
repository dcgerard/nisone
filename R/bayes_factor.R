## Bayes Factor calculations from t statistics

#' Bayes factor from t-statistic
#'
#' Given a user-provided prior on the standardized effect size, will calculate
#' the Bayes factor.
#'
#' Let \eqn{\delta} be the standardized effect size, let \eqn{\sigma^2} be
#' the variance (assumed equal in two-sample case). In the one-sample case,
#' \eqn{\delta = \frac{\mu - \mu_0}{\sigma}}, where \eqn{\mu_0} is the null
#' value. In the two-sample case, \eqn{\delta = \frac{\mu_1 - \mu_2}{\sigma}},
#' where \eqn{\mu_1} and \eqn{\mu_2} are the means of the two-samples. We place
#' the prior \eqn{1/\sigma^2} under the null and \eqn{\pi(\delta)/\sigma^2}
#' under the alternative, for some arbitrary density \eqn{\pi(\cdot)}. Given
#' this setting, the Bayes factor is a function of the \eqn{t}-statistic.
#' We calculate it via numeric integration.
#'
#' @param t The observed t-statistic.
#' @param nu The degrees of freedom of the t-statistic. This should be
#'    \eqn{n-1} in the one-sample case. This should be \eqn{n_1+n_2-2} in the two-sample case.
#' @param nd Effective sample size. This should be \eqn{n} in the one-sample case.
#'    This should be \eqn{(\frac{1}{n_1} + \frac{1}{n_2})^{-1}} in the two-sample case.
#' @param prior The prior on the standardized effect size, \eqn{\delta}. In the one-sample case
#'     this is the prior over \eqn{\delta = \frac{\mu - \mu_0}{\sigma}}, where \eqn{\mu_0}
#'     is the null value. In the two-sample case this is the prior over
#'     \eqn{\delta = \frac{\mu_1 - \mu_2}{\sigma}}, where \eqn{\mu_1} and \eqn{\mu_2} are
#'     the means of the two-samples. The default prior is a Cauchy, but you
#'     can put in any density function you want. If your prior has pointmasses
#'     in it, this function won't work.
#'
#' @returns The Bayes factor to a corresponding t-statistic.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Gronau, Q. F., Ly, A., & Wagenmakers, E. J. (2020). Informed Bayesian t-tests. \emph{The American Statistician}. \doi{10.1080/00031305.2018.1562983}}
#' }
#'
#' @examples
#' # One sample t, n = 10, t-statistic = 2
#' bft(t = 2, nu = 10 - 1, nd = 10)
#'
#' # Two sample t, n1 = 10, n2 = 8, t-statistic = 2
#' bft(t = 2, nu = 10 + 8 - 2, nd = 1 / (1 / 10 + 1 / 8))
#'
#' @export
bft <- function(t, nu, nd, prior = NULL) {
  if (is.null(prior)) {
    prior <- stats::dcauchy
  }

  f <- function(delta) {
    stats::dt(x = t, df = nu, ncp = delta * sqrt(nd)) * prior(delta) /
      stats::dt(x = t, df = nu)
  }

  ## I don't think this is an issue because it only shows up in the
  ## extreme tails, where the area is negligible
  suppressWarnings(
    bfout <- stats::integrate(f = f, lower = -Inf, upper = Inf)
  )

  return(bfout[[1]])
}
