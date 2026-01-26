## Bayes Factor calculations from t statistics

#' Bayes factor from t-statistic
#'
#' Given a user-provided prior on the standardized effect size, will calculate
#' the Bayes factor.
#'
#' Let \eqn{\delta} be the standardized effect size, let \eqn{\sigma^2} be the variance
#' (assumed equal in two-sample case). For a two-sample t-test, we place
#' the prior \eqn{1/\sigma^2} under the null and \eqn{\pi(\delta)/\sigma^2}
#' under the alternative, for some arbitrary density \eqn{\pi(\cdot)}. Given
#' this setting, the Bayes factor is a function of the t-statistic. We calculate
#' it via numeric integration.
#'
#' @param t The observed t-statistic.
#' @param nu The degrees of freedom of the t-statistic.
#' @param prior The prior on the standardized effect size. In the one-sample case
#'     this is the prior over \eqn{\frac{\mu - \mu_0}{\sigma}}, where \eqn{\mu_0}
#'     is the null value. In the two-sample case this is the prior over
#'     \eqn{\frac{\mu_1 - \mu_2}{\sigma}}, where \eqn{\mu_1} and \eqn{\mu_2} are
#'     the means of the two-samples. The default prior is a Cauchy, but you
#'     can put in any function you want.
#'
#' @return The Bayes factor to a corresponding t-statistic.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Gronau, Q. F., Ly, A., & Wagenmakers, E. J. (2020). Informed Bayesian t-tests. \emph{The American Statistician}. \doi{10.1080/00031305.2018.1562983}}
#' }
#'
#' @export
bft <- function(t, nu, prior = NULL) {
  if (is.null(prior)) {
    prior <- stats::dcauchy
  }

  f <- function(delta) {
    stats::dt(x = t, df = nu, ncp = delta * sqrt(nu + 1)) * prior(delta) /
      stats::dt(x = t, df = nu)
  }

  ## I don't think this is an issue because it only shows up in the
  ## extreme tails, where the area is negligible
  suppressWarnings(
    bfout <- stats::integrate(f = f, lower = -Inf, upper = Inf)
  )

  return(bfout[[1]])
}
