#' Different ways to calculate Kummer's confluent hypergeometric functions.
#'
#' This is often denoted by M(a,b,z) or 1F1(a,b,z). I use the NIST
#' computational strategies \url{https://dlmf.nist.gov/13.29}.
#'
#' @param a first value
#' @param b Second value
#' @param z Third value
#' @param method Either Maclaurin series or Integral solution
#' @param nterms The number of terms of the Macllaurin series
#'
#' @author David Gerard
#'
#' @noRd
hypergeo1f1 <- function(a, b, z, method = c("series", "integral"), nterms = 10) {
  method <- match.arg(method)
  if (method == "integral") {
    stopifnot(b > a)
    f <- function(t) exp(z * t) * t ^ (a - 1) * (1 - t) ^ (b - a - 1)
    intout <- stats::integrate(f, lower = 0, upper = 1)
    val <- gamma(b) * intout$value / (gamma(a) * gamma(b - a))
  } else if (method == "series") {
    val <- 1
    for (i in seq_len(nterms)) {
      val <- val + prod(seq(a, a + i - 1, by = 1)) / (factorial(i) * prod(seq(b, b + i - 1, by = 1))) * z ^ i
    }
  } else if (method == "recurrance") {
    # https://dlmf.nist.gov/13.3#i
    stopifnot(a %% 1 == 0)
    m0 <- 1
    m1 <- (b - 1) * exp(z) * z ^ (1 - b) * (gamma(b - 1)) ## not done, need incomplete gamma
  }
  return(val)
}
