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

#' log-sum-exp trick
#'
#' Calculates log(w_1e^{x_1} + w_2e^{x_2} + ... + w_ne^{x_n})
#' for given values of x and w
#'
#' @param x The log of the summands
#' @param weights The weights. Possibly negative, as long as sum is positive,
#'     but I don't check this.
#'
#'
#' @author David Gerard
#'
#' @noRd
log_sum_exp <- function(x, weights = rep(1, length(x)), na.rm = FALSE) {
  stopifnot(length(weights) == length(x))
  if (na.rm) {
    weights <- weights[!is.na(x)]
    x <- x[!is.na(x)]
  } else if (any(is.na(x))) {
    return(NA)
  }
  if (all(x == -Inf)) {
    return(-Inf)
  } else if (any(x == Inf)) {
    return(Inf)
  } else {
    z <- max(x)
    return(log(sum(weights * exp(x-z))) + z)
  }
}


#' Hypergeometric1F1(a/2, 1/2, z) for positive integer a and positive number z
#'
#' Uses recurrence relation of https://dlmf.nist.gov/13.3#i
#' 13.3.1
#'
#' `Hypergeometric1F1[1/2, 1/2, z]` = `exp(z)`
#' `Hypergeometric1F1[3/2, 1/2, z]` = `exp(z) * (1 + 2 * z)`
#' `Hypergeometric1F1[-1, 1/2, z]` = `1 - 2 * z`
#' `Hypergeometric1F1[0, 1/2, z]` = `1`
#' `Hypergeometric1F1[1, 1/2, z]` = `1 + exp(z) * sqrt(pi) * sqrt(z) * (2 * pnorm(sqrt(2 * z)) - 1)`
#'
#' Recurrence relation:
#' `(b - a) * M(a - 1, b, z) + (2 * a - b + z) * M(a, b, z) - a * M(a + 1, b, z) = 0`
#'
#'
#' @param a an integer
#' @param z a positive number
#'
#' @author David Gerard
#'
#' @noRd
hg1f1_special <- function(a, z, log = TRUE) {
  stopifnot(z >= 0)
  if (z == 0) {
    if (log) {
      return(0)
    } else {
      return(1)
    }
  }
  if (a %% 2  == 1) {
    ## First two elements on log-scale
    m1_2 <- z
    m3_2 <- z + log(1 + 2 * z)

    if (a == 1) {
      sval <- m1_2
    } else if (a == 3) {
      sval <- m3_2
    } else {
      m_nm1 <- m1_2
      m_n <- m3_2
      n <- 1.5
      maxit <- (a - 3) / 2
      for (i in 1:maxit) {
        m_np1 <- log_sum_exp(x = c(m_nm1, m_n), weights = c(0.5 - n, 2 * n - 0.5 + z)) - log(n)
        n <- n + 1
        m_nm1 <- m_n
        m_n <- m_np1
      }
      sval <- m_np1
    }
  } else if (a %% 2 == 0) {
    ## First two elements on log-scale
    m0_2 <- 0
    m2_2 <- log_sum_exp(c(0, z + 0.5 * log(pi) + 0.5 * log(z) + log((2 * stats::pnorm(sqrt(2 * z)) - 1))))
    if (a == 0) {
      sval <- m0_2
    } else if (a == 2) {
      sval_ <- m2_2
    } else {
      m_nm1 <- m0_2
      m_n <- m2_2
      n <- 1
      maxit <- (a - 2) / 2
      for (i in 1:maxit) {
        m_np1 <- log_sum_exp(x = c(m_nm1, m_n), weights = c(0.5 - n, 2 * n - 0.5 + z)) - log(n)
        n <- n + 1
        m_nm1 <- m_n
        m_n <- m_np1
      }
      sval <- m_np1
    }
  } else {
    stop("a is not an integer")
  }
  if (!log) {
    sval <- exp(sval)
  }
  return(sval)
}


dginvnorm <- function(x, alpha = 2, mu = 0, tau = 1, log = FALSE) {
  stopifnot(alpha > 1, tau > 0)

  ## normalizing constant
  lK <- (alpha - 1) * log(tau) +
    -mu^2 / (2 * tau^2) +
    (alpha - 1) * log(2) / 2 +
    lgamma((alpha - 1) / 2) +
    hg1f1_special(a = alpha - 1, z = mu^2 / (2 * tau^2), log = TRUE)

  dval <- -1/(2 * tau^2) * (1 / x - mu)^2 - alpha * log(abs(x)) - lK

  if (!log) {
    dval <- exp(dval)
  }
  return(dval)
}



