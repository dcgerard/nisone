## Methods to calculate Kummer's confluent hypergeometric function

#' CHGM computes the confluent hypergeometric function M(a,b,x).
#'
#' This is converted to R from the original Fortran routine in
#' Zhang and Jin (1996).
#'
#' Known to not work well when a << -1 and x > 0, and when a >> 1 and x < 0
#'
#' @param a a
#' @param b b
#' @param x x
#'
#' @section Licensing:
#' This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
#' they give permission to incorporate this routine into a user program
#' provided that the copyright is acknowledged.
#'
#' @references
#' Shanjie Zhang, Jianming Jin,
#' Computation of Special Functions,
#' Wiley, 1996,
#' ISBN: 0-471-11963-6,
#' LC: QA351.C45.
#'
#' Original Fortran code available at: \url{https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.f}
#'
#' @author Shanjie Zhang, Jianming Jin, David Gerard (just converted it to R)
#'
#' @noRd
chgm <- function(a, b, x) {
  a0 <- a
  a1 <- a
  x0 <- x
  hg <- 0.0

  if (b == 0.0 || b == -abs(as.integer(b))) {
    hg <- 1.0e300
  } else if (a == 0.0 || x == 0.0) {
    hg <- 1.0
  } else if (a == -1.0) {
    hg <- 1.0 - x / b
  } else if (a == b) {
    hg <- exp(x)
  } else if (a - b == 1.0) {
    hg <- (1.0 + x / b) * exp(x)
  } else if (a == 1.0 && b == 2.0) {
    hg <- (exp(x) - 1.0) / x
  } else if (a == as.integer(a) && a < 0.0) {
    m <- as.integer(-a)
    r <- 1.0
    hg <- 1.0
    for (k in 1:m) {
      r <- r * (a + k - 1.0) / k / (b + k - 1.0) * x
      hg <- hg + r
    }
  }

  if (hg != 0.0) {
    return(hg)
  }

  if (x < 0.0) {
    a <- b - a
    a0 <- a
    x <- abs(x)
  }

  nl <- if (a < 2.0) 0 else 1

  if (a >= 2.0) {
    la <- as.integer(a)
    a <- a - la - 1.0
  }

  for (n in 0:nl) {
    if (a0 >= 2.0) {
      a <- a + 1.0
    }

    if (x <= 30.0 + abs(b) || a < 0.0) {
      hg <- 1.0
      rg <- 1.0
      for (j in 1:500) {
        rg <- rg * (a + j - 1.0) / (j * (b + j - 1.0)) * x
        hg <- hg + rg
        if (abs(rg / hg) < 1.0e-15) {
          break
        }
      }
    } else {
      ta <- gamma(a)
      tb <- gamma(b)
      xg <- b - a
      tba <- gamma(xg)
      sum1 <- 1.0
      sum2 <- 1.0
      r1 <- 1.0
      r2 <- 1.0
      for (i in 1:8) {
        r1 <- -r1 * (a + i - 1.0) * (a - b + i) / (x * i)
        r2 <- -r2 * (b - a + i - 1.0) * (a - i) / (x * i)
        sum1 <- sum1 + r1
        sum2 <- sum2 + r2
      }
      hg1 <- tb / tba * x^(-a) * cos(pi * a) * sum1
      hg2 <- tb / ta * exp(x) * x^(a - b) * sum2
      hg <- hg1 + hg2
    }

    if (n == 0) {
      y0 <- hg
    } else if (n == 1) {
      y1 <- hg
    }
  }

  if (a0 >= 2.0) {
    for (i in 1:(la - 1)) {
      hg <- ((2.0 * a - b + x) * y1 + (b - a) * y0) / a
      y0 <- y1
      y1 <- hg
      a <- a + 1.0
    }
  }

  if (x0 < 0.0) {
    hg <- hg * exp(x0)
  }

  a <- a1
  x <- x0

  return(hg)
}

#' Different ways to calculate Kummer's confluent hypergeometric functions.
#'
#' I didn't end up using this since hg1f1_special() worked better in my use case.
#' You can also try the gsl::hyperg_1F1() function for good general use, but
#' bad use case in my application.
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

#' Hypergeometric1F1(a/2, 1/2, z) for positive integer a and positive number z
#'
#' AKA M(a/2, 1/2, z)
#'
#' Uses recurrence relation of https://dlmf.nist.gov/13.3#i
#' 13.3.1
#'
#' `Hypergeometric1F1[1/2, 1/2, z]` = `exp(z)`
#' `Hypergeometric1F1[3/2, 1/2, z]` = `exp(z) * (1 + 2 * z)`1z
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
      sval <- m2_2
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

#' Kummer's confluent hypergeometric function 1F1(a; b; z)
#'
#' Evaluates the series  sum_{n=0}^{Inf} (a)_n / (b)_n * z^n / n!
#' which converges for all finite z.
#' <https://dlmf.nist.gov/13.2#E2>
#'
#' Terms are tracked as (log|term|, sign), accumulated via log-sum-exp
#' into separate positive and negative buckets, combined at the end.
#' This avoids intermediate overflow when z is large.
#'
#' @param log logical; if TRUE return log(1F1) instead of 1F1.
#'
#' @noRd
.hyp1f1 <- function(a, b, z, log = FALSE, tol = 1e-14, max_iter = 1000L) {

  # Log-scale summation via log-sum-exp.
  # Terms are tracked as (log|term|, sign) and accumulated into separate
  # positive and negative buckets, so intermediate overflow is impossible.
  # In our application z = mu^2/(2*tau^2) >= 0, so all terms are positive
  # and log_neg stays at -Inf throughout; the two-bucket logic is retained
  # for correctness if z < 0 is ever passed.
  log_term  <- 0     # log|term_n|, initialised to log(1) = 0
  sign_term <- 1     # sign of term_n

  log_pos <- 0       # log of running sum of positive terms  (starts at 1)
  log_neg <- -Inf    # log of |running sum of negative terms| (none yet)

  for (n in seq_len(max_iter)) {
    ratio     <- (a + n - 1L) / ((b + n - 1L) * n) * z
    log_term  <- log_term + base::log(abs(ratio))
    sign_term <- sign_term * sign(ratio)

    if (sign_term > 0) {
      log_pos <- log_pos + log1p(exp(log_term - log_pos))
    } else {
      log_neg <- if (is.infinite(log_neg))
        log_term
      else
        log_neg + log1p(exp(log_term - log_neg))
    }

    # Convergence: |term_n| < tol * |partial sum|, evaluated in log space
    if (log_pos >= log_neg) {
      log_abs_s <- log_pos + log1p(-exp(log_neg - log_pos))
    } else {
      log_abs_s <- log_neg + log1p(-exp(log_pos - log_neg))
    }

    if (log_term < base::log(tol) + log_abs_s) {
      if (log_pos >= log_neg) {
        result_sign <- 1
      } else {
        result_sign <- -1
      }
      if (result_sign < 0 && log) {
        stop(".hyp1f1: log requested but 1F1 is negative")
      }
      if (log) {
        return(log_abs_s)
      } else {
        return(result_sign * exp(log_abs_s))
      }
    }
  }

  ## If it got here, it failed to converge, so use hg1f1_special() as a backup
  if (b == 0.5 && z >= 0) {
    return(hg1f1_special(a = 2 * a, z = z, log = log))
  }

  ## If my special conditions are not met, throw an error.
  stop(sprintf(
    ".hyp1f1: failed to converge after %d iterations (a=%g, b=%g, z=%g)",
    max_iter, a, b, z))
}

#' Log of normalizing constant K(alpha, mu, tau)
#'
#' @noRd
.gin_log_K <- function(alpha, mu, tau) {
  a     <- alpha
  log_h <- .hyp1f1(0.5 * (a - 1), 0.5, mu^2 / (2 * tau^2), log = TRUE)
  log_inv_K <- (a - 1) * log(tau) +
    (-mu^2 / (2 * tau^2)) +
    ((a - 1) / 2) * log(2) +
    lgamma(0.5 * (a - 1)) +
    log_h
  return(-log_inv_K)
}
