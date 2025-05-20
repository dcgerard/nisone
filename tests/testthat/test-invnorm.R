test_that("invnorm works", {
  pvec <- seq(0.1, 0.9, by = 0.1)
  qvec <- qinvnorm(p = pvec, imean = 0, isd = 1)
  pvec2 <- pinvnorm(q = qvec, imean = 0, isd = 1)
  expect_equal(pvec, pvec2)

  qvec <- seq(-5, 5, by = 0.1)
  pvec <- pinvnorm(q = qvec, imean = 1, isd = 2)
  qvec2 <- qinvnorm(p = pvec, imean = 1, isd = 2)
  expect_equal(qvec, qvec2)

  x <- -3
  p1 <- stats::integrate(f = dinvnorm, lower = -Inf, upper = x)
  p2 <- pinvnorm(q = x)
  expect_equal(p1[[1]], p2[[1]])
})


test_that("median correct", {
  x <- 2
  med1 <- x / (1 - stats::qnorm(-0.5 + stats::pnorm(1)))
  med2 <- qinvnorm(p = 0.5, imean = 1/x, isd = 1/abs(x))
  expect_equal(med1, med2)

  x <- -3
  med1 <- x / (1 - stats::qnorm(-0.5 + stats::pnorm(1)))
  med2 <- qinvnorm(p = 0.5, imean = 1/x, isd = 1/abs(x))
  expect_equal(med1, med2)
})

test_that("invnorm and ginvnorm for alpha = 2 are the same", {
  expect_equal(
    dinvnorm(x = 1, imean = 2, isd = 3, log = TRUE),
    dginvnorm(x = 1, mu = 2, tau = 3, alpha = 2, log = TRUE)
  )

  expect_equal(
    pinvnorm(q = 1, imean = 2, isd = 3),
    pginvnorm(q = 1, mu = 2, tau = 3, alpha = 2),
    tolerance = 1e-6
  )

  expect_equal(
    dinvnorm(x = 10, imean = -50, isd = 20, log = TRUE),
    dginvnorm(x = 10, mu = -50, tau = 20, alpha = 2, log = TRUE),
    tolerance = 1e-6
  )
})


test_that("mean is same as in Robert (1991) Generalized inverse normal distributions", {
  # See Lemma 2 of that paper
  alpha <- 6
  mu <- 3
  tau <- 2

  rmean <- mginvnorm(alpha = alpha, mu = mu, tau = tau)

  f <- function(z) {
    z * dginvnorm(x = z, alpha = alpha, mu = mu, tau = tau, log = FALSE)
  }
  nmean <- stats::integrate(f = f, lower = -Inf, upper = Inf)[[1]]

  expect_equal(rmean, nmean)
})
