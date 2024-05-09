test_that("invnorm works", {
  pvec <- seq(0.1, 0.9, by = 0.1)
  qvec <- qinvnorm(p = pvec, imean = 0, isd = 1)
  pvec2 <- pinvnorm(x = qvec, imean = 0, isd = 1)
  expect_equal(pvec, pvec2)

  qvec <- seq(-5, 5, by = 0.1)
  pvec <- pinvnorm(x = qvec, imean = 1, isd = 2)
  qvec2 <- qinvnorm(p = pvec, imean = 1, isd = 2)
  expect_equal(qvec, qvec2)

  x <- -3
  p1 <- stats::integrate(f = dinvnorm, lower = -Inf, upper = x)
  p2 <- pinvnorm(x = x)
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
