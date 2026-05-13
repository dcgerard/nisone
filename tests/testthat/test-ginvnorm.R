test_that("1F1 works", {
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 1, method = "integral"), exp(1) - 1)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "integral"), (exp(2) - 1) / 2)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "integral"), (exp(2) - 1) / 2)
  expect_equal(hypergeo1f1(a = 2, b = 3, z = 1, method = "integral"), 2)
  expect_equal(hypergeo1f1(a = 3, b = 5, z = 3, method = "integral"), (exp(3) - 4) * 4 / 9)

  expect_equal(hypergeo1f1(a = 1, b = 2, z = 1, method = "series"), exp(1) - 1, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "series"), (exp(2) - 1) / 2, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "series"), (exp(2) - 1) / 2, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 3, b = 5, z = 3, method = "series"), (exp(3) - 4) * 4 / 9, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 4, b = 3, z = 1, method = "series"), exp(1) * 4 / 3, tolerance = 1e-3)


  expect_equal(
    exp(hg1f1_special(3, 2)),
    hypergeo1f1(a = 3/2, b = 1/2, z = 2, method = "series", nterms = 100)
  )

  expect_equal(
    exp(hg1f1_special(4, 2)),
    hypergeo1f1(a = 4/2, b = 1/2, z = 2, method = "series", nterms = 100)
  )

  expect_equal(
    exp(hg1f1_special(11, 2)),
    hypergeo1f1(a = 11/2, b = 1/2, z = 2, method = "series", nterms = 100)
  )

  expect_equal(
    exp(hg1f1_special(12, 2)),
    hypergeo1f1(a = 12/2, b = 1/2, z = 2, method = "series", nterms = 100)
  )

  expect_equal(
    exp(hg1f1_special(12, 0)),
    hypergeo1f1(a = 0, b = 1/2, z = 2, method = "series", nterms = 100)
  )

})

test_that("lse works", {
  x <- c(1, 3, -1)
  expect_equal(
    log(sum(exp(x))),
    log_sum_exp(x)
  )
  expect_equal(log_sum_exp(c(-Inf, 1)), 1)
  expect_equal(log_sum_exp(c(-Inf, -Inf)), -Inf)
  expect_equal(log_sum_exp(c(1, Inf)), Inf)
  expect_equal(log_sum_exp(c(1, 2, NA), na.rm = FALSE), NA)
  expect_equal(log_sum_exp(c(1, 2, NA), na.rm = TRUE), log(sum(exp(c(1, 2)))))
})


test_that("dginvnorm works", {
  alpha <- 10
  mu <- -10
  tau <- 1
  f <- function(x) dginvnorm(x = x, alpha = alpha, mu = mu, tau = tau, log = FALSE)
  ival <- stats::integrate(f = f, lower = -Inf, upper = Inf)
  expect_equal(ival$value, 1)

  x <- seq(-10, 10, length.out = 500)
  y <- dginvnorm(x = x, alpha = 2, mu = 1/10, tau = 1/10)
  graphics::plot(x, y, type = "l")

  # bench::mark(
  #   ginormal::dgin(z = 100, alpha = 100, mu = 5, tau = 1, log = FALSE),
  #   dginvnorm(x = 100, alpha = 100, mu = 5, tau = 1, log = FALSE)
  # )
})

test_that("pginvnorm and qginvnorm are inverses", {
  set.seed(1)
  pvec <- stats::runif(10)

  mu <- rnorm(1)
  tau <- rchisq(1, df = 1)
  alpha <- 4
  expect_equal(
    pginvnorm(qginvnorm(p = pvec, alpha = alpha, mu = mu, tau = tau), alpha = alpha, mu = mu, tau = tau),
    pvec,
    tolerance = 1e-4
  )

  expect_equal(
    pginvnorm(q = 1000, alpha = 10, mu = 0, tau = 1),
    1
  )
  expect_equal(
    pginvnorm(q = -1000, alpha = 10, mu = 0, tau = 1),
    0
  )
})

test_that("ginvnorm is the same as invnorm for alpha = 2", {
  expect_equal(
    pinvnorm(q = 1, imean = 1, isd = 2),
    pginvnorm(q = 1, alpha = 2, mu= 1, tau = 2)
  )

  expect_equal(
    qinvnorm(p = 0.975, imean = 2, isd = 4),
    qginvnorm(p = 0.975, alpha = 2, mu= 2, tau = 4),
    tolerance = 1e-4
  )

  expect_equal(
    qinvnorm(p = 0.975, imean = 10, isd = 6),
    qginvnorm(p = 0.975, alpha = 2, mu= 10, tau = 6),
    tolerance = 1e-4
  )

  expect_equal(
    qinvnorm(p = 0.01, imean = 10, isd = 6),
    qginvnorm(p = 0.01, alpha = 2, mu= 10, tau = 6),
    tolerance = 1e-4
  )

})


test_that("pginvnorm_old() and pginvnorm() are the same", {
  skip("Not ready yet")
  # I use modal information in pginvnorm() that should make it better behaved for more extreme tau

  expect_equal(
    pginvnorm(q = 0.001, alpha = 3, mu = 10, tau = 100),
    pginvnorm_old(q = 0.001, alpha = 3, mu = 10, tau = 100)
  )

  alpha <- sample(2:100, size = 1)
  mu <- rnorm(n = 1, sd = 100)
  tau <- rchisq(n = 1, df = 1)
  modes <- xginvnorm(alpha = alpha, mu = mu, tau = tau)
  # expect_equal(
  #   pginvnorm(q = modes[[2]], alpha = alpha, mu = mu, tau = tau),
  #   pginvnorm_old(q = modes[[2]], alpha = alpha, mu = mu, tau = tau)
  # )

  cred <- qginvnorm(p = c(0.025, 0.975), alpha = alpha, mu = mu, tau = tau)
  pginvnorm(q = cred[[1]], alpha = alpha, mu = mu, tau = tau)
  pginvnorm(q = cred[[2]], alpha = alpha, mu = mu, tau = tau)

  alpha <- 50
  mu <- -20
  tau <- 1e-3
  modes <- xginvnorm(alpha = alpha, mu = mu, tau = tau)
  expect_equal(
    pginvnorm(q = modes[[2]], alpha = alpha, mu = mu, tau = tau),
    1,
    tolerance = 1e-6
  )

})
