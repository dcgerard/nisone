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
  ival <- integrate(f = f, lower = -Inf, upper = Inf)
  expect_equal(ival$value, 1)
})
