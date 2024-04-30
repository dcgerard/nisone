test_that("get same alpha as Blachman and Machol", {
  expect_equal(wc_alpha(width = 1, center = "X"), 0.5)
  expect_equal(wc_alpha(width = 2.42, center = "X"), 0.2, tolerance = 1e-3)
  expect_equal(wc_alpha(width = 4.84, center = "X"), 0.1, tolerance = 1e-3)
  expect_equal(wc_alpha(width = 9.68, center = "X"), 0.05, tolerance = 1e-3)
  expect_equal(wc_alpha(width = 48.39, center = "X"), 0.01, tolerance = 1e-3)

  expect_equal(wc_alpha(width = 0.5, center = "X/2"), 0.5)
  expect_equal(wc_alpha(width = 2.31, center = "X/2"), 0.2, tolerance = 1e-3)
  expect_equal(wc_alpha(width = 4.79, center = "X/2"), 0.1, tolerance = 1e-3)
  expect_equal(wc_alpha(width = 9.65, center = "X/2"), 0.05, tolerance = 1e-3)
  expect_equal(wc_alpha(width = 48.39, center = "X/2"), 0.01, tolerance = 1e-3)
})


test_that("Simulations at actual confidence level", {
  set.seed(1)
  # X center
  n <- 100000
  mu <- 100
  width <- 2.42
  level <- 1 - wc_alpha(width = width, center = "X")
  lambda <- wc_cov(width = width, center = "X")
  sigma <- mu / lambda
  x <- stats::rnorm(n = n, mean = mu, sd = sigma)
  lower <- x - width * abs(x)
  upper <- x + width * abs(x)
  expect_equal(mean(lower <= mu & mu <= upper), level, tolerance = 1e-2)

  # X/2 center
  n <- 100000
  mu <- 100
  width <- 2.31
  level <- 1 - wc_alpha(width = width, center = "X/2")
  lambda <- wc_cov(width = width, center = "X/2")
  sigma <- mu / lambda
  x <- stats::rnorm(n = n, mean = mu, sd = sigma)
  lower <- x/2 - width * abs(x)
  upper <- x/2 + width * abs(x)
  expect_equal(mean(lower <= mu & mu <= upper), level, tolerance = 1e-2)
})

test_that("cdf and density are part of same family", {
  x <- -3
  val1 <- stats::integrate(d_wc, lower = -Inf, upper = x)
  val2 <- p_wc(q = x)
  expect_equal(val1[[1]], val2)

  x <- 1.8
  val1 <- stats::integrate(d_wc, lower = x, upper = Inf)
  val2 <- 1 - p_wc(q = x)
  expect_equal(val1[[1]], val2)

  x <- -3
  val1 <- stats::integrate(d_wc, lower = -Inf, upper = x, center = "X/2")
  val2 <- p_wc(q = x, center = "X/2")
  expect_equal(val1[[1]], val2)

  x <- 1.8
  val1 <- stats::integrate(d_wc, lower = x, upper = Inf, center = "X/2")
  val2 <- 1 - p_wc(q = x, center = "X/2")
  expect_equal(val1[[1]], val2)
})

test_that("cdf and CI work out", {
  expect_equal(p_wc(q = -2) * 2, wc_alpha(width = 2))
  expect_equal((1 - p_wc(q = 2)) * 2, wc_alpha(width = 2))

  expect_equal(p_wc(q = -2, center = "X/2") * 2, wc_alpha(width = 2, center = "X/2"))
  expect_equal((1 - p_wc(q = 2, center = "X/2")) * 2, wc_alpha(width = 2, center = "X/2"))
})
