test_that("idist agrees with other implementations", {

  ## Density ----
  expect_equal(
    didist(x = 1, center = 2, scale = 3, fam = "normal"),
    dinvnorm(x = 1, imean = 2, isd = 3)
  )

  expect_equal(
    didist(x = 1, center = 2, scale = 3, fam = "cauchy"),
    didist(x = 1, center = 2, scale = 3, ddist = dt, df = 1)
  )

  ## CDF ----
  expect_equal(
    pidist(q = 1, center = 2, scale = 3, fam = "normal"),
    pinvnorm(q = 1, imean = 2, isd = 3)
  )

  expect_equal(
    pidist(q = 1, center = 2, scale = 3, fam = "cauchy"),
    pidist(q = 1, center = 2, scale = 3, pdist = pt, df = 1)
  )

  expect_equal(
    pidist(q = -1, center = 2, scale = 3, fam = "normal"),
    pinvnorm(q = -1, imean = 2, isd = 3)
  )

  expect_equal(
    pidist(q = -1, center = 2, scale = 3, fam = "cauchy"),
    pidist(q = -1, center = 2, scale = 3, pdist = pt, df = 1)
  )

  ## Quantile ----
  expect_equal(
    qidist(p = 0.1, center = 2, scale = 3, fam = "normal"),
    qinvnorm(p = 0.1, imean = 2, isd = 3)
  )

  expect_equal(
    qidist(p = 0.1, center = 2, scale = 3, fam = "cauchy"),
    qidist(p = 0.1, center = 2, scale = 3, qdist = qt, pdist = pt, df = 1)
  )
})


test_that("Inverse of central Cauchy is Cauchy", {

  pvec <- seq(0.1, 0.9, by = 0.1)

  expect_equal(
    qt(p = pvec, df = 1),
    qidist(p = pvec, center = 0, scale = 1, fam = "cauchy")
  )

  qvec <- seq(-3, 3)

  expect_equal(
    pt(q = qvec, df = 1),
    pidist(q = qvec, center = 0, scale = 1, fam = "cauchy")
  )
})

test_that("idist gives credible intervals about same as n=1 confidence intervals", {
  A <- 0
  X <- 4
  alpha <- 0.011

  ## In uniform case, worst-case nu is 0.5
  ## Asymptotics are worse for uniform case
  ci1(x = X, A = A, family = "uniform", level = 1 - alpha)
  qidist(p = c(alpha / 2, 1 - alpha / 2), center = 1 / X, scale = 2 / abs(X), fam = "uniform")

  ## In Cauchy case, worst-case nu is 1
  ci <- ci1(x = X, A = A, family = "normal", level = 1 - alpha)[, c("lower", "upper")]
  attributes(ci) <- NULL
  expect_equal(
    ci,
    qidist(p = c(alpha / 2, 1 - alpha / 2), center = 1 / X, scale = 1 / abs(X), fam = "normal"),
    tolerance = 1e-5
  )

  ## In normal case, worst-case nu is 1
  ci <- ci1(x = X, A = A, family = "cauchy", level = 1 - alpha)[, c("lower", "upper")]
  attributes(ci) <- NULL
  expect_equal(
    ci,
    qidist(p = c(alpha / 2, 1 - alpha / 2), center = 1 / X, scale = 1 / abs(X), fam = "cauchy"),
    tolerance = 1e-6
  )
})

test_that("qidist and pidist are inverses", {
  pvec <- seq(0.01, 0.99, by = 0.01)
  qvec <- qidist(p = pvec, center = 4, scale = 3, fam = "uniform")
  expect_equal(
    pidist(q = qvec, center = 4, scale = 3, fam = "uniform"),
    pvec
  )

  qvec <- qidist(p = pvec, center = 4, scale = 3, fam = "normal")
  expect_equal(
    pidist(q = qvec, center = 4, scale = 3, fam = "normal"),
    pvec
  )

  qvec <- qidist(p = pvec, center = 4, scale = 3, fam = "cauchy")
  expect_equal(
    pidist(q = qvec, center = 4, scale = 3, fam = "cauchy"),
    pvec
  )
})

test_that("rdist agrees with pdist and qdist", {

  ## Uniform ----
  x <- ridist(n = 10000, center = 1/4, scale = 1/2, fam = "uniform")
  pvec <- seq(0.1, 0.9, by = 0.1)
  qvec <- quantile(x = x, probs = pvec)
  names(qvec) <- NULL
  expect_equal(
    qvec,
    qidist(p = pvec, center = 1/4, scale = 1/2, fam = "uniform"),
    tolerance = 0.1
  )

  mean(x < -10)
  pidist(q = -10, center = 1/4, scale = 1/2, fam = "uniform")

  ## Normal ---
  x <- ridist(n = 100000, center = 1/4, scale = 1/2, fam = "normal")
  pvec <- seq(0.1, 0.9, by = 0.1)
  qvec <- quantile(x = x, probs = pvec)
  names(qvec) <- NULL
  expect_equal(
    qvec,
    qidist(p = pvec, center = 1/4, scale = 1/2, fam = "normal"),
    tolerance = 0.1
  )

  mean(x < -10)
  pidist(q = -10, center = 1/4, scale = 1/2, fam = "normal")
})

test_that("Integral of didist is 1 over R", {
  X <- 4
  f <- function(x) {
    didist(x = x, center = 1 / X, scale = 2 / abs(X), fam = "uniform")
  }
  expect_equal(
    integrate(f = f, lower = -Inf, upper = Inf)[[1]],
    1
  )

  f <- function(x) {
    didist(x = x, center = 1 / X, scale = 2 / abs(X), fam = "normal")
  }
  expect_equal(
    integrate(f = f, lower = -Inf, upper = Inf)[[1]],
    1
  )

  f <- function(x) {
    didist(x = x, center = 1 / X, scale = 2 / abs(X), fam = "cauchy")
  }
  expect_equal(
    integrate(f = f, lower = -Inf, upper = Inf)[[1]],
    1
  )


})


