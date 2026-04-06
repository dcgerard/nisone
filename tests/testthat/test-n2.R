test_that("n2 coverage is true", {
  skip("simulations")
  TOL <- 1e-6
  n <- 2
  nu <- augtbounds$nu[augtbounds$n == n & abs(augtbounds$alpha - 0.05) < TOL]
  eta <- augtbounds$eta[augtbounds$n == n & abs(augtbounds$alpha - 0.05) < TOL]

  mu <- 1
  sig <- mu / nu

  rout <- replicate(n = 100000, expr = {
    x <- stats::rnorm(n = n, mean = mu, sd = sig)
    newt <- aug_t(x = x)
    oldt <- t.test(x = x)$conf.int
    new_cover <- newt[[1]] < mu && mu < newt[[2]]
    old_cover <- oldt[[1]] < mu && mu < oldt[[2]]
    new_width <- newt[[2]] - newt[[1]]
    old_width <- oldt[[2]] - oldt[[1]]
    c(new_cover = new_cover, old_cover = old_cover, new_width = new_width, old_width = old_width)
  })

  rownames(rout)

  ## Coverage
  mean(rout[1, ])
  mean(rout[2, ])

  ## Width
  mean(rout[3, ]^2 / 4)
  (n + nu^2) / (n + 1)^2 * eta^2 * sig^2

  mean(rout[4, ]^2 / 4)
  qt(p = 0.975, df = n - 1)^2 * sig^2 / n
})

test_that("worst alpha is correct", {
  set.seed(1)
  n <- 2
  alpha <- 0.05
  eout <- eta_alpha(alpha = alpha, n = n)
  nu <- eout[["nu"]]
  eta <- eout[["eta"]]
  nsim <- 100000
  z <- stats::rnorm(n = nsim)
  w <- stats::rchisq(n = nsim, df = n - 1)
  t2 <- (n * z - sqrt(2) * nu)^2 / ((n + 1) * w + (z + sqrt(n) * nu)^2)
  expect_equal(mean(t2 > eta^2), alpha, tolerance = 1e-2)
})

test_that("w_z is correct", {
  z <- 1
  n <- 3
  eta2 <- 2
  nu <- 1.5

  expect_equal(
    w_z(z = z, n = n, eta2 = eta2, nu = nu),
    (n * z - sqrt(n) * nu)^2 / (eta2 * (n + 1)) - (z + sqrt(n) * nu)^2 / (n + 1)
  )
})

test_that("old and new augtbounds are same", {
  load("./augtbounts_old.RData")
  expect_equal(augtbounds$eta, augtbounds_old$eta, tolerance = 1e-3)
})

test_that("w_z hasn't changed with weight when wt = 1", {
  w_z_old <- function(z, n, eta2, nu) {
    ((n^2 - eta2) * z^2 - 2 * (n^(1.5) + sqrt(n) * eta2) * nu * z + (1 - eta2) * n * nu^2) / ((n + 1) * eta2)
  }
  expect_equal(
    w_z(z = 1.3, n = 3, eta2 = 1.1, nu = 0.9, wt = 1),
    w_z_old(z = 1.3, n = 3, eta2 = 1.1, nu = 0.9)
  )

  expect_equal(
    w_z(z = 0.9, n = 10, eta2 = 0.1, nu = 1, wt = 1),
    w_z_old(z = 0.9, n = 10, eta2 = 0.1, nu = 1)
  )
})

test_that("mean_aug() and var_aug() work", {
  set.seed(10)
  x <- rnorm(4)
  A <- 3
  wt <- 2
  expect_equal(
    mean_aug(x = x, A = A, wt = 2),
    mean(c(x, rep(A, wt)))
  )
  expect_equal(
    var_aug(x = x, A = A, wt = 2),
    stats::var(c(x, rep(A, wt)))
  )
})

test_that("wt = 1 gives same augbounds as before", {
  expect_equal(
    eta_alpha(alpha = 0.05, n = 4, wt = 1)[["eta"]],
    augtbounds$eta[augtbounds$alpha == 0.05 & augtbounds$n == 4]
  )

  expect_equal(
    eta_alpha(alpha = 0.01, n = 20, wt = 1)[["eta"]],
    augtbounds$eta[augtbounds$alpha == 0.01 & augtbounds$n == 20]
  )

  expect_equal(
    eta_alpha(alpha = 0.1, n = 80, wt = 1)[["eta"]],
    augtbounds$eta[augtbounds$alpha == 0.1 & augtbounds$n == 80]
  )

  expect_equal(
    eta_alpha(alpha = 0.2, n = 3, wt = 1)[["eta"]],
    augtbounds$eta[augtbounds$alpha == 0.2 & augtbounds$n == 3]
  )
})


test_that("Coverage is correct for larger weight", {
  set.seed(1)
  n <- 2
  wt <- 10
  alpha <- 0.05
  mu <- 3
  A <- 1
  nrep <- 1000000

  eout <- eta_alpha(alpha = alpha, n = n, wt = wt)
  eta <- eout[["eta"]]
  nu <- eout[["nu"]]

  sigma <- abs(mu - A) / nu

  x <- matrix(stats::rnorm(n = n * nrep, mean = mu, sd = sigma), ncol = n)

  mean_vec <- apply(X = x, MARGIN = 1, FUN = mean_aug, A = A, wt = wt)
  var_vec <- apply(X = x, MARGIN = 1, FUN = var_aug, A = A, wt = wt)

  lower <- mean_vec - eta * sqrt(var_vec) / sqrt(n + wt)
  upper <- mean_vec + eta * sqrt(var_vec) / sqrt(n + wt)

  nfail <- sum(mu < lower | mu > upper)

  binom.test(x = nfail, n = nrep)
})




