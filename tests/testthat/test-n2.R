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


test_that("old and new augtbounds are same", {
  load("./augtbounts_old.RData")
  expect_equal(augtbounds$eta, augtbounds_old$eta, tolerance = 1e-3)
})
