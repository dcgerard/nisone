test_that("n2 coverage is true", {
  skip("simulations")
  nu <- 1.2
  n <- 2

  mu <- 1
  sig <- mu / nu

  rout <- replicate(n = 100000, expr = {
    x <- stats::rnorm(n = n, mean = mu, sd = sig)
    newt <- aug_t(x = x)
    oldt <- t.test(x = x)$conf.int
    ## fab <- fabCI::fabtCI(y = x)

    new_cover <- newt[[1]] < mu && mu < newt[[2]]
    old_cover <- oldt[[1]] < mu && mu < oldt[[2]]
    ## fab_cover <- fab[[1]] < mu && mu < fab[[2]]
    new_width <- newt[[2]] - newt[[1]]
    old_width <- oldt[[2]] - oldt[[1]]
    ## fab_width <- fab[[2]] - fab[[1]]
    c(new_cover = new_cover, old_cover = old_cover, new_width = new_width, old_width = old_width)
  })

  rownames(rout)

  ## Coverage
  mean(rout[1, ])
  mean(rout[2, ])

  ## Width
  mean(rout[3, ])
  mean(rout[4, ])
})

test_that("worst alpha is correct", {
  skip("simulations")
  n <- 2
  alpha <- 0.05
  eout <- eta_alpha(alpha = alpha, n = n)
  nu <- eout[["nu"]]
  eta <- eout[["eta"]]  ## should be about 5.2
  nsim <- 100000
  z <- stats::rnorm(n = nsim)
  w <- stats::rchisq(n = nsim, df = n - 1)
  t2 <- (n * z + nu)^2 / ((n + 1) * w + (z - nu)^2)
  mean(t2 > eta^2) ## should be alpha
})
