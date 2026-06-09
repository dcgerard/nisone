test_that("tail approximation is accurate", {
  alpha <- 2
  mu <- 0
  tau <- 1
  qginvnorm(p = 0.99, alpha = alpha, mu = mu, tau = tau)

  pginvnorm(q = 20, alpha = alpha, mu = mu, tau = tau, lower.tail = FALSE)
  ginv_tail(q = 20, alpha = alpha, mu = mu, tau = tau)
})
