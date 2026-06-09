## Final optimal nu
f <- function(x, n) {
  (0.5 * n) * log(x) - .hyp1f1(a = 0.5 * n, b = 0.5, z = x, log = TRUE)
}

nuvec <- rep(NA_real_, length.out = 9)
for (i in seq_along(nuvec)) {
  cat(i, "\n")
  n <- i + 1
  oout <- stats::optim(
    par = 1,
    fn = f,
    lower = 0,
    upper = Inf,
    method = "L-BFGS-B",
    control = list(fnscale = -1),
    n = n
  )
  nuvec[[i]] <- oout$par
}

n <- 10
nuvec[[n - 1]]
n / 4 + 1/12 + 2 / (81 * n) ## Asymptotic approximation
