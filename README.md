
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nisone

<!-- badges: start -->

[![R-CMD-check](https://github.com/dcgerard/nisone/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dcgerard/nisone/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/dcgerard/nisone/graph/badge.svg)](https://app.codecov.io/gh/dcgerard/nisone)
<!-- badges: end -->

Provides different interval estimates of a location parameter when the
sample size is one or more. These include classical methods when n=1,
new Bayesian analogues of these classical methods, and extensions of
these methods to larger sample sizes in the normal case. Other functions
calculate Bayes factors based on t-statistics, and implement the
(generalized) inverse normal distribution (and other inverse
distributions). See Gerard (2026) for details of these methods.

## Installation

You can install the development version of nisone from
[GitHub](https://github.com/dcgerard/nisone) with:

``` r
# install.packages("pak")
pak::pak("dcgerard/nisone", build_vignettes = TRUE)
```

# n=1 confidence interval

If we observe $X = 10$ and we have prior guess that the mean is around
5, then the resulting 95% CI based on a normal model is

``` r
library(nisone)
ci1(x = 10, A = 5)
#>       x center  lower upper
#> [1,] 10    7.5 -40.76 55.76
```

If we just want to assume that $X$ comes from *any* symmetric unimodal
distribution, a 95% CI for the mode is

``` r
ci1(x = 10, A = 5, family = "uniform")
#>       x center  lower upper
#> [1,] 10    7.5 -87.43 102.4
```

See more functionality by running

``` r
browseVignettes("nisone")
```

# References

- Gerard, D. (2026). Constructing and extending *n* = 1 Bayesian
  confidence intervals of the normal mean. *In preparation*.
