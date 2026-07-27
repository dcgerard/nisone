# Augmented t-interval

Calculates a t-interval using augmented data `c(x,A)`. The multiplier of
this interval bounds the level above `level`, so these intervals are
typically conservative. This method is very fast for `n <= 100` and
`level %in% c(0.8, 0.9, 0.95, 0.99)` because I saved those multipliers
in an internal dataset. But it can be slow (on the order of a second)
for other confidence levels or larger sample sizes.

## Usage

``` r
aug_t(x, A = 0, level = 0.95, wt = 1)
```

## Arguments

- x:

  The vector of data

- A:

  The prior mean

- level:

  The level of the interval

- wt:

  Weight for A. Don't touch this unless you know what you are doing.

## Value

The augmented t-interval of the specified level.

## Details

Suppose \\X_1,\ldots,X_n \sim N(\mu,\sigma^2)\\. Suppose we have prior
value \\A\\. This provides intervals of the form \$\$\hat{\mu} \pm \eta
\hat{\sigma}/\sqrt{n + 1},\$\$ where \\\hat{\mu}\\ and \\\hat{\sigma}\\
are the sample mean and sample standard deviation of the augmented data
\\X_1,\ldots,X_n,A\\, and \\\eta\\ is chosen large enough to maintain
the confidence level at all parameter values.

The `wt` argument allows for more copies of \\A\\ to be included in the
data augmentation. But it doesn't work well with more data augmentation
so you should not set it above 1. Though, you can set `wt` to be between
0 and 1 (to have less data augmentation) and this does seem to work
pretty well, but I haven't studied it extensively, so use at your own
risk.

Note that you have to choose \\A\\ *before* seeing the data. If you
choose it based on the data, then you no longer maintain the confidence
level.

## Author

David Gerard

## Examples

``` r
set.seed(1)
## A is exactly correct
x <- rnorm(10)
aug_t(x)
#> [1] -0.3828999  0.6232686
stats::t.test(x)$conf.int
#> [1] -0.4261948  0.6906003
#> attr(,"conf.level")
#> [1] 0.95

## A is very wrong
x <- rnorm(10)
aug_t(x, A = 100)
#> [1] -11.09738  29.73165
stats::t.test(x)$conf.int
#> [1] -0.5162399  1.0139298
#> attr(,"conf.level")
#> [1] 0.95
```
