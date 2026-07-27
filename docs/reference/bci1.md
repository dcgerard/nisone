# Bayesian credible interval when n=1

Suppose \\X \sim N(\mu, \sigma^2)\\. This gives you credible intervals
of a specified level when we use prior \\\sigma = \|\mu\|\\ with
probability 1, and \\\pi(\mu) = \|\mu - A\|^{-1}\\. See the
[`n1post()`](https://dcgerard.github.io/nisone/reference/n1post.md) for
access to the full posterior distribution.

## Usage

``` r
bci1(x, A = 0, level = 0.95)
```

## Arguments

- x:

  a double

- A:

  The prior center.

- level:

  The level of the credible interval

## Value

The credible interval of the specified level.

## See also

[`n1post()`](https://dcgerard.github.io/nisone/reference/n1post.md):
Functions to access the full posterior distribution.

## Author

David Gerard

## Examples

``` r
bci1(x = 1, A = 0)
#>      x       med     lower    upper
#> [1,] 1 0.7098261 -9.151018 10.15462
```
