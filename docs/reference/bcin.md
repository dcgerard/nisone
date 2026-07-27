# Bayesian credible interval when n \>= 1

Let \\X_i \sim N(\mu, \sigma^2)\\. Let \\Y_i = X_i - A\\, \\\beta =
\mu - A\\, and \\\nu = \|\mu - A\|/\sigma\\. Then, for a fixed value of
\\\nu\\, this comes up with a \\(1-\alpha)100\\\\ posterior credible
interval of \\\mu\\ where the prior is \\\pi(\beta) = \|\beta\|^{-1}\\
or \\\pi(\mu) = \|\mu - A\|^{-1}\\. The full posterior distribution is
generalized inverse normal (implemented by
[`ginvnorm`](https://dcgerard.github.io/nisone/reference/ginvnorm.md))
with shape parameter \\n + 1\\, inverse mean of `sum(y) / sum(y^2)`, and
inverse variance `1 / (nu^2 * sum(y^2))`.

## Usage

``` r
bcin(x, A = 0, nu = 1, level = 0.95)
```

## Arguments

- x:

  a double

- A:

  The prior center.

- nu:

  The fixed value of nu.

- level:

  The level of the credible interval

## Value

The credible interval of the provided level.

## See also

[`ginvnorm()`](https://dcgerard.github.io/nisone/reference/ginvnorm.md):
The generalized inverse normal distribution.

## Author

David Gerard

## Examples

``` r
set.seed(1)
x <- stats::rnorm(4, mean = 1, sd = 1)
bcin(x = x, A = 3)
#>      lower      upper 
#> -0.1186979  2.0748793 

```
