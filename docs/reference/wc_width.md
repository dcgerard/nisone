# Given confidence level, provide half-width of CI

\\(1-\alpha)100\\\\ confidence intervals are of the form \\X \pm
\eta\|X-A\|\\ or \\(X + A)/2 \pm \eta\|X-A\|\\. You give \\\alpha\\ and
this function will provide \\\eta\\ if X follows a Cauchy, normal, or
uniform distribution.

## Usage

``` r
wc_width(
  alpha,
  center = c("X", "ave", "approx"),
  family = c("normal", "cauchy", "uniform")
)
```

## Arguments

- alpha:

  The confidence error probability. We produce a (1-alpha)100 percent
  confidence interval.

- center:

  Either X, ave, or the asymptotic width (approx).

- family:

  Either the normal, Cauchy, or uniform distribution.

## Value

Returns half-width of worst-case confidence interval in units of \|X -
A\|.

## Author

David Gerard

## Examples

``` r
wc_width(alpha = 0.2, center = "X")
#> [1] 2.420555
wc_width(alpha = 0.2, center = "ave")
#> [1] 2.312472
wc_width(alpha = 0.2, center = "approx")
#> [1] 2.419707

wc_width(alpha = 0.05, center = "X")
#> [1] 9.67885
wc_width(alpha = 0.05, center = "ave")
#> [1] 9.652953
wc_width(alpha = 0.05, center = "approx")
#> [1] 9.678829
```
