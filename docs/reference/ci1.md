# n=1 confidence interval

When you have one observation, \\X\\ from some symmetric distribution
with center \\\mu\\, for a pre-specified \\A\\, produces confidence
intervals for the center form \\X \pm \eta \|X - A\|\\ (`type = "x"`) or
\\(X + A)/2 \pm \eta \|X - A\|\\ (`type = "ave"`). The average intervals
have smaller width on average, and so are the default. We allow \\X\\ to
either come from a normal distribution, a Cauchy, or a uniform. Note
that if you use `family = "uniform"` then these intervals are valid
confidence intervals for the median/mode for *any* symmetric unimodal
density, which I think is really amazing.

## Usage

``` r
ci1(
  x,
  A = 0,
  type = c("ave", "x"),
  family = c("normal", "cauchy", "uniform"),
  level = 0.95
)
```

## Arguments

- x:

  A vector of single observations.

- A:

  Your "prior knowledge" or "augmented data value".

- type:

  Either centered at x or the average of x and A

- family:

  Either the normal distribution, Cauchy, or uniform.

- level:

  The level of the confidence interval.

## Value

A matrix of with 4 columns: the data, the center, the lower bound, the
upper bound.

## References

- Blachman, N., & Machol, R. (1987). Confidence intervals based on one
  or more observations. *IEEE Transactions on Information Theory*,
  33(3), 373-382.
  [doi:10.1109/TIT.1987.1057306](https://doi.org/10.1109/TIT.1987.1057306)

## Author

David Gerard

## Examples

``` r
ci1(c(1, 2, 10), type = "x")
#>       x center     lower     upper
#> [1,]  1      1  -8.67885  10.67885
#> [2,]  2      2 -17.35770  21.35770
#> [3,] 10     10 -86.78850 106.78850
ci1(c(1, 2, 10), type = "ave")
#>       x center      lower     upper
#> [1,]  1    0.5  -9.152953  10.15295
#> [2,]  2    1.0 -18.305907  20.30591
#> [3,] 10    5.0 -91.529534 101.52953
```
