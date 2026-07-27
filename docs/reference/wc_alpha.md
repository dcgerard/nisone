# Worst case exclusion probability of a n=1 confidence interval

This assumes the random variable X follows a normal distribution.

## Usage

``` r
wc_alpha(width, center = c("X", "ave"))
```

## Arguments

- width:

  The half-width of the interval, in units of \|X\|. This needs to be at
  least 1 if `center = "X"`, or at least 0.5 if `center = "ave"`.

- center:

  What is the center of the interval? Either `"X"` or `"ave"`.

## Value

The worst case probability (a priori) that a confidence interval of
half-width `width` will not capture the mean.

## Author

David Gerard

## Examples

``` r
wc_alpha(width = 2, center = "X")
#> [1] 0.242164
wc_alpha(width = 2, center = "ave")
#> [1] 0.2278774
```
