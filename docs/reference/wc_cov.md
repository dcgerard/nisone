# Worst-case coefficient of variation

This assumes the random variable X is normal. Intervals are of the form
center +/- width \* \|X\|.

## Usage

``` r
wc_cov(width, center = c("X", "ave"))
```

## Arguments

- width:

  The half-width of the interval, in units of \|X\|. This needs to be at
  least 1 if `center = "X"`, or at least 0.5 if `center = "ave"`.

- center:

  What is the center of the interval? Either `"X"` or `"ave"`.

## Value

Worst case coefficient of variation \\\nu = \|X - A\|/\sigma\\.

## Details

We use A = 0 for calculations, but the worst case COV doesn't change if
A is non-zero. So you can think of intervals of the form center +/-
width \* \|X - A\| where center is either X or (X + A)/2.

## Author

David Gerard

## Examples

``` r
wc_cov(width = 2, center = "X")
#> [1] 0.7861103
wc_cov(width = 2, center = "ave")
#> [1] 0.9193019
```
