# Confidence density based on n=1 confidence interval

Let \\X \sim N(\mu, \sigma^2)\\. Consider intervals of the form \\X \pm
\eta \|X\|\\ or \\X / 2 \pm \eta \|X\|\\. These intervals define a
confidence distribution of \\\mu\\. This is the density of the
confidence distribution of \\(\mu - X)/\|X\|\\ or \\(\mu - X/2)/\|X\|\\.
Please note that \\\mu\\ is random here not \\X\\. Also note that this
confidence distribution does not exist between the 25th and 75th
percentiles.

## Usage

``` r
d_wc(x, center = c("X", "ave"))
```

## Arguments

- x:

  The value at which to calculate the density. Only defined for
  `abs(x) >= 1` when `center = "X"`, and for `abs(x) >= 0.5` when
  `center = "ave"`.

- center:

  What is the center of the interval? Either `"X"` or `"X/2"` (`"ave"`).

## Value

n=1 confidence density.

## See also

[`p_wc()`](https://dcgerard.github.io/nisone/reference/p_wc.md)

## Author

David Gerard

## Examples

``` r
xseq <- seq(-10, 10, length.out = 500)
dseq <- d_wc(xseq)
plot(xseq, dseq, type = "l")


dseq <- d_wc(xseq, center = "ave")
plot(xseq, dseq, type = "l")

```
