# Confidence distribution CDF based on n=1 confidence interval

Let \\X \sim N(\mu, \sigma^2)\\. Consider intervals of the form \\X \pm
\eta \|X\|\\ or \\X / 2 \pm \eta \|X\|\\. These intervals define a
confidence distribution of \\\mu\\. This is the CDF of the confidence
distribution of \\(\mu - X)/\|X\|\\ or \\(\mu - X/2)/\|X\|\\. Please
note that \\\mu\\ is random here not \\X\\. Also note that this
confidence distribution does not exist between the 25th and 75th
percentiles.

## Usage

``` r
p_wc(q, center = c("X", "ave"))
```

## Arguments

- q:

  The quantile. Only defined for `abs(q) >= 1` when `center = "X"`, and
  for `abs(q) >= 0.5` when `center = "ave"`.

- center:

  What is the center of the interval? Either `"X"` or `X/2`, (`"ave"`).

## Value

n=1 confidence distribution CDF.

## See also

[`d_wc()`](https://dcgerard.github.io/nisone/reference/d_wc.md)

## Author

David Gerard

## Examples

``` r
qseq <- seq(-10, 10, length.out = 100)
pseq <- p_wc(qseq)
plot(qseq, pseq, type = "l", xlab = "q", ylab = "CDF")

```
