# General inverse distribution

Let \\X\\ come from some location-scale family \$\$f(x) =
\frac{1}{\sigma}\rho\left(\frac{X - \mu}{\sigma}\right),\$\$ for some
standard distribution \\\rho()\\. You provide \\\rho()\\, \\\mu\\ and
\\\sigma\\ and these functions will provide the density, distribution,
quantile, and random generation for \\Z = 1/X\\. Convenience arguments
for the normal, Cauchy, and uniform are provided.

## Usage

``` r
didist(
  x,
  center = 0,
  scale = 1,
  fam = c("normal", "cauchy", "uniform"),
  ddist = NULL,
  log = FALSE,
  ...
)

pidist(
  q,
  center = 0,
  scale = 1,
  fam = c("normal", "cauchy", "uniform"),
  pdist = NULL,
  lower.tail = TRUE,
  log.p = FALSE,
  ...
)

qidist(
  p,
  center = 0,
  scale = 1,
  fam = c("normal", "cauchy", "uniform"),
  qdist = NULL,
  pdist = NULL,
  ...
)

ridist(
  n,
  center = 0,
  scale = 1,
  fam = c("normal", "cauchy", "uniform"),
  rdist = NULL,
  ...
)
```

## Arguments

- x, q:

  vector of quantiles

- center:

  The center parameter \\\mu\\.

- scale:

  The scale parameter \\\sigma\\.

- fam:

  One of `"normal"`, `"cauchy"`, or `"uniform"`. If `ddist`, `pdist`,
  `qdist`, or `rdist` are specified then this argument is ignored.

- ddist:

  Density function of standard distribution \\\rho()\\. Should have a
  `log` argument.

- log, log.p:

  logical; if `TRUE`, probabilities p are given as log(p).

- ...:

  Additional arguments for `ddist`, `pdist`, `qdist`, and `rdist`.

- pdist:

  Cumulative distribution function of standard distribution.

- lower.tail:

  logical; if `TRUE` (default), probabilities are P(X\<=x) otherwise,
  P(X\>x).

- p:

  vector of probabilities

- qdist:

  Quantile function of standard distribution.

- n:

  sample size

- rdist:

  Random generation of standard distribution.

## Value

Either a random sample (`ridist`), the density (`didist`), the tail
probability (`pidist`), or the quantile (`qidist`) of the inverse
distribution.

## Functions

- `didist()`: Density function.

- `pidist()`: Distribution function.

- `qidist()`: Quantile function.

- `ridist()`: Random generation.

## Author

David Gerard

## Examples

``` r
didist(1, fam = "normal")
#> [1] 0.2419707
didist(1, center = 3, scale = 2, ddist = dt, df = 1)
#> [1] 0.07957747

pidist(1, fam = "normal")
#> [1] 0.6586553
pidist(1, center = 3, scale = 2, pdist = pt, df = 1)
#> [1] 0.937167

qidist(0.1, fam = "normal")
#> [1] -3.947154
qidist(0.1, center = 3, scale = 2, qdist = qt, pdist = pt, df = 1)
#> [1] -0.2427205
```
