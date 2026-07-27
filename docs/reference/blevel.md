# True level of the Bayesian credible interval at a given coefficient of variation

True level of the Bayesian credible interval at a given coefficient of
variation

## Usage

``` r
blevel(level, nu)
```

## Arguments

- level:

  The credible interval level.

- nu:

  (mu - A)/sigma

## Value

The true level of the (1-alpha) credible interval.

## Author

David Gerard

## Examples

``` r
blevel(0.95, 0.5)
#> [1] 0.9634901
blevel(0.95, 1)
#> [1] 0.95
blevel(0.95, 1.5)
#> [1] 0.9599223

```
