# Augmented variance

If wt is an integer, this is var(c(x, rep(A, wt))). But we generalize it
to non-integer wt as well.

## Usage

``` r
var_aug(x, A, wt = 1)
```

## Arguments

- x:

  Data

- A:

  Prior value

- wt:

  Weight for A. Don't touch this unless you know what you are doing.

## Value

The augmented variance

## Author

David Gerard

## Examples

``` r
var_aug(c(1, 2, 3), 4, wt = 10)
#> [1] 0.9358974
```
