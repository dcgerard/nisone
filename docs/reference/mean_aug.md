# Augmented mean

If wt is an integer, this is mean(c(x, rep(A, wt))). But we generalize
it to non-integer wt as well.

## Usage

``` r
mean_aug(x, A, wt = 1)
```

## Arguments

- x:

  Data

- A:

  Prior value

- wt:

  Weight for A. Don't touch this unless you know what you are doing.

## Value

The augmented mean

## Author

David Gerard

## Examples

``` r
mean_aug(c(1, 2, 3), 4, wt = 10)
#> [1] 3.538462
```
