# Bayes factor from t-statistic

Given a user-provided prior on the standardized effect size, will
calculate the Bayes factor.

## Usage

``` r
bft(t, nu, nd, prior = NULL)
```

## Arguments

- t:

  The observed t-statistic.

- nu:

  The degrees of freedom of the t-statistic. This should be \\n-1\\ in
  the one-sample case. This should be \\n_1+n_2-2\\ in the two-sample
  case.

- nd:

  Effective sample size. This should be \\n\\ in the one-sample case.
  This should be \\(\frac{1}{n_1} + \frac{1}{n_2})^{-1}\\ in the
  two-sample case.

- prior:

  The prior on the standardized effect size, \\\delta\\. In the
  one-sample case this is the prior over \\\delta = \frac{\mu -
  \mu_0}{\sigma}\\, where \\\mu_0\\ is the null value. In the two-sample
  case this is the prior over \\\delta = \frac{\mu_1 - \mu_2}{\sigma}\\,
  where \\\mu_1\\ and \\\mu_2\\ are the means of the two-samples. The
  default prior is a Cauchy, but you can put in any density function you
  want. If your prior has pointmasses in it, this function won't work.

## Value

The Bayes factor to a corresponding t-statistic.

## Details

Let \\\delta\\ be the standardized effect size, let \\\sigma^2\\ be the
variance (assumed equal in two-sample case). In the one-sample case,
\\\delta = \frac{\mu - \mu_0}{\sigma}\\, where \\\mu_0\\ is the null
value. In the two-sample case, \\\delta = \frac{\mu_1 -
\mu_2}{\sigma}\\, where \\\mu_1\\ and \\\mu_2\\ are the means of the
two-samples. We place the prior \\1/\sigma^2\\ under the null and
\\\pi(\delta)/\sigma^2\\ under the alternative, for some arbitrary
density \\\pi(\cdot)\\. Given this setting, the Bayes factor is a
function of the \\t\\-statistic. We calculate it via numeric
integration.

## References

- Gronau, Q. F., Ly, A., & Wagenmakers, E. J. (2020). Informed Bayesian
  t-tests. *The American Statistician*.
  [doi:10.1080/00031305.2018.1562983](https://doi.org/10.1080/00031305.2018.1562983)

## Author

David Gerard

## Examples

``` r
# One sample t, n = 10, t-statistic = 2
bft(t = 2, nu = 10 - 1, nd = 10)
#> [1] 1.113217

# Two sample t, n1 = 10, n2 = 8, t-statistic = 2
bft(t = 2, nu = 10 + 8 - 2, nd = 1 / (1 / 10 + 1 / 8))
#> [1] 1.393369
```
