# Evaluate Survival Quantiles

Find the time t at which the survival probability equals a given value.
That is, find t such that S(t) = p.

## Usage

``` r
surv_quantile(x, probs, ...)
```

## Arguments

- x:

  A `surv_dist` object

- probs:

  A numeric vector of survival probabilities in \[0, 1\]

- ...:

  additional arguments passed to methods

## Value

A numeric vector of times

## Examples

``` r
dist1 <- define_surv_param('exp', rate = 0.12)
surv_quantile(dist1, c(0.9, 0.5, 0.1))
#> [1]  0.8780043  5.7762265 19.1882091
```
