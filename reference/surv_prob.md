# Evaluate Survival Probabilities

Generate survival probabilities for a survival distribution at the
specified times.

## Usage

``` r
surv_prob(x, time, ...)
```

## Arguments

- x:

  A `surv_dist` object

- time:

  A numeric vector of times

- ...:

  additional arguments passed to methods

## Value

A numeric vector of survival probabilities

## Examples

``` r
dist1 <- define_surv_param('exp', rate = 0.12)
surv_prob(dist1, c(0, 1, 2, 3))
#> [1] 1.0000000 0.8869204 0.7866279 0.6976763
  
```
