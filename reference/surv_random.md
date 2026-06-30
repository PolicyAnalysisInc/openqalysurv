# Generate Random Survival Times

Sample random event times from any survival distribution using
inverse-transform sampling via `surv_quantile`.

## Usage

``` r
surv_random(x, n, ...)
```

## Arguments

- x:

  A `surv_dist` object

- n:

  Number of random samples to generate

- ...:

  additional arguments passed to `surv_quantile`

## Value

A numeric vector of random event times

## Examples

``` r
dist1 <- define_surv_param('exp', rate = 0.05)
samples <- surv_random(dist1, 100)
```
