# Generating Predictions

Any openqalysurv survival object, or supported third-party distribution
object, can be used to generate predicted values

## Survival Probabilities

Predicted survival probabilities can be generated for a given survival
distribution and set of times using the `surv_prob` function.

``` r

  dist1 <- define_surv_param('exp', rate = 0.04)
  surv_prob(dist1, c(0, 1, 2, 3))
#> [1] 1.0000000 0.9607894 0.9231163 0.8869204
```

## Event Probabilities

Event probabilities can be generated for a given survival distributon
and set of time intervals using the `event_prob` function.

``` r

  dist1 <- define_surv_cure('exp', theta = 0.21, rate = 0.04, mixture = TRUE)
  event_prob(dist1, c(0, 1, 2, 3), c(1, 2, 3, 4))
#> [1] 0.03097634 0.03071312 0.03044387 0.03016860
```
