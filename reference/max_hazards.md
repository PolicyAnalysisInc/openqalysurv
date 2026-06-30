# Maximum Hazards

Combine two or more survival distributions by taking the pointwise
maximum of their discrete hazards at each time step.

## Usage

``` r
max_hazards(dist1, dist2, ..., cycle_length = 1)
```

## Arguments

- dist1:

  first survival distribution

- dist2:

  second survival distribution

- ...:

  additional survival distributions

- cycle_length:

  step size used by `surv_quantile` for grid search (default 1). Has no
  effect on `surv_prob`.

## Value

A `surv_max_haz` object.

## Examples

``` r

dist1 <- define_surv_param("exp", rate = 0.01)
dist2 <- define_surv_param("exp", rate = 0.02)
max_dist <- max_hazards(dist1, dist2)
```
