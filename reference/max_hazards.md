# Maximum Hazards

Combine two or more survival distributions by taking the pointwise
maximum of their discrete hazards at each time step. The steps form a
fixed grid from time 0 with step `cycle_length / 10`, so the survival at
a time does not depend on which other times are requested.

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

  the time step the distribution is used with (default 1). `surv_prob`
  takes the maximum hazard over steps of `cycle_length / 10`, and
  `surv_quantile` searches in steps of `cycle_length`.

## Value

A `surv_max_haz` object.

## Examples

``` r

dist1 <- define_surv_param("exp", rate = 0.01)
dist2 <- define_surv_param("exp", rate = 0.02)
max_dist <- max_hazards(dist1, dist2)
```
