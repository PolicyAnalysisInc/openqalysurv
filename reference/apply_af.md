# Apply Acceleration Factor

Apply an acceleration factor to proportionally increase or reduce the
time-to-event of a survival distribution.

## Usage

``` r
apply_af(dist, af, log_af = FALSE)
```

## Arguments

- dist:

  a survival distribution

- af:

  an acceleration factor to be applied to survival distribution

- log_af:

  optional argument (defaults to `FALSE`) to indicate that provided
  acceleration factor is on log scale

## Value

A `surv_aft` object.

## Examples

``` r

dist1 <- define_surv_param("exp", rate = 0.25)
aft_dist <- apply_af(dist1, 1.5)
```
