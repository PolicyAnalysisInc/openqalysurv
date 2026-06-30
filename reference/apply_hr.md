# Apply Hazard Ratio

Apply hazard ratio to a distribution to proportionally reduce or
increase the hazard rate.

## Usage

``` r
apply_hr(dist, hr, log_hr = FALSE)
```

## Arguments

- dist:

  a survival distribution

- hr:

  a hazard ratio to be applied to survival distribution

- log_hr:

  optional argument (defaults to `FALSE`) to indicate that provided
  hazard ratio is on log scale

## Value

A `surv_ph` object.

## Examples

``` r

ph_dist <- apply_hr(
 define_surv_param("exp", rate = 0.25),
 0.5
)
```
