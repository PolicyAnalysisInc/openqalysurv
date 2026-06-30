# Add Hazards

Combine two or more survival distributions as independent risks by
adding their hazards.

## Usage

``` r
add_hazards(dist1, dist2, ...)
```

## Arguments

- dist1:

  survival distribution to add

- dist2:

  second survival distribution to add

- ...:

  additional survival distributions to add

## Value

a `surv_add_haz` object

## Examples

``` r

dist1 <- define_surv_param("exp", rate = .125)
dist2 <- define_surv_param("weibull", shape = 1.2, scale = 50)
combined_dist <- add_hazards(dist1, dist2)
```
