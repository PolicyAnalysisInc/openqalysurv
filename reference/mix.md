# Mix Two or More Survival Distributions

Mix a set of survival distributions using the specified weights.

## Usage

``` r
mix(dist1, weight1, dist2, weight2, ...)
```

## Arguments

- dist1:

  first survival distribution to mix

- weight1:

  probability weight for first distribution

- dist2:

  second survival distribution to mix

- weight2:

  probability weight for second distribution

- ...:

  additional distributions and weights

## Value

A `surv_mix` object.

## Examples

``` r

dist1 <- define_surv_param("exp", rate = .5)
dist2 <- define_surv_param("gompertz", rate = .5, shape = 1)
pooled_dist <- mix(dist1, 0.25, dist2, 0.75)
```
