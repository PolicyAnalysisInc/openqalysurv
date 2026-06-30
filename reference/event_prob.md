# Evaluate Event Probabilities

Generate the conditional probability of an even during an interval of
time.

## Usage

``` r
event_prob(x, start, end, ...)
```

## Arguments

- x:

  A `surv_dist` object

- start:

  A numeric vector of interval start times

- end:

  A numeric vector of interval end times

- ...:

  additional arguments passed to methods

## Examples

``` r
dist1 <- define_surv_param('exp', rate = 0.12)
surv_prob(dist1, c(0, 1, 2, 3))
#> [1] 1.0000000 0.8869204 0.7866279 0.6976763
  
```
