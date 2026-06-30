# Apply Odds Ratio

Apply an odds ratio to proportionally increase or reduce the odds of
survival

## Usage

``` r
apply_or(dist, or, log_or = FALSE)
```

## Arguments

- dist:

  a survival distribution

- or:

  an odds ratio to be applied to survival distribution

- log_or:

  optional argument (defaults to `FALSE`) to indicate that provided odds
  ratio is on log scale

## Value

A `surv_po` object.

## Examples

``` r

dist1 <- define_surv_param("exp", rate = 0.25)
po_dist <- apply_or(dist1, 1.12)
```
