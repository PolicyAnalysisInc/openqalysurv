# Check for a Survival Distribution

Checks whether an object can be evaluated by
[`surv_prob()`](https://policyanalysisinc.github.io/openqalysurv/reference/surv_prob.md),
i.e. it is a survival distribution from openqalysurv, a `flexsurvreg`
fit, or a `survfit` object.

## Usage

``` r
is_surv_dist(x)
```

## Arguments

- x:

  an object to check

## Value

`TRUE` if `x` is a valid survival distribution, otherwise `FALSE`.
