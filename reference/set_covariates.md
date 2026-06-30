# Set Covariates of a Survival Model

Generate a survival distribution representing model predictions for a
specified cohort. The cohort can be defined by providing a data frame of
covariate values (for multiple subjects) or by providing covariate
values as named arguments (for a single subject).

## Usage

``` r
set_covariates(dist, data, ...)
```

## Arguments

- dist:

  a survfit or flexsurvreg object

- data:

  a data.frame representing subjets for which predictions will be
  generated

- ...:

  optional argument representing covariate values to generate
  predictions for, can be used instead of data argument

## Value

a `surv_model` object

## Examples

``` r
library(flexsurv)
#> Loading required package: survival
fs1 <- flexsurvreg(
  Surv(rectime, censrec)~group,
  data=flexsurv::bc,
  dist = "llogis"
)
good_model <- set_covariates(fs1, group = "Good")
cohort <- data.frame(group=c("Good", "Good", "Medium", "Poor"))
mixed_model <- set_covariates(fs1, data = cohort)
```
