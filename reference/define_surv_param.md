# Define Parametric Distribution

Alias for
[`define_surv_flexsurv()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_flexsurv.md).
Define a parametric survival distribution using flexsurv's native
parameterization.

## Usage

``` r
define_surv_param(distribution, ...)
```

## Arguments

- distribution:

  a parametric survival distribution (see
  [`define_surv_flexsurv()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_flexsurv.md)
  for supported distributions).

- ...:

  distribution parameters (see
  [`define_surv_flexsurv()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_flexsurv.md)
  for details).

## Value

a `surv_parametric` object.

## Examples

``` r

define_surv_param(distribution = "exp", rate = .5)
#> An exponential distribution (rate = 0.5).
define_surv_param(distribution = "gompertz", rate = .5, shape = 1)
#> A Gompertz distribution (shape = 1.0, rate = 0.5).
```
