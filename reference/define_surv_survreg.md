# Define Parametric Distribution from survreg Parameters

Define a parametric survival distribution using parameter estimates from
R's
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
output. Parameters are converted from the AFT parameterization to
flexsurv's native parameterization.

## Usage

``` r
define_surv_survreg(distribution, intercept, scale = 1)
```

## Arguments

- distribution:

  one of `"exponential"`, `"weibull"`, `"lognormal"`, or
  `"loglogistic"`.

- intercept:

  the intercept from the survreg model (`coef(model)["(Intercept)"]`).

- scale:

  the scale parameter from the survreg model (`model$scale`). Defaults
  to 1 (used for exponential).

## Value

a `surv_parametric` object.

## Details

The following conversions are applied:

|  |  |  |
|----|----|----|
| **survreg distribution** | **flexsurv distribution** | **Conversion** |
| exponential | exp | `rate = exp(-intercept)` |
| weibull | weibull | `shape = 1/scale`, `scale = exp(intercept)` |
| lognormal | lnorm | `meanlog = intercept`, `sdlog = scale` |
| loglogistic | llogis | `shape = 1/scale`, `scale = exp(intercept)` |

### Survival function formulations

In terms of survreg parameters:

- **Exponential**: `S(t) = exp(-exp(-intercept) * t)`

- **Weibull**: `S(t) = exp(-(t / exp(intercept))^(1/scale))`

- **Lognormal**: `S(t) = 1 - Phi((log(t) - intercept) / scale)` where
  `Phi()` is the standard normal CDF

- **Log-logistic**: `S(t) = 1 / (1 + (t / exp(intercept))^(1/scale))`

## Examples

``` r

define_surv_survreg("exponential", intercept = 3)
#> An exponential distribution (rate = 0.0498).
define_surv_survreg("weibull", intercept = 3, scale = 0.8)
#> A Weibull (AFT) distribution (shape = 1.25, scale = 20.09).
```
