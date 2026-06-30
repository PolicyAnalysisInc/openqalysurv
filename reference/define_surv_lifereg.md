# Define Parametric Distribution from SAS PROC LIFEREG Parameters

Define a parametric survival distribution using parameter estimates from
SAS PROC LIFEREG output. Parameters are converted from the AFT
parameterization to flexsurv's native parameterization.

## Usage

``` r
define_surv_lifereg(distribution, intercept, scale = 1, shape = NULL)
```

## Arguments

- distribution:

  one of `"exponential"`, `"weibull"`, `"lnormal"`, `"llogistic"`, or
  `"gamma"`.

- intercept:

  the Intercept from the LIFEREG Parameter Estimates table.

- scale:

  the Scale parameter from the LIFEREG output. Defaults to 1 (used for
  exponential).

- shape:

  the Shape parameter from the LIFEREG output. Only used for
  `distribution = "gamma"`.

## Value

a `surv_parametric` object.

## Details

The following conversions are applied:

|  |  |  |
|----|----|----|
| **SAS DIST=** | **flexsurv distribution** | **Conversion** |
| EXPONENTIAL | exp | `rate = exp(-intercept)` |
| WEIBULL | weibull | `shape = 1/scale`, `scale = exp(intercept)` |
| LNORMAL | lnorm | `meanlog = intercept`, `sdlog = scale` |
| LLOGISTIC | llogis | `shape = 1/scale`, `scale = exp(intercept)` |
| GAMMA | gengamma | `mu = intercept`, `sigma = scale`, `Q = shape` |

### Survival function formulations

In terms of LIFEREG parameters:

- **Exponential**: `S(t) = exp(-exp(-intercept) * t)`

- **Weibull**: `S(t) = exp(-(t / exp(intercept))^(1/scale))`

- **Lognormal**: `S(t) = 1 - Phi((log(t) - intercept) / scale)` where
  `Phi()` is the standard normal CDF

- **Log-logistic**: `S(t) = 1 / (1 + (t / exp(intercept))^(1/scale))`

- **Generalized gamma**: No simple closed form; `mu = intercept`,
  `sigma = scale`, `Q = shape`; see Prentice (1974)

## References

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

## Examples

``` r

define_surv_lifereg("exponential", intercept = 3)
#> An exponential distribution (rate = 0.0498).
define_surv_lifereg("weibull", intercept = 3, scale = 0.8)
#> A Weibull (AFT) distribution (shape = 1.25, scale = 20.09).
define_surv_lifereg("gamma", intercept = 2.3, scale = 0.4, shape = -0.03)
#> A generalized gamma (stable) distribution (mu = 2.30, sigma = 0.40, Q = -0.03).
```
