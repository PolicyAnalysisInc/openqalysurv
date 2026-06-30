# Define Parametric Distribution (flexsurv Parameterization)

Define a parametric survival distribution using flexsurv's native
parameterization. A complete listing of supported distributions is
provided in the details section.

## Usage

``` r
define_surv_flexsurv(distribution, ...)
```

## Arguments

- distribution:

  a parametric survival distribution.

- ...:

  additional distribution parameters (see details section below)

## Value

a `surv_parametric` object.

## Details

Supported distributions are listed in the table below.

|  |  |  |  |
|----|----|----|----|
| **Distribution** | **Description** | **Parameters** | **Notes** |
| "exp" | Exponential | rate |  |
| "lnorm" | Lognormal | meanlog, sdlog |  |
| "llogis" | Log-Logistic | shape, scale |  |
| "weibull" | Weibull (AFT) | shape, scale |  |
| "weibullPH" | Weibull (PH) | shape, scale |  |
| "gompertz" | Gompertz | shape, rate |  |
| "gamma" | Gamma | shape, scale |  |
| "gengamma" | Generalized Gamma (stable) | mu, sigma, Q | Described in Prentice (1974) |
| "gengamma.orig" | Generalized Gamma (original) | shape, scale, k | Described in Stacy (1962) |
| "genf" | Generalized F (stable) | mu, sigma, Q, P | Described in Prentice (1975) |
| "genf.org" | Generalized F (original) | mu, sigma, s1, s2 | Described in Prentice (1975) |

### Survival function formulations

- **Exponential**: `S(t) = exp(-rate * t)`

- **Weibull (AFT)**: `S(t) = exp(-(t / scale)^shape)`

- **Weibull (PH)**: `S(t) = exp(-scale * t^shape)`

- **Lognormal**: `S(t) = 1 - Phi((log(t) - meanlog) / sdlog)` where
  `Phi()` is the standard normal CDF

- **Log-logistic**: `S(t) = 1 / (1 + (t / scale)^shape)`

- **Gompertz**: hazard `h(t) = rate * exp(shape * t)`, survival
  `S(t) = exp(-(rate / shape) * (exp(shape * t) - 1))`

- **Gamma**: No closed-form survival function; computed numerically via
  [`pgamma()`](https://rdrr.io/r/stats/GammaDist.html)

- **Generalized gamma (stable/original)**: No simple closed form; see
  Prentice (1974) and Stacy (1962)

- **Generalized F (stable/original)**: No simple closed form; see
  Prentice (1975)

## References

Stacy, E. W. (1962). A generalization of the gamma distribution. Annals
of Mathematical Statistics 33:1187-92.

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

R. L. Prentice (1975). Discrimination among some parametric models.
Biometrika 62(3):607-614.

## Examples

``` r

define_surv_flexsurv(distribution = "exp", rate = .5)
#> An exponential distribution (rate = 0.5).
define_surv_flexsurv(distribution = "gompertz", rate = .5, shape = 1)
#> A Gompertz distribution (shape = 1.0, rate = 0.5).
```
