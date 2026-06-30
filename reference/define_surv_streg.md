# Define Parametric Distribution from Stata streg Parameters

Define a parametric survival distribution using parameter estimates from
Stata's `streg` output. Parameters are converted from the AFT or PH
parameterization to flexsurv's native parameterization.

## Usage

``` r
define_surv_streg(distribution, ..., metric = c("aft", "ph"))
```

## Arguments

- distribution:

  one of `"exponential"`, `"weibull"`, `"lognormal"`, `"loglogistic"`,
  `"ggamma"` (AFT metric), or `"gompertz"` (PH metric only).

- ...:

  distribution parameters from streg output. The `_cons` coefficient
  should be passed as `cons`. See details for parameter names by
  distribution.

- metric:

  either `"aft"` or `"ph"`, indicating which metric was used to fit the
  model.

## Value

a `surv_parametric` object.

## Details

### AFT metric parameters

|  |  |  |
|----|----|----|
| **Distribution** | **Parameters** | **Conversion** |
| exponential | `cons` | `rate = exp(-cons)` |
| weibull | `cons`, `ln_p` | `shape = exp(ln_p)`, `scale = exp(cons)` |
| lognormal | `cons`, `sigma` | `meanlog = cons`, `sdlog = sigma` |
| loglogistic | `cons`, `gamma` | `shape = gamma`, `scale = exp(cons)` |
| ggamma | `cons`, `sigma`, `kappa` | `mu = cons`, `sigma = sigma`, `Q = kappa` |

### PH metric parameters

|                  |                 |                                          |
|------------------|-----------------|------------------------------------------|
| **Distribution** | **Parameters**  | **Conversion**                           |
| exponential      | `cons`          | `rate = exp(cons)`                       |
| weibull          | `cons`, `ln_p`  | `shape = exp(ln_p)`, `scale = exp(cons)` |
| gompertz         | `cons`, `gamma` | `shape = gamma`, `rate = exp(cons)`      |

### Survival function formulations

#### AFT metric

- **Exponential**: `S(t) = exp(-exp(-cons) * t)`

- **Weibull**: `S(t) = exp(-(t / exp(cons))^exp(ln_p))`

- **Lognormal**: `S(t) = 1 - Phi((log(t) - cons) / sigma)` where `Phi()`
  is the standard normal CDF

- **Log-logistic**: `S(t) = 1 / (1 + (t / exp(cons))^gamma)`

- **Generalized gamma**: No simple closed form; `mu = cons`,
  `sigma = sigma`, `Q = kappa`

#### PH metric

- **Exponential**: `S(t) = exp(-exp(cons) * t)`

- **Weibull**: `S(t) = exp(-exp(cons) * t^exp(ln_p))`

- **Gompertz**: hazard `h(t) = exp(cons) * exp(gamma * t)`, survival
  `S(t) = exp(-(exp(cons) / gamma) * (exp(gamma * t) - 1))`

## Examples

``` r

define_surv_streg("exponential", cons = 3, metric = "aft")
#> An exponential distribution (rate = 0.0498).
define_surv_streg("weibull", cons = 3, ln_p = 0.405, metric = "aft")
#> A Weibull (AFT) distribution (shape = 1.5, scale = 20.1).
define_surv_streg("gompertz", cons = -3, gamma = 0.05, metric = "ph")
#> A Gompertz distribution (shape = 0.0500, rate = 0.0498).
```
