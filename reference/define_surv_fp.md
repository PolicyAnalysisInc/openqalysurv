# Define Fractional Polynomial Survival Distribution

Define a fractional polynomial (FP) survival distribution on the
log-hazard scale, as used in network meta-analysis of survival data.
Supports FP1 (one power) and FP2 (two powers) models.

## Usage

``` r
define_surv_fp(betas, powers)
```

## Arguments

- betas:

  numeric vector of coefficients. Length must equal `length(powers) + 1`
  (intercept plus one per power).

- powers:

  numeric vector of fractional polynomial powers (length 1 for FP1,
  length 2 for FP2).

## Value

a `surv_fp` object.

## Details

The log-hazard is modeled as:

- FP1: `log h(t) = beta0 + beta1 * t^(p1)`

- FP2: `log h(t) = beta0 + beta1 * t^(p1) + beta2 * t^(p2)`

where `t^0 = log(t)` by convention, and repeated powers use the
Royston-Altman convention of multiplying by `log(t)`.

Survival probabilities are computed as `S(t) = exp(-H(t))` where `H(t)`
is the cumulative hazard obtained by numerical integration of the hazard
function.

The standard power set is {-2, -1, -0.5, 0, 0.5, 1, 2, 3}, but any
real-valued powers are accepted.

## References

Jansen, J. P. (2011). Network meta-analysis of survival data with
fractional polynomials. BMC Medical Research Methodology, 11, 61.

Royston, P. and Altman, D.G. (1994). Regression using fractional
polynomials of continuous covariates: parsimonious parametric modelling.
Journal of the Royal Statistical Society Series C, 43, 429-467.

## Examples

``` r

# FP1 model
define_surv_fp(betas = c(-3, 0.5), powers = c(1))
#> An FP1 survival distribution (beta0 = -3.0, beta1 = 0.5, powers = [1]).

# FP2 model
define_surv_fp(betas = c(-3, 0.5, -0.2), powers = c(0, 1))
#> An FP2 survival distribution (beta0 = -3.0, beta1 = 0.5, beta2 = -0.2, powers = [0, 1]).
```
