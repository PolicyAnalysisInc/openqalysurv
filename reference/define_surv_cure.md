# Define Parametric Mixture or Non-Mixture Cure Distribution

Define a parametric mixture or non-mixture cure distribution with given
parameter values. Uses the `flexsurvcure` package
([flexsurvcure::pmixsurv](https://rdrr.io/pkg/flexsurvcure/man/mixsurv.html)/[flexsurvcure::pnmixsurv](https://rdrr.io/pkg/flexsurvcure/man/nmixsurv.html))
for evaluation.

## Usage

``` r
define_surv_cure(distribution, theta, ..., mixture = TRUE)
```

## Arguments

- distribution:

  A parametric survival distribution. See details for a listing of valid
  distributions.

- theta:

  A numeric value in \[0, 1\] representing the cure fraction, i.e. the
  long-term probability of never experiencing the event.

- ...:

  Additional distribution parameters (see details section below)

- mixture:

  a logical determining whether a mixture or non-mixture model is being
  defined.

## Value

a `surv_dist_cure` object.

## Details

**Mixture cure model** (`mixture = TRUE`): \$\$S(t) = \theta + (1 -
\theta) \\ S_u(t)\$\$

**Non-mixture cure model** (`mixture = FALSE`): \$\$S(t) =
\theta^{F_u(t)}\$\$

where \\\theta\\ is the cure fraction (long-term survival probability),
\\S_u(t)\\ is the survival function for the uncured subgroup, and
\\F_u(t) = 1 - S_u(t)\\ is the corresponding CDF.

Supported distributions for the uncured subgroup are listed in the table
below.

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
| "gengamma" | Generalized Gamma (stable) | mu, sigma, Q | Parameterization from Prentice (1974) |
| "gengamma.orig" | Generalized Gamma (original) | shape, scale, k | Original parameterization from Stacy (1962) |
| "genf" | Generalized F (stable) | mu, sigma, Q, P | Stable reparameterization from Prentice (1975) |
| "genf.orig" | Generalized F (original) | mu, sigma, s1, s2 | Original parameterization described by Prentice (1975) |

## References

Lambert, P.C. (2007). Modeling of the cure fraction in survival studies.
The Stata Journal, 7(3), 351-375.

Stacy, E. W. (1962). A generalization of the gamma distribution. Annals
of Mathematical Statistics 33:1187-92.

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

R. L. Prentice (1975). Discrimination among some parametric models.
Biometrika 62(3):607-614.

## Examples

``` r

define_surv_cure(distribution = "exp", theta = 0.34, rate = .5)
#> An exponential mixture cure distribution (theta = 0.34, rate = 0.50).
define_surv_cure(
 distribution = "weibull",
 theta = 0.5, shape = 1.5,
 scale = 34.43,
 mixture = TRUE
)
#> A Weibull (AFT) mixture cure distribution (theta = 0.5, shape = 1.5, scale = 34.4).
```
