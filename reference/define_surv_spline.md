# Define Royston & Parmar Spline Survival Distribution

Define Royston & Parmar restricted cubic spline parametric survival
distribution.

## Usage

``` r
define_surv_spline(scale, ...)

define_spline_survival(scale, ...)
```

## Arguments

- scale:

  "hazard", "odds", or "normal", as described in flexsurvspline. With
  the default of no knots in addition to the boundaries, these models
  reduce to the Weibull, log-logistic and log-normal respectively. The
  scale must be common to all times.

- ...:

  parameters and knot log times of spline distribution, which can be
  provided either in order starting with spline parameters followed by
  knot log times, or by names (e.g gamma1, gamma2, ... gammaN, knots1,
  knots2, ... knotsN). See examples below for named and unnamed calls.

## Value

a `surv_spline` object.

## Details

These models use restricted cubic splines to flexibly model a
transformation of the survival function as a function of log time, as
proposed by Royston & Parmar (2002) and implemented in
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).
Restricted cubic splines are piecewise cubic polynomials constrained to
be linear beyond the boundary knots, which avoids overfitting in the
tails where data are sparse. The `scale` parameter determines which
transformation of the survival function is modeled as the restricted
cubic spline: `"hazard"` fits the log cumulative hazard (proportional
hazards), `"odds"` fits the log cumulative odds (proportional odds), and
`"normal"` fits the probit (inverse normal) transformation. With no
internal knots beyond the boundaries, these restricted cubic spline
models reduce to the Weibull, log-logistic, and lognormal distributions
respectively.

## References

Royston, P. and Parmar, M. (2002). Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of
treatment effects. Statistics in Medicine 21(1):2175-2197.

## Examples

``` r

define_surv_spline(
 scale = 'hazard',
 -2.08, 2.75, 0.23, # parameters
 -1.62, 0.57, 1.191 # knot times
)
#> A Royston & Parmar spline model of log cumulative hazard with 3 knots (gamma = [-2.08, 2.75, 0.23], knots = [-1.62, 0.57, 1.19]).
```
