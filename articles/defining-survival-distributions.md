# Defining Survival Distributions

The openqalysurv package provides functions for defining a wide variety
of survival distributions of different types.

## Parametric

The `define_surv_param` function allows you to define a parametric
survival model with specified distribution and parameter values.
Parameterization is based on the
[flexsurv](https://cran.r-project.org/web/packages/flexsurv/index.html)
package.

``` r

  param_surv <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  print(param_surv)
#> A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
```

## Mixture/Non-Mixture Cure

The `define_surv_cure` function allows you to define a parametric
mixture/non-mixture cure model with specified distribution, cure
fraction, and parameter values.

``` r

  cure_surv <- define_surv_cure('exp', theta = 0.21, rate = 0.04, mixture = TRUE)
  print(cure_surv)
#> An exponential mixture cure distribution (theta = 0.21, rate = 0.04).
```

## Royston & Parmar Spline

The `define_surv_spline` function allows you to define a Royston &
Parmar restricted cubic spline model with the specified scale, parameter
values, and knot times.

``` r

  spline_dist <- define_surv_spline(
    scale = 'hazard', # spline used to model log cumulative hazards
    -2.08, 2.75, 0.23, # parameters
    -1.62, 0.57, 1.191 # knot times
  )
  print(spline_dist)
#> A Royston & Parmar spline model of log cumulative hazard with 3 knots (gamma = [-2.08, 2.75, 0.23], knots = [-1.62, 0.57, 1.19]).
```

## Kaplan-Meier

The `define_surv_km` function allows you to define a survival
distribution based on the times and survival probabilities from a
Kaplan-Meier.

``` r

  km_output <- data.frame(
    times = c(0, 2, 5, 10, 17, 21),
    surv_prob = c(1, 0.98, 0.94, 0.89, 0.78, 0.78)
  )
  km_dist <- define_surv_km(km_output, "times", "surv_prob")
  print(km_dist)
#> A Kaplan-Meier distribution:
#>    time survival
#>   <dbl>    <dbl>
#> 1     0     1   
#> 2     2     0.98
#> 3     5     0.94
#> 4    10     0.89
#> 5    17     0.78
#> 6    21     0.78
```

## Life-Table

The `define_surv_lifetable` function allows you to define a survival
distribution based on a life-table, starting age, and gender mix.

``` r

  life_table <- data.frame(
    age_in_years = c(50, 51, 52, 53, 54, 55),
    p_death_male = c(0.0021, 0.0022, 0.0022, 0.0023, 0.0022, 0.0024),
    p_death_female = c(0.0018, 0.0018, 0.0018, 0.0019, 0.0019, 0.0020)
  )
  lifetable_dist <- define_surv_lifetable(
    life_table, # data.frame containing life-table
    51, # starting age
    0.48, # percent male
    age_col = 'age_in_years', # column in life-table containing age
    male_col = 'p_death_male', # column in life-table containing male p(death)
    female_col = 'p_death_female' # column in life-table containing female p(death)
  )
  print(lifetable_dist)
#> A life-table survival distribution (48.0% male, 52.0% female):
#>     age   male female
#>   <dbl>  <dbl>  <dbl>
#> 1    51 0.0022 0.0018
#> 2    52 0.0022 0.0018
#> 3    53 0.0023 0.0019
#> 4    54 0.0022 0.0019
#> 5    55 0.0024 0.002
```

## Custom Function

The `define_surv_func` function allows you to define a survival
distribution based on a custom function.

``` r

  custom_dist <- define_surv_func(function(t) pweibull(t, 1.2, 20.1, lower.tail = FALSE))
  print(custom_dist)
#> A survival distribution based on a custom function
#>   Function:
#>     function (t) 
#>     pweibull(t, 1.2, 20.1, lower.tail = FALSE)
```

## Third-Party Support

In addition to the functions provided by openqalysurv, you can use
survival distributions created from supported third-party packages.

``` r

  library(survival)
  library(flexsurv)
  survfit_dist <- survfit(Surv(recyrs, censrec)~1, data = bc)
  flexsurv_dist <- flexsurvreg(Surv(recyrs, censrec)~1, data = bc, dist = 'lnorm')
```
