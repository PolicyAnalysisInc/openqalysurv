# Modifying & Combining Distributions

Once you’ve created survival distributions, openqalysurv provides
various methods for combining modifying them.

## Apply Hazard Ratio:

Hazard ratios can be applied to a survival distribution using the
`apply_hr` function.

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- apply_hr(dist1, 0.85)
  print(dist2)
#> A proportional hazards survival distribution:
#>     * Hazard Ratio: 0.85
#>     * Baseline Distribution: A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
```

## Apply Odds Ratio

Odds ratios can be applied to a survival distribution using the
`apply_or` function.

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- apply_or(dist1, 0.45)
  print(dist2)
#> A proportional odds survival distribution:
#>     * Odds Ratio: 0.45
#>     * Baseline Distribution: A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
```

## Apply Accleration Factor

Accleration factors, which proportionally scale event times, can be
applied using the `apply_af` function.

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- apply_af(dist1, 0.75)
  print(dist2)
#> An accelerated failure time survival distribution:
#>     * Acceleration Factor: 0.75
#>     * Baseline Distribution: A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
```

## Apply Time Shift

Hazards for a distribution can be shifted forward or backward by a fixed
amount of time using the `apply_shift` function. Note that a positive
value of shift will move hazards forward and a negative value will shift
hazards backward (resulting in an initial period with no event risk).

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- apply_shift(dist1, 2.1)
  print(dist2)
#> A shifted survival distribution:
#>     * Shift: 2.1
#>     * Baseline Distribution: A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
```

## Mix Distributions

Survival distributions can be mixed according to specified weights using
the `mix` function.

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- define_surv_param('weibull', shape = 1.4, scale = 4.1)
  dist3 <- mix(dist1, 0.8, dist2, 0.2) # Mix of 80% dist1 and 20% dist2
  print(dist3)
#> A mixed survival distribution:
#>     * Distribution 1 (80%): A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
#>     * Distribution 2 (20%): A Weibull (AFT) distribution (shape = 1.4, scale = 4.1).
```

## Join Distributions

Survival distributions can be joined together at specified cutpoints
using the `join` function.

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- define_surv_param('weibull', shape = 1.4, scale = 4.1)
  dist3 <- join(dist1, 3, dist2) # Join dist2 and dist3 at time 3
  print(dist3)
#> A joined survival distribution:
#>     * Segment 1 (t = 0 - 3): A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
#>     * Segment 2 (t = 3 - ∞): A Weibull (AFT) distribution (shape = 1.4, scale = 4.1).
```

## Add Hazards

Survival distributions can be combined as independent risks using the
`add_hazards` function. This results in a survival distribution with
survival function equal to the product of the survival functions of
combined distributions.

``` r

  dist1 <- define_surv_param('weibull', shape = 1.2, scale = 5.3)
  dist2 <- define_surv_param('weibull', shape = 1.4, scale = 4.1)
  dist3 <- add_hazards(dist1, dist2) # combined hazards of dist1 and dist2
  print(dist3)
#> A survival distribution combining the hazards of:
#>     * Distribution 1: A Weibull (AFT) distribution (shape = 1.2, scale = 5.3).
#>     * Distribution 2: A Weibull (AFT) distribution (shape = 1.4, scale = 4.1).
```

## Set Covariates

Models fitted using `survfit` or `flexsurvreg` with predictors can be
used to generate a survival distribution representing a subject or set
of subjects using `set_covariates`.

``` r

library(flexsurv)
#> Loading required package: survival
library(tibble)

# Fit a model using bc dataset using good as predictor
print(as_tibble(bc))
#> # A tibble: 686 × 4
#>    censrec rectime group recyrs
#>      <int>   <int> <fct>  <dbl>
#>  1       0    1342 Good    3.68
#>  2       0    1578 Good    4.32
#>  3       0    1760 Good    4.82
#>  4       0    1152 Good    3.16
#>  5       0     967 Good    2.65
#>  6       0     629 Good    1.72
#>  7       0     461 Good    1.26
#>  8       1     991 Good    2.72
#>  9       0     758 Good    2.08
#> 10       0    1617 Good    4.43
#> # ℹ 676 more rows
fs1 <- flexsurvreg(Surv(recyrs, censrec)~group, data = bc, dist = 'weibull')

# Survival distribution representing patient w/ group = "Good"
good_model <- set_covariates(fs1, data.frame(group = 'Good'))
print(good_model)
#> Predicted survival from fitted model for specified cohort:
#>     * Cohort: 
#>           group
#>           <chr>
#>         1 Good 
#>     * Model: 
#>         Call:
#>         flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc, 
#>             dist = "weibull")
#>         
#>         Estimates: 
#>                      data mean  est      L95%     U95%     se       exp(est)  L95%   
#>         shape             NA     1.3797   1.2548   1.5170   0.0668       NA        NA
#>         scale             NA    11.4229   9.1818  14.2110   1.2728       NA        NA
#>         groupMedium   0.3338    -0.6136  -0.8623  -0.3649   0.1269   0.5414    0.4222
#>         groupPoor     0.3324    -1.2122  -1.4583  -0.9661   0.1256   0.2975    0.2326
#>                      U95%   
#>         shape             NA
#>         scale             NA
#>         groupMedium   0.6943
#>         groupPoor     0.3806
#>         
#>         N = 686,  Events: 299,  Censored: 387
#>         Total time at risk: 2113.425
#>         Log-likelihood = -811.9419, df = 4
#>         AIC = 1631.884
#> 

# Survival distribution with patient representing cohort w/
# equal mix of group = "Good" and group = "Poor"
good_and_poor_model <- set_covariates(fs1, data.frame(group = c('Good', 'Poor')))
print(good_and_poor_model)
#> Predicted survival from fitted model for specified cohort:
#>     * Cohort: 
#>           group
#>           <chr>
#>         1 Good 
#>         2 Poor 
#>     * Model: 
#>         Call:
#>         flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc, 
#>             dist = "weibull")
#>         
#>         Estimates: 
#>                      data mean  est      L95%     U95%     se       exp(est)  L95%   
#>         shape             NA     1.3797   1.2548   1.5170   0.0668       NA        NA
#>         scale             NA    11.4229   9.1818  14.2110   1.2728       NA        NA
#>         groupMedium   0.3338    -0.6136  -0.8623  -0.3649   0.1269   0.5414    0.4222
#>         groupPoor     0.3324    -1.2122  -1.4583  -0.9661   0.1256   0.2975    0.2326
#>                      U95%   
#>         shape             NA
#>         scale             NA
#>         groupMedium   0.6943
#>         groupPoor     0.3806
#>         
#>         N = 686,  Events: 299,  Censored: 387
#>         Total time at risk: 2113.425
#>         Log-likelihood = -811.9419, df = 4
#>         AIC = 1631.884
#> 
```
