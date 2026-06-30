# Define Kaplan-Meier Distribution

Define a survival distribution based on a table of Kaplan-Meier output
containing times and survival probabilities.

## Usage

``` r
define_surv_km(x, time_col = "time", surv_col = "survival")
```

## Arguments

- x:

  a data frame with columns for time and survival probability. By
  default, these columns are assumed to be named `time` and `survival`,
  but these can be configured by via the optional parameters.

- time_col:

  the name of the time column (defaults to `time`)

- surv_col:

  the name of the time column (defaults to `survival`)

## Value

a `surv_km` object.

## Examples

``` r
df <- data.frame(
     time = c(0, 1, 5, 10),
     survival = c(1, 0.9, 0.7, 0.5)
)
define_surv_km(df)
#> A Kaplan-Meier distribution:
#>    time survival
#>   <dbl>    <dbl>
#> 1     0      1  
#> 2     1      0.9
#> 3     5      0.7
#> 4    10      0.5
 
```
