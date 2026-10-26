# Define Survival Table Distribution

Another name for
[`define_surv_km()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_km.md):
a survival distribution from a table of times and survival
probabilities.

## Usage

``` r
define_surv_table(x, time_col = "time", surv_col = "survival")
```

## Arguments

- x:

  a data frame with columns for time and survival probability. By
  default, these columns are assumed to be named `time` and `survival`,
  but these can be configured by via the optional parameters.

- time_col:

  the name of the time column (defaults to `time`)

- surv_col:

  the name of the survival column (defaults to `survival`)

## Value

a `surv_km` object.
