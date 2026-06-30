# Define Life-Table Distribution

Define a survival distribution based on a life-table containing
mortality rates by age and gender.

## Usage

``` r
define_surv_lifetable(
  x,
  start_age,
  percent_male,
  output_unit = "years",
  age_col = "age",
  male_col = "male",
  female_col = "female",
  dpy = get_dpy(),
  percent_female = NULL
)
```

## Arguments

- x:

  a data frame with columns for the starting age in each age band, the
  conditional probability of death at each age for men, and the
  conditional probability of death at each for women. Default column
  names are "age", "male", and "female", but can be set via optional
  arguments.

- start_age:

  starting age of the population.

- percent_male:

  percent of population that is male. Calculated as 1 - percent_female
  if not provided.

- output_unit:

  optional arguemnt for time unit resulting survival distribution will
  be defined in terms of. Valid options are
  `c("days", "weeks", "months", "years")`. Defaults to `"years"`.

- age_col:

  optional argument to change name of the `age` column accepted by the
  first argument.

- male_col:

  optional argument to change name of the `male` column accepted by the
  first argument.

- female_col:

  optional argument to change name of the `female` column accepted by
  the first argument.

- dpy:

  optional argument specifying the number of days per year used in time
  unit conversion. This argument will be populated automatically in an
  openqaly model.

- percent_female:

  percent of population that is female. Calculated as 1 - percent_male
  if not provided.

## Value

a `surv_lifetable` object.

## Examples

``` r
 x <- data.frame(
     age = c(0, 1, 2, 3),
     male = c(0.011, 0.005, 0.003, 0.002),
     female = c(0.010, 0.005, 0.004, 0.002)
 )
 define_surv_lifetable(x, 1, 0.45)
#> A life-table survival distribution (45.0% male, 55.0% female):
#>     age  male female
#>   <dbl> <dbl>  <dbl>
#> 1     1 0.005  0.005
#> 2     2 0.003  0.004
#> 3     3 0.002  0.002
 
```
