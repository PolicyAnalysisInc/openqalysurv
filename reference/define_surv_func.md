# Define survival distribution based on a function

Define survival distribution based on a function

## Usage

``` r
define_surv_func(f, ...)
```

## Arguments

- f:

  a function that takes a vector of times and returns a vector of
  corresponding survival probabilities

- ...:

  additional arguments to be passed to f

## Value

a `surv_function` object.
