# Extract Product-Limit Tables

Extracts the product-limit table from a survfit object for all strata.
Only `survfit` and unstratified `survfit.coxph` objects are supported.

## Usage

``` r
extract_strata(sf)
```

## Arguments

- sf:

  A survit object.

## Value

A tidy data.frame of the product-limit tables for all strata.
