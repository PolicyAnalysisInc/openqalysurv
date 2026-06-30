# Extract Product-Limit Table for Stratum

Extracts the product-limit table from a survfit object for a given
stratum. Only
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
and unstratified
[`survival::survfit.coxph()`](https://rdrr.io/pkg/survival/man/survfit.coxph.html)
objects are supported.

## Usage

``` r
extract_stratum(sf, index)
```

## Arguments

- sf:

  A survit object.

- index:

  The index number of the strata to extract.

## Value

A data frame of the product-limit table for the given stratum.
