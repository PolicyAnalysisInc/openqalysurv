# Apply Time-Dependent Hazard Ratio

Apply a time-dependent hazard ratio to a distribution where the hazard
ratio can vary across time points.

## Usage

``` r
apply_td_hr(dist, time, hr, log_hr = FALSE)
```

## Arguments

- dist:

  a survival distribution

- time:

  a numeric vector of time points corresponding to each hazard ratio

- hr:

  a numeric vector of hazard ratios (must be same length as time, or
  length 1)

- log_hr:

  optional argument (defaults to `FALSE`) to indicate that provided
  hazard ratios are on log scale

## Value

A `surv_td_ph` object if length(hr) \> 1, otherwise delegates to
apply_hr.

## Details

The `time` and `hr` vectors are parallel: `hr[i]` is the hazard ratio
that applies at `time[i]`. Duplicate (time, hr) pairs are allowed and
will be deduplicated. However, conflicting HRs for the same time point
(i.e., the same time with different HR values) will cause an error.

**Interaction with apply_af:**

\(1\) Composition order affects which HR applies at a given query time:

- `apply_af(apply_td_hr(dist, time, hr), af)`: HR intervals are in
  accelerated time

- `apply_td_hr(apply_af(dist, af), time, hr)`: HR intervals are in query
  time

\(2\) Acceleration factors affect HR granularity. When wrapping an AFT
distribution with af \< 1, each unit HR interval spans a larger portion
of the baseline survival curve, reducing the effective precision of
time-dependent effects.

## Examples

``` r

# Apply time-dependent HR: 0.5 for t in [0,5), 0.8 for t in [5,10)
td_dist <- apply_td_hr(
 define_surv_param("exp", rate = 0.1),
 time = 0:9,
 hr = c(rep(0.5, 5), rep(0.8, 5))
)
```
