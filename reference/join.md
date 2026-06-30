# Join Distributions

Join two or more distributions together at the specified cut points.

## Usage

``` r
join(dist1, cut1, dist2, ...)
```

## Arguments

- dist1:

  survival distribution to use from time `0` to `cut`

- cut1:

  cut point between `dist1` and `dist2`

- dist2:

  survival distribution to use from `cut`

- ...:

  Additional cutpoints and distributions

## Value

A `surv_join` object

## Details

`join` ensures the overall survival curve is continuous at each cut
point by rescaling the next distribution's survival values. For a single
cut point \\c\\ joining distributions \\S_1\\ and \\S_2\\, the joined
survival function is:

\$\$S(t) = \begin{cases} S_1(t) & t \le c \\ \frac{S_1(c)}{S_2(c)} \cdot
S_2(t) & t \> c \end{cases}\$\$

Importantly, the second distribution is evaluated at absolute time
\\t\\, not at \\t - c\\. This means `join` only rescales survival
probabilities for continuity — it does not shift or remap the timescale
of any distribution.

All distributions passed to `join` must therefore share the same time
origin (typically time zero = study baseline). If a distribution was
estimated from data that begins at the cut point — so that its own time
zero corresponds to the cut point — it must first be shifted with
[`apply_shift`](https://policyanalysisinc.github.io/openqalysurv/reference/apply_shift.md)
to realign the timescales:


    # dist_from_cutpoint was fit to data starting at t = c
    # (its time zero = the cut point)
    shifted <- apply_shift(dist_from_cutpoint, c)
    joined  <- join(km, c, shifted)

Without the shift, `join` would evaluate the distribution at absolute
time \\t\\ instead of the intended \\t - c\\, producing incorrect
survival estimates beyond the cut point.

## Examples

``` r

dist1 <- define_surv_param("exp", rate = 0.05)
dist2 <- define_surv_param("gompertz", rate = .5, shape = 1)
dist3 <- define_surv_param("exp", rate = 0.25)
join_dist <- join(dist1, 20, dist2)
join_dist2 <- join(dist1, 20, dist2, 50, dist3)
```
