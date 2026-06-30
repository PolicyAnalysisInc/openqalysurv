# Apply Fixed Shift in Time

Apply a fixed shift in time to move the hazards of a survival
distribution forwards or backwards.

## Usage

``` r
apply_shift(dist, shift)
```

## Arguments

- dist:

  a survival distribution

- shift:

  amount of time to shift the distribution. A positive value delays
  hazards (shifts the survival curve to the right), and a negative value
  accelerates hazards (shifts the survival curve to the left).

## Value

A `surv_shift` object.

## Details

The shifted survival function is defined as:

\$\$S\_{\mathrm{shifted}}(t) = \begin{cases} 1 & t \le \delta \\ S(t -
\delta) & t \> \delta \end{cases}\$\$

where \\\delta\\ is the shift amount and \\S\\ is the baseline survival
function.

A positive shift delays all hazards by \\\delta\\ time units (the curve
moves to the right); a negative shift accelerates them (the curve moves
to the left).

A common use of `apply_shift` is to realign the timescale of a
distribution before passing it to
[`join`](https://policyanalysisinc.github.io/openqalysurv/reference/join.md).
When a distribution is estimated from data that begins at a cut point
\\c\\ — so that its own time zero corresponds to \\c\\ — applying
`apply_shift(dist, c)` maps its internal time to the absolute study
timescale, which `join` requires.

## Examples

``` r

shifted_dist <- apply_shift(
 define_surv_param("exp", rate = 0.25),
 3
)
```
