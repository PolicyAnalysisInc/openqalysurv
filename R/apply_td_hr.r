#' Apply Time-Dependent Hazard Ratio
#'
#' Apply a time-dependent hazard ratio to a distribution where the hazard ratio
#' can vary across time points.
#'
#' @name apply_td_hr
#' @rdname apply_td_hr
#' @export
#'
#' @param dist a survival distribution
#' @param time a numeric vector of time points corresponding to each hazard ratio
#' @param hr a numeric vector of hazard ratios (must be same length as time, or length 1)
#' @param log_hr optional argument (defaults to `FALSE`) to indicate that provided
#'   hazard ratios are on log scale
#'
#' @return A `surv_td_ph` object if length(hr) > 1, otherwise delegates to apply_hr.
#'
#' @details
#' The `time` and `hr` vectors are parallel: `hr[i]` applies on the interval
#' `(time[i - 1], time[i]]`, with `time[0] = 0`, and the last hazard ratio
#' continues after the last time point. With `time` set to the end of each
#' cycle (for example `cycle`), cycle k is given `hr[k]`. A hazard ratio whose
#' interval has zero length (for example `hr[1]` when `time[1] = 0`) covers no
#' time. Duplicate (time, hr) pairs are allowed and will be deduplicated.
#' However, conflicting HRs for the same time point (i.e., the same time with
#' different HR values) will cause an error.
#'
#' **Interaction with apply_af:**
#'
#' (1) Composition order affects which HR applies at a given query time:
#' \itemize{
#'   \item `apply_af(apply_td_hr(dist, time, hr), af)`: HR intervals are in accelerated time
#'   \item `apply_td_hr(apply_af(dist, af), time, hr)`: HR intervals are in query time
#' }
#'
#' (2) Acceleration factors affect HR granularity. When wrapping an AFT
#' distribution with af < 1, each unit HR interval spans a larger portion of
#' the baseline survival curve, reducing the effective precision of
#' time-dependent effects.
#'
#' @examples
#'
#' # Apply time-dependent HR: 0.5 in cycles 1-5 (t in (0, 5]), 0.8 after
#' td_dist <- apply_td_hr(
#'  define_surv_param("exp", rate = 0.1),
#'  time = 1:10,
#'  hr = c(rep(0.5, 5), rep(0.8, 5))
#' )
#'
#' @tests
#' dist1 <- define_surv_param("exp", rate = 0.1)
#'
#' # Delegation to apply_hr when length(hr) == 1
#' expect_equal(
#'  apply_td_hr(dist1, 0, 0.5),
#'  apply_hr(dist1, 0.5)
#' )
#'
#' # Log scale delegation
#' expect_equal(
#'  apply_td_hr(dist1, 0, log(0.5), TRUE),
#'  apply_hr(dist1, 0.5)
#' )
#'
#' # Basic structure for td_hr
#' td_dist <- apply_td_hr(dist1, 0:1, c(0.5, 0.8))
#' expect_true(inherits(td_dist, 'surv_td_ph'))
#' expect_true(inherits(td_dist, 'surv_dist'))
#' expect_equal(td_dist$hr, c(0.5, 0.8))
#' expect_equal(td_dist$time, c(0, 1))
#'
#' # Duplicate (time, hr) pairs are deduplicated
#' td_dup <- apply_td_hr(dist1, c(0, 1, 0, 1), c(0.5, 0.8, 0.5, 0.8))
#' expect_equal(td_dup$time, c(0, 1))
#' expect_equal(td_dup$hr, c(0.5, 0.8))
#'
#' # Error conditions
#' expect_error(
#'  apply_td_hr('foo', 0:1, c(0.5, 0.8)),
#'  'Error applying time-dependent hazard ratio, invalid survival distribution provided.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_td_hr(dist1, 0:1, 'foo'),
#'  'Error applying time-dependent hazard ratio, "hr" must be numeric.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_td_hr(dist1, 0:1, c(0.5, NA_real_)),
#'  'Error applying time-dependent hazard ratio, "hr" cannot contain NA values.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_td_hr(dist1, 0:1, c(0.5, -0.2)),
#'  'Error applying time-dependent hazard ratio, "hr" cannot contain negative values.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_td_hr(dist1, c(0, 1, 2), c(0.5, 0.8)),
#'  'Error applying time-dependent hazard ratio, "time" and "hr" must have the same length.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_td_hr(dist1, c(0, 1, 0), c(0.5, 0.8, 0.7)),
#'  'Error applying time-dependent hazard ratio, conflicting hazard ratios found for the same time value.',
#'  fixed = TRUE
#' )
apply_td_hr <- function(dist, time, hr, log_hr = FALSE) {

    # Check that dist is a valid type
    if (!is_surv_dist(dist)) {
        err <- get_and_populate_message('apply_td_hr_wrong_type_dist')
        stop(err, call. = show_call_error())
    }

    # Check that hr is numeric
    if (!any(c('integer', 'numeric') %in% class(hr))) {
        err <- get_and_populate_message('apply_td_hr_wrong_type_hr')
        stop(err, call. = show_call_error())
    }

    # If log_hr is specified then exponentiate it
    if (log_hr) {
        hr <- exp(hr)
    }

    # Check that hr doesn't contain missing values
    if (any(is.na(hr))) {
        err <- get_and_populate_message('apply_td_hr_missing_hr')
        stop(err, call. = show_call_error())
    }

    # Check that hr >= 0 for all values
    if (any(hr < 0)) {
        err <- get_and_populate_message('apply_td_hr_invalid_hr')
        stop(err, call. = show_call_error())
    }

    # If length(hr) == 1, delegate to apply_hr
    if (length(hr) == 1) {
        return(apply_hr(dist, hr, log_hr = FALSE))
    }

    # Validate time
    check_times(time, 'applying time-dependent hazard ratio', 'time')

    # Validate lengths match
    if (length(time) != length(hr)) {
        err <- get_and_populate_message('apply_td_hr_length_mismatch')
        stop(err, call. = show_call_error())
    }

    # Check for conflicting HRs at same time
    time_hr_df <- data.frame(time = time, hr = hr)
    unique_pairs <- unique(time_hr_df)
    if (anyDuplicated(unique_pairs$time)) {
        err <- get_and_populate_message('apply_td_hr_conflicting_hr')
        stop(err, call. = show_call_error())
    }

    # Sort by time
    unique_pairs <- unique_pairs[order(unique_pairs$time), ]

    create_list_object(
        c('surv_td_ph', 'surv_dist'),
        dist = dist,
        time = unique_pairs$time,
        hr = unique_pairs$hr
    )

}

#' @export
#'
#' @tests
#'
#' dist1 <- define_surv_param("exp", rate = 0.1)
#'
#' # Constant HR should match apply_hr
#' td_dist <- apply_td_hr(dist1, 0:99, rep(0.5, 100))
#' ph_dist <- apply_hr(dist1, 0.5)
#' expect_equal(
#'  surv_prob(td_dist, 0:99),
#'  surv_prob(ph_dist, 0:99),
#'  tolerance = 1e-10
#' )
#'
#' # Identity: HR = 1 should match baseline
#' td_dist_id <- apply_td_hr(dist1, 0:49, rep(1, 50))
#' expect_equal(
#'  surv_prob(td_dist_id, 0:49),
#'  surv_prob(dist1, 0:49),
#'  tolerance = 1e-10
#' )
#'
#' # Time = 0 should always be 1
#' td_dist2 <- apply_td_hr(dist1, 0:2, c(0.5, 0.8, 1.2))
#' expect_equal(surv_prob(td_dist2, 0), 1)
#'
#' # Cycle example: hr[k] applies on (k - 1, k], so cycle k gets hr[k]
#' cyc_hr <- 0.3 + 0.01 * (0:9)
#' cyc <- apply_td_hr(dist1, 1:10, cyc_hr)
#' expect_equal(
#'  surv_prob(cyc, 1:10),
#'  exp(-0.1 * cumsum(cyc_hr)),
#'  tolerance = 1e-10
#' )
#'
#' # Query times need not be sorted or unique
#' expect_equal(
#'  surv_prob(cyc, c(7, 2, 7, 0.5)),
#'  surv_prob(cyc, sort(unique(c(7, 2, 0.5))))[c(3, 2, 3, 1)],
#'  tolerance = 1e-12
#' )
#'
#' # Quantile round-trip on the cycle example
#' expect_equal(
#'  surv_quantile(cyc, surv_prob(cyc, c(2.5, 6.5))),
#'  c(2.5, 6.5),
#'  tolerance = 1e-6
#' )
#'
#' # A single HR matches apply_hr at any time
#' expect_equal(
#'  surv_prob(apply_td_hr(dist1, 100, 2), c(1, 3, 50, 120)),
#'  surv_prob(apply_hr(dist1, 2), c(1, 3, 50, 120)),
#'  tolerance = 1e-10
#' )
#'
#' # time = c(5, 10), hr = c(2, 0.5): hr 2 on (0, 5], then 0.5 after 5 (continues past 10)
#' expect_equal(
#'  surv_prob(apply_td_hr(dist1, c(5, 10), c(2, 0.5)), c(1, 5, 7, 10, 12)),
#'  surv_prob(join(apply_hr(dist1, 2), 5, apply_hr(dist1, 0.5)), c(1, 5, 7, 10, 12)),
#'  tolerance = 1e-10
#' )
#'
#' # Past the last stored time the last HR continues indefinitely
#' # hr[i] applies on (time[i-1], time[i]] with time[0] = 0, so hr[1] covers no time here
#' join_past <- join(apply_hr(dist1, 0.8), 1, apply_hr(dist1, 1.2))
#' test_past <- c(2, 3, 5, 10)
#' expect_equal(
#'  surv_prob(td_dist2, test_past),
#'  surv_prob(join_past, test_past),
#'  tolerance = 1e-10
#' )
#' expect_equal(
#'  surv_prob(apply_td_hr(dist1, 0, 0.5), 3),
#'  surv_prob(apply_hr(dist1, 0.5), 3),
#'  tolerance = 1e-10
#' )
#'
#' # CORE VALIDATION: td_hr should match join + apply_hr for piecewise constant HR
#' # HR=0.5 for (0,1], HR=0.8 after 1 (hr[2] covers (0,1], hr[3] covers (1,2])
#' td_dist3 <- apply_td_hr(dist1, 0:3, c(0.5, 0.5, 0.8, 0.8))
#' join_dist <- join(apply_hr(dist1, 0.5), 1, apply_hr(dist1, 0.8))
#' test_times <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
#' expect_equal(
#'  surv_prob(td_dist3, test_times),
#'  surv_prob(join_dist, test_times),
#'  tolerance = 1e-10
#' )
#'
#' # Multiple HR transitions: HR=0.5 for (0,1], HR=1.0 for (1,3], HR=1.5 after 3
#' td_dist4 <- apply_td_hr(dist1, 0:5, c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5))
#' join_dist2 <- join(apply_hr(dist1, 0.5), 1, apply_hr(dist1, 1.0), 3, apply_hr(dist1, 1.5))
#' test_times2 <- c(0, 1, 2, 3, 4, 5)
#' expect_equal(
#'  surv_prob(td_dist4, test_times2),
#'  surv_prob(join_dist2, test_times2),
#'  tolerance = 1e-10
#' )
#'
#' # Non-exponential: Weibull distribution
#' dist_weib <- define_surv_param("weibull", shape = 1.5, scale = 10)
#' td_weib <- apply_td_hr(dist_weib, 0:3, c(0.5, 0.5, 0.8, 0.8))
#' join_weib <- join(apply_hr(dist_weib, 0.5), 1, apply_hr(dist_weib, 0.8))
#' expect_equal(
#'  surv_prob(td_weib, test_times),
#'  surv_prob(join_weib, test_times),
#'  tolerance = 1e-10
#' )
#'
#' # Composite baseline: apply_hr on top of existing surv_ph
#' dist_ph <- apply_hr(dist1, 0.7)
#' td_composite <- apply_td_hr(dist_ph, 0:3, c(0.5, 0.5, 0.8, 0.8))
#' join_composite <- join(apply_hr(dist_ph, 0.5), 1, apply_hr(dist_ph, 0.8))
#' expect_equal(
#'  surv_prob(td_composite, test_times),
#'  surv_prob(join_composite, test_times),
#'  tolerance = 1e-10
#' )
#'
#' # NON-ALIGNED QUERY TIMES: query times between stored HR time points
#'
#' # Test: sparse stored times, queries fall between them
#' td_sparse <- apply_td_hr(dist1, c(0, 2, 4), c(0.5, 0.8, 1.0))
#' join_sparse <- join(apply_hr(dist1, 0.8), 2, apply_hr(dist1, 1.0))
#' expect_equal(
#'  surv_prob(td_sparse, c(0, 1, 2, 3, 4)),
#'  surv_prob(join_sparse, c(0, 1, 2, 3, 4)),
#'  tolerance = 1e-10
#' )
#'
#' # Test: fractional query times
#' expect_equal(
#'  surv_prob(td_sparse, c(0, 0.5, 1.5, 2, 2.5, 3.5, 4)),
#'  surv_prob(join_sparse, c(0, 0.5, 1.5, 2, 2.5, 3.5, 4)),
#'  tolerance = 1e-10
#' )
#'
#' # Test: sparse stored, dense queries
#' td_dense <- apply_td_hr(dist1, c(0, 5, 10), c(0.5, 0.8, 1.0))
#' join_dense <- join(apply_hr(dist1, 0.8), 5, apply_hr(dist1, 1.0))
#' expect_equal(
#'  surv_prob(td_dense, seq(0, 10, by = 0.5)),
#'  surv_prob(join_dense, seq(0, 10, by = 0.5)),
#'  tolerance = 1e-10
#' )
#'
#' # Test: Weibull with non-aligned times
#' dist_weib2 <- define_surv_param("weibull", shape = 1.5, scale = 10)
#' td_weib2 <- apply_td_hr(dist_weib2, c(0, 2, 5), c(0.5, 0.8, 1.2))
#' join_weib2 <- join(apply_hr(dist_weib2, 0.8), 2, apply_hr(dist_weib2, 1.2))
#' expect_equal(
#'  surv_prob(td_weib2, c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)),
#'  surv_prob(join_weib2, c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)),
#'  tolerance = 1e-10
#' )
#'
#' # Test: single query point mid-interval
#' expect_equal(surv_prob(td_dense, 3), surv_prob(join_dense, 3), tolerance = 1e-10)
#' expect_equal(surv_prob(td_dense, 7), surv_prob(join_dense, 7), tolerance = 1e-10)
#'
surv_prob.surv_td_ph <- function(x, time, ...) {
    check_times(time, 'calculating survival probabilities', 'time')

    # Interval starts: hr[i] applies from the previous stored time
    iv <- td_hr_intervals(x)

    # Get baseline log-survival at the interval starts
    log_surv_stored <- log(surv_prob(x$dist, iv$start, ...))

    # Get baseline log-survival at query times
    log_surv_query <- log(surv_prob(x$dist, time, ...))

    # Call Rcpp function for efficient computation (each query time is
    # evaluated independently, so query times need not be sorted)
    compute_td_surv_prob(time, iv$start, iv$hr, log_surv_stored, log_surv_query)
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.1)
#'
#' # Roundtrip for td_ph (use high probs that resolve within stored range)
#' td_dist <- apply_td_hr(dist1, 0:99, rep(c(0.5, 0.8), each = 50))
#' probs <- c(0.9, 0.5, 0.2)
#' times <- surv_quantile(td_dist, probs)
#' expect_equal(surv_prob(td_dist, times), probs, tolerance = 1e-4)
#'
#' # Roundtrip across HR change points, including past the last stored time
#' dist_cp <- define_surv_param('weibull', shape = 1.5, scale = 10)
#' td_cp <- apply_td_hr(dist_cp, c(0, 2, 5), c(0.5, 0.8, 1.2))
#' query_cp <- c(1, 3, 8, 20)
#' probs_cp <- surv_prob(td_cp, query_cp)
#' expect_equal(surv_quantile(td_cp, probs_cp), query_cp, tolerance = 1e-6)
#' expect_equal(surv_prob(td_cp, surv_quantile(td_cp, probs_cp)), probs_cp, tolerance = 1e-8)
#'
#' # Edge cases
#' expect_equal(surv_quantile(td_dist, 1), 0)
#' expect_equal(surv_quantile(td_dist, 0), Inf)
#'
surv_quantile.surv_td_ph <- function(x, probs, ...) {
    check_probs(probs)

    # Interval starts: hr[k] applies from start[k] to start[k + 1]
    iv <- td_hr_intervals(x)
    start <- iv$start
    hr <- iv$hr
    n <- length(start)

    # Compute S_td at each interval start
    stored_surv <- surv_prob(x, start, ...)

    vapply(probs, function(p) {
        if (p == 1) return(0)
        if (p == 0) return(Inf)

        # If p >= S_td(start_1), it's before the first interval start
        if (p >= stored_surv[1]) {
            if (hr[1] == 0) return(Inf)
            p_base <- p^(1/hr[1])
            return(surv_quantile(x$dist, p_base, ...))
        }

        # Find the interval
        for (k in seq_len(n - 1)) {
            if (p >= stored_surv[k + 1]) {
                # p is in interval [start_k, start_{k+1}), where hr[k] applies
                if (hr[k] == 0) next
                # S_td(t) = S_td(t_k) * (S_base(t)/S_base(t_k))^hr_k
                # Invert: S_base(t) = S_base(t_k) * (p/S_td(t_k))^(1/hr_k)
                s_base_at_tk <- surv_prob(x$dist, start[k], ...)
                p_base <- s_base_at_tk * (p / stored_surv[k])^(1/hr[k])
                return(surv_quantile(x$dist, p_base, ...))
            }
        }

        # Extrapolate beyond the last interval start with the last HR
        if (hr[n] == 0) return(Inf)
        s_base_at_tn <- surv_prob(x$dist, start[n], ...)
        p_base <- s_base_at_tn * (p / stored_surv[n])^(1/hr[n])
        surv_quantile(x$dist, p_base, ...)
    }, numeric(1))
}

#' @export
#'
#' @tests
#' dist1 <- apply_td_hr(define_surv_param('exp', rate = 0.1), 0:2, c(0.5, 0.8, 1.2))
#' expect_output(print(dist1), 'time-dependent proportional hazards')
#'
print.surv_td_ph <- function(x, ...) {
    bl_dist_output <- to_list_item_output(x$dist)
    hr_min <- min(x$hr)
    hr_max <- max(x$hr)
    n_points <- length(x$hr)
    time_range <- paste0('[', min(x$time), ', ', max(x$time), ']')
    output <- paste0(
        c(
            'A time-dependent proportional hazards survival distribution:',
            glue('    * Number of time points: {n_points}'),
            glue('    * Time range: {time_range}'),
            glue('    * Hazard ratio range: [{hr_min}, {hr_max}]'),
            glue('    * Baseline Distribution: {bl_dist_output}')
        ),
        collapse = '\n'
    )
    cat(output)
}

# Intervals of a time-dependent hazard ratio: hr[i] applies on
# (time[i - 1], time[i]] with time[0] = 0, so interval i starts at time[i - 1].
# Of several intervals sharing a start, all but the last have zero length and
# are dropped; the last interval always remains.
td_hr_intervals <- function(x) {
    n <- length(x$time)
    start <- c(0, x$time[-n])
    keep <- !duplicated(start, fromLast = TRUE)
    list(start = start[keep], hr = x$hr[keep])
}
