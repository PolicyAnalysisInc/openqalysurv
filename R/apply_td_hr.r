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
#' The `time` and `hr` vectors are parallel: `hr[i]` is the hazard ratio that
#' applies at `time[i]`. Duplicate (time, hr) pairs are allowed and will be
#' deduplicated. However, conflicting HRs for the same time point (i.e., the same
#' time with different HR values) will cause an error.
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
#' # Apply time-dependent HR: 0.5 for t in [0,5), 0.8 for t in [5,10)
#' td_dist <- apply_td_hr(
#'  define_surv_param("exp", rate = 0.1),
#'  time = 0:9,
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
#' # Error when time exceeds max stored time
#' expect_error(
#'  surv_prob(td_dist2, 5),
#'  'Error calculating survival probabilities, query time exceeds maximum stored time.',
#'  fixed = TRUE
#' )
#'
#' # CORE VALIDATION: td_hr should match join + apply_hr for piecewise constant HR
#' # HR=0.5 for [0,2), HR=0.8 for [2,4)
#' td_dist3 <- apply_td_hr(dist1, 0:3, c(0.5, 0.5, 0.8, 0.8))
#' join_dist <- join(apply_hr(dist1, 0.5), 2, apply_hr(dist1, 0.8))
#' test_times <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
#' expect_equal(
#'  surv_prob(td_dist3, test_times),
#'  surv_prob(join_dist, test_times),
#'  tolerance = 1e-10
#' )
#'
#' # Multiple HR transitions: HR=0.5 for [0,2), HR=1.0 for [2,4), HR=1.5 for [4,6)
#' td_dist4 <- apply_td_hr(dist1, 0:5, c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5))
#' join_dist2 <- join(apply_hr(dist1, 0.5), 2, apply_hr(dist1, 1.0), 4, apply_hr(dist1, 1.5))
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
#' join_weib <- join(apply_hr(dist_weib, 0.5), 2, apply_hr(dist_weib, 0.8))
#' expect_equal(
#'  surv_prob(td_weib, test_times),
#'  surv_prob(join_weib, test_times),
#'  tolerance = 1e-10
#' )
#'
#' # Composite baseline: apply_hr on top of existing surv_ph
#' dist_ph <- apply_hr(dist1, 0.7)
#' td_composite <- apply_td_hr(dist_ph, 0:3, c(0.5, 0.5, 0.8, 0.8))
#' join_composite <- join(apply_hr(dist_ph, 0.5), 2, apply_hr(dist_ph, 0.8))
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
#' join_sparse <- join(apply_hr(dist1, 0.5), 2, apply_hr(dist1, 0.8), 4, apply_hr(dist1, 1.0))
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
#' join_dense <- join(apply_hr(dist1, 0.5), 5, apply_hr(dist1, 0.8), 10, apply_hr(dist1, 1.0))
#' expect_equal(
#'  surv_prob(td_dense, seq(0, 10, by = 0.5)),
#'  surv_prob(join_dense, seq(0, 10, by = 0.5)),
#'  tolerance = 1e-10
#' )
#'
#' # Test: Weibull with non-aligned times
#' dist_weib2 <- define_surv_param("weibull", shape = 1.5, scale = 10)
#' td_weib2 <- apply_td_hr(dist_weib2, c(0, 2, 5), c(0.5, 0.8, 1.2))
#' join_weib2 <- join(apply_hr(dist_weib2, 0.5), 2, apply_hr(dist_weib2, 0.8), 5, apply_hr(dist_weib2, 1.2))
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

    # Check that query times don't exceed max stored time
    max_query <- max(time)
    max_stored <- max(x$time)
    if (max_query > max_stored) {
        err <- get_and_populate_message('apply_td_hr_time_exceeds_stored')
        stop(err, call. = show_call_error())
    }

    # Get baseline log-survival at stored time points
    log_surv_stored <- log(surv_prob(x$dist, x$time, ...))

    # Get baseline log-survival at query times
    log_surv_query <- log(surv_prob(x$dist, time, ...))

    # Call Rcpp function for efficient computation
    compute_td_surv_prob(time, x$time, x$hr, log_surv_stored, log_surv_query)
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
