#' Maximum Hazards
#'
#' Combine two or more survival distributions by taking the pointwise
#' maximum of their discrete hazards at each time step. The steps form a fixed
#' grid from time 0 with step `cycle_length / 10`, so the survival at a time
#' does not depend on which other times are requested.
#'
#' @name max_hazards
#' @rdname max_hazards
#' @export
#'
#' @param dist1 first survival distribution
#' @param dist2 second survival distribution
#' @param ... additional survival distributions
#' @param cycle_length the time step the distribution is used with (default
#'   1). `surv_prob` takes the maximum hazard over steps of
#'   `cycle_length / 10`, and `surv_quantile` searches in steps of
#'   `cycle_length`.
#'
#' @return A `surv_max_haz` object.
#'
#' @examples
#'
#' dist1 <- define_surv_param("exp", rate = 0.01)
#' dist2 <- define_surv_param("exp", rate = 0.02)
#' max_dist <- max_hazards(dist1, dist2)
#'
#' @tests
#'
#' dist1 <- define_surv_param("exp", rate = 0.01)
#' dist2 <- define_surv_param("exp", rate = 0.02)
#' dist3 <- define_surv_param("weibull", shape = 1.5, scale = 50)
#'
#' expect_equal(
#'  max_hazards(dist1, dist2),
#'  create_list_object(
#'      c('surv_max_haz', 'surv_combined', 'surv_dist'),
#'      dists = list(dist1, dist2),
#'      cycle_length = 1
#'  )
#' )
#'
#' # Custom cycle_length
#' expect_equal(
#'  max_hazards(dist1, dist2, cycle_length = 0.5)$cycle_length,
#'  0.5
#' )
#'
#' # Three distributions
#' expect_equal(
#'  length(max_hazards(dist1, dist2, dist3)$dists),
#'  3
#' )
#'
#' # Error: invalid dist
#' expect_error(
#'  max_hazards(dist1, 'foo'),
#'  'Error computing maximum hazards, invalid survival distribution provided.',
#'  fixed = TRUE
#' )
#'
#' # Error: invalid cycle_length
#' expect_error(
#'  max_hazards(dist1, dist2, cycle_length = -1),
#'  'Error computing maximum hazards, "cycle_length" must be a single positive numeric value.',
#'  fixed = TRUE
#' )
#'
#' expect_error(
#'  max_hazards(dist1, dist2, cycle_length = 0),
#'  'Error computing maximum hazards, "cycle_length" must be a single positive numeric value.',
#'  fixed = TRUE
#' )
#'
max_hazards <- function(dist1, dist2, ..., cycle_length = 1) {

    # Collect distributions (filter out named args like cycle_length)
    dots <- list(...)
    dists <- append(list(dist1, dist2), dots)

    # Validate distributions
    walk(dists, function(x) {
        valid <- is_surv_dist(x)
        if (!valid) {
            err <- get_and_populate_message('max_hazards_wrong_type_dist')
            stop(err, call. = show_call_error())
        }
    })

    # Validate cycle_length
    if (!is.numeric(cycle_length) || length(cycle_length) != 1 ||
        is.na(cycle_length) || cycle_length <= 0) {
        err <- get_and_populate_message('max_hazards_invalid_cycle_length')
        stop(err, call. = show_call_error())
    }

    create_list_object(
        c('surv_max_haz', 'surv_combined', 'surv_dist'),
        dists = dists,
        cycle_length = cycle_length
    )
}

#' @export
#'
#' @tests
#' # Two exponentials: max_hazards(exp(0.01), exp(0.02)) ~ exp(0.02)
#' dist1 <- define_surv_param("exp", rate = 0.01)
#' dist2 <- define_surv_param("exp", rate = 0.02)
#' mh <- max_hazards(dist1, dist2)
#' times <- seq(1, 50, by = 1)
#' expect_equal(
#'  surv_prob(mh, times),
#'  surv_prob(dist2, times),
#'  tolerance = 1e-6
#' )
#'
#' # S(0) = 1 and monotonically decreasing
#' expect_equal(surv_prob(mh, 0), 1)
#' s_vals <- surv_prob(mh, 0:100)
#' expect_true(all(diff(s_vals) <= 0))
#'
#' # The value at t does not depend on which other times are queried, and
#' # query times need not be sorted or unique
#' mh_cross <- max_hazards(define_surv_param("exp", rate = 0.1),
#'                         define_surv_param("weibull", shape = 2, scale = 10))
#' expect_equal(surv_prob(mh_cross, c(5, 10))[2], surv_prob(mh_cross, 10), tolerance = 1e-6)
#' expect_equal(
#'  surv_prob(mh_cross, c(10, 3, 10, 0)),
#'  c(surv_prob(mh_cross, 10), surv_prob(mh_cross, 3), surv_prob(mh_cross, 10), 1)
#' )
#'
#' # The hazards cross at t = 5, so S(10) = exp(-(0.5 + 0.75)) up to the grid step
#' expect_equal(surv_prob(mh_cross, 10), exp(-1.25), tolerance = 1e-3)
#'
surv_prob.surv_max_haz <- function(x, time, ...) {
    check_times(time, 'calculating survival probabilities', 'time')

    n_times <- length(time)
    if (n_times == 0) return(numeric(0))

    # Fixed grid from 0 with step cycle_length / 10, long enough to cover the
    # latest requested time. Results at a time depend only on the grid, not on
    # the other requested times, which may be in any order.
    step <- x$cycle_length / 10
    max_time <- max(time)
    n_steps <- ceiling(max_time / step)
    if (n_steps * step < max_time) n_steps <- n_steps + 1
    if (n_steps == 0) return(rep(1, n_times))
    grid <- (0:n_steps) * step

    # Cumulative hazard of each component on the grid
    cum_haz <- vapply(
        x$dists,
        function(d) -log(surv_prob(d, grid, ...)),
        numeric(length(grid))
    )
    cum_haz <- matrix(cum_haz, nrow = length(grid))

    # Hazard in each step: the largest component increment. A component whose
    # survival has reached 0 has infinite hazard from then on.
    incr <- apply(cum_haz, 2, diff)
    incr <- matrix(incr, nrow = length(grid) - 1)
    incr[is.nan(incr)] <- Inf
    max_incr <- do.call(pmax, as.data.frame(incr))
    grid_haz <- c(0, cumsum(max_incr))

    # Interpolate the cumulative hazard linearly within each step
    idx <- findInterval(time, grid, rightmost.closed = TRUE)
    w <- (time - grid[idx]) / step
    lower <- grid_haz[idx]
    upper <- grid_haz[idx + 1]
    haz <- ifelse(w == 0 | is.infinite(lower), lower, lower + w * (upper - lower))
    exp(-haz)
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param("exp", rate = 0.01)
#' dist2 <- define_surv_param("exp", rate = 0.02)
#' mh <- max_hazards(dist1, dist2)
#'
#' # Roundtrip
#' probs <- c(0.9, 0.5)
#' times <- surv_quantile(mh, probs)
#' expect_equal(surv_prob(mh, times), probs, tolerance = 0.01)
#'
#' # Edge cases
#' expect_equal(surv_quantile(mh, 1), 0)
#' expect_equal(surv_quantile(mh, 0), Inf)
#'
surv_quantile.surv_max_haz <- function(x, probs, ...) {
    check_probs(probs)
    vapply(probs, function(p) {
        if (p == 1) return(0)
        if (p == 0) return(Inf)

        # Upper bound: min_i Q_i(p) (S_max <= min_i S_i)
        component_quantiles <- map_dbl(x$dists, function(d) surv_quantile(d, p, ...))
        upper <- min(component_quantiles)
        if (is.infinite(upper)) return(Inf)

        # Build grid from 0 to upper at cycle_length resolution
        grid <- seq(0, upper, by = x$cycle_length)
        if (grid[length(grid)] < upper) grid <- c(grid, upper)

        # Compute S_max on grid in a single forward pass
        s_max <- surv_prob(x, grid, ...)

        # Find first grid point where S_max(t) <= p
        idx <- which(s_max <= p)[1]
        if (is.na(idx)) return(upper)
        grid[idx]
    }, numeric(1))
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param("exp", rate = 0.01)
#' dist2 <- define_surv_param("exp", rate = 0.02)
#' expect_output(
#'  print(max_hazards(dist1, dist2)),
#'  'A survival distribution taking the maximum hazard of:',
#'  fixed = TRUE
#' )
#'
print.surv_max_haz <- function(x, ...) {
    items_lines <- map_chr(seq_along(x$dists), function(i) {
        dist_output <- to_list_item_output(x$dists[[i]])
        glue('    * Distribution {i}: {dist_output}')
    })
    output <- paste0(
        c(
            'A survival distribution taking the maximum hazard of:',
            items_lines,
            glue('    * Cycle length: {x$cycle_length}')
        ),
        collapse = '\n'
    )
    cat(output)
}
