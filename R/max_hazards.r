#' Maximum Hazards
#'
#' Combine two or more survival distributions by taking the pointwise
#' maximum of their discrete hazards at each time step.
#'
#' @name max_hazards
#' @rdname max_hazards
#' @export
#'
#' @param dist1 first survival distribution
#' @param dist2 second survival distribution
#' @param ... additional survival distributions
#' @param cycle_length step size used by `surv_quantile` for grid search
#'   (default 1). Has no effect on `surv_prob`.
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
surv_prob.surv_max_haz <- function(x, time, ...) {
    check_times(time, 'calculating survival probabilities', 'time')

    n_times <- length(time)
    if (n_times == 0) return(numeric(0))

    # Sort time and track original order
    ord <- order(time)
    sorted_time <- time[ord]

    # Compute S_k(t) for each component at all query times
    all_surv <- map(x$dists, function(d) surv_prob(d, sorted_time, ...))

    # Initialize result
    result <- numeric(n_times)
    result[1] <- 1  # S(t_0) if t_0 = 0, handle below

    if (sorted_time[1] == 0) {
        result[1] <- 1
        start_idx <- 2
    } else {
        # Need to compute from 0 to first time point
        # Get S_k at a "previous" time of 0
        prev_surv <- map_dbl(x$dists, function(d) surv_prob(d, 0, ...))
        max_haz <- 0
        for (k in seq_along(x$dists)) {
            h_k <- -log(all_surv[[k]][1] / prev_surv[k])
            max_haz <- max(max_haz, h_k)
        }
        result[1] <- exp(-max_haz)
        start_idx <- 2
    }

    if (n_times >= start_idx) {
        for (i in start_idx:n_times) {
            # Discrete hazard for each component in step [t_{i-1}, t_i]
            max_haz <- 0
            for (k in seq_along(x$dists)) {
                s_prev <- all_surv[[k]][i - 1]
                s_curr <- all_surv[[k]][i]
                if (s_prev > 0) {
                    h_k <- -log(s_curr / s_prev)
                } else {
                    h_k <- Inf
                }
                max_haz <- max(max_haz, h_k)
            }
            result[i] <- result[i - 1] * exp(-max_haz)
        }
    }

    # Restore original order
    final <- numeric(n_times)
    final[ord] <- result
    final
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
