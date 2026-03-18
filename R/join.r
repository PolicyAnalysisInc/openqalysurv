#' Join Distributions
#' 
#' Join two or more distributions together at the specified cut points.
#' 
#' @name join
#' @rdname join
#' @export
#' 
#' @param dist1 survival distribution to use from time `0` to `cut`
#' @param cut1 cut point between `dist1` and `dist2`
#' @param dist2 survival distribution to use from `cut`
#' @param ... Additional cutpoints and distributions
#'   
#' @return A `surv_join` object
#'
#' @details
#' \code{join} ensures the overall survival curve is continuous at each cut
#' point by rescaling the next distribution's survival values. For a single
#' cut point \eqn{c} joining distributions \eqn{S_1} and \eqn{S_2}, the
#' joined survival function is:
#'
#' \deqn{S(t) = \begin{cases} S_1(t) & t \le c \\ \frac{S_1(c)}{S_2(c)}
#' \cdot S_2(t) & t > c \end{cases}}{S(t) = S_1(t) for t <= c;
#' S(t) = [S_1(c) / S_2(c)] * S_2(t) for t > c}
#'
#' Importantly, the second distribution is evaluated at absolute time
#' \eqn{t}, not at \eqn{t - c}. This means \code{join} only rescales
#' survival probabilities for continuity — it does not shift or remap the
#' timescale of any distribution.
#'
#' All distributions passed to \code{join} must therefore share the same
#' time origin (typically time zero = study baseline). If a distribution was
#' estimated from data that begins at the cut point — so that its own
#' time zero corresponds to the cut point — it must first be shifted with
#' \code{\link{apply_shift}} to realign the timescales:
#'
#' \preformatted{
#' # dist_from_cutpoint was fit to data starting at t = c
#' # (its time zero = the cut point)
#' shifted <- apply_shift(dist_from_cutpoint, c)
#' joined  <- join(km, c, shifted)
#' }
#'
#' Without the shift, \code{join} would evaluate the distribution at
#' absolute time \eqn{t} instead of the intended \eqn{t - c}, producing
#' incorrect survival estimates beyond the cut point.
#'
#' @examples
#'
#' dist1 <- define_surv_param("exp", rate = 0.05)
#' dist2 <- define_surv_param("gompertz", rate = .5, shape = 1)
#' dist3 <- define_surv_param("exp", rate = 0.25)
#' join_dist <- join(dist1, 20, dist2)
#' join_dist2 <- join(dist1, 20, dist2, 50, dist3)
#' 
#' @tests
#'
#' dist1 <- define_surv_param("exp", rate = 0.05)
#' dist2 <- define_surv_param("gompertz", rate = .5, shape = 1)
#' dist3 <- define_surv_param("exp", rate = 0.25)
#' expect_equal(
#'  join(dist1, 5, dist2),
#'  create_list_object(
#'      c('surv_join', 'surv_combined', 'surv_dist'),
#'      dists = list(dist1, dist2),
#'      cuts = 5
#'  )
#' )
#' expect_equal(
#'  join(dist1, 5, dist2, 10, dist3),
#'  create_list_object(
#'      c('surv_join', 'surv_combined', 'surv_dist'),
#'      dists = list(dist1, dist2, dist3),
#'      cuts = c(5, 10)
#'  )
#' )
#' expect_error(
#'  join(dist1, 5, 'foo'),
#'  'Error joining distributions, invalid survival distribution provided.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  join(dist1, 'foo', dist2),
#'  'Error joining distributions, cuts times must be numeric.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  join(dist1, NA_real_, dist2),
#'  'Error joining distributions, cuts times cannot be NA.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  join(dist1, -1, dist2),
#'  'Error joining distributions, cut times cannot be negative.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  join(dist1, -1, dist2, 3),
#'  'Error joining distributions, must provide an odd number of arguments corresponding to n distributions and n - 1 cut points.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  join(dist1, 5, dist2, 3, dist3),
#'  'Error joining distributions, distributions and cutpoints must be provided in order.',
#'  fixed = TRUE
#' )
join <- function(dist1, cut1, dist2, ...) {

    # Extract, check, and prepare extra arguments
    dots <- list(...)
    n_dots <- length(dots)
    if (n_dots %% 2 != 0) {
        err <- get_and_populate_message('join_wrong_n_args')
        stop(err, call. = show_call_error())
    }
    extra_dists <- NULL
    extra_cuts <- NULL
    if (n_dots > 0) {
        extra_cuts <- dots[seq(from = 1, to = n_dots, by = 2)]
        extra_dists <- dots[seq(from = 2, to = n_dots, by = 2)]
    }

    # Compile and check distributions
    dists <- append(list(dist1, dist2), extra_dists)
    walk(dists, function(x) {
        valid <- is_surv_dist(x)
        if (!valid) {
            err <- get_and_populate_message('join_wrong_type_dist')
            stop(err, call. = show_call_error())
        }
    })

    # Compile and check cutpoints
    cuts <- imap_dbl(
        append(list(cut1), extra_cuts),
        function(x, i) {

            # Check that cut is numeric
            is_numeric <- any(c('integer', 'numeric') %in% class(x))
            if (!is_numeric) {
                err <- get_and_populate_message('join_wrong_type_cut')
                stop(err, call. = show_call_error())
            }

            # Check that cut isn't missing
            missing_cut <- any(is.na(x))
            if (missing_cut) {
                err <- get_and_populate_message('join_missing_cut')
                stop(err, call. = show_call_error())
            }

            # Check that cut >= 0
            invalid_cut <- any(x < 0)
            if (invalid_cut) {
                err <- get_and_populate_message('join_invalid_cut')
                stop(err, call. = show_call_error())
            }

            as.numeric(truncate_param(paste0('cut', i), x))
        }
    )

    # Check that cuts are in increasing order
    diffs <- diff(cuts)
    if(any(diffs < 0)) {
        err <- get_and_populate_message('join_cuts_order')
        stop(err, call. = show_call_error())
    }


    create_list_object(
        c('surv_join', 'surv_combined', 'surv_dist'),
        dists = dists,
        cuts = cuts
    )
}

#' @export 
#' 
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.05)
#' dist2 <- define_surv_param('exp', rate = 0.1)
#' dist3 <- define_surv_param('exp', rate = 0.2)
#' dist4 <- join(dist1, 20, dist2)
#' dist5 <- join(dist1, 20, dist2, 40, dist3)
#' dist6 <- join(dist1, 0.3, dist2, 2.1, dist3)
#' expect_equal(
#'  surv_prob(dist4, seq_len(100)),
#'  ppexp(
#'      seq_len(100),
#'      rate = c(0.05, 0.1),
#'      t = c(0, 20),
#'      lower.tail = FALSE
#'  )
#' )
#' expect_equal(
#'  surv_prob(dist5, seq_len(100)),
#'  ppexp(
#'      seq_len(100),
#'      rate = c(0.05, 0.1, 0.2),
#'      t = c(0, 20, 40),
#'      lower.tail = FALSE
#'  )
#' )
#' expect_equal(
#'  surv_prob(dist6, seq_len(100)),
#'  ppexp(
#'      seq_len(100),
#'      rate = c(0.05, 0.1, 0.2),
#'      t = c(0, 0.3, 2.1),
#'      lower.tail = FALSE
#'  )
#' )
surv_prob.surv_join <- function(x, time, ...) {
    check_times(time, 'calculating survival probabilities', 'time')
    n_times <- length(time)
    ret <- surv_prob(x$dists[[1]], time, ...)
    n_cuts <- length(x$cuts)
    surv_at_cut <- surv_prob(x$dists[[1]], x$cut[1], ...)
    for (i in seq_len(n_cuts)) {
        cut <- x$cut[i]
        dist <- x$dist[[i + 1]]
        surv_at_cut_this_dist <- surv_prob(dist, cut, ...)
        indicator <- time > cut
        scalar <- surv_at_cut / surv_at_cut_this_dist
        ret[indicator] <- scalar * surv_prob(dist, time[indicator])
        if (i < n_cuts) {
            next_cut <- x$cut[i + 1]
            surv_at_cut <- scalar * surv_prob(dist, next_cut, ...)
        }
    }

    ret
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.05)
#' dist2 <- define_surv_param('exp', rate = 0.1)
#'
#' # Roundtrip for 2-segment join
#' jdist <- join(dist1, 20, dist2)
#' probs <- c(0.9, 0.5, 0.1)
#' times <- surv_quantile(jdist, probs)
#' expect_equal(surv_prob(jdist, times), probs, tolerance = 1e-6)
#'
#' # Roundtrip for 3-segment join
#' dist3 <- define_surv_param('exp', rate = 0.2)
#' jdist3 <- join(dist1, 10, dist2, 30, dist3)
#' times3 <- surv_quantile(jdist3, probs)
#' expect_equal(surv_prob(jdist3, times3), probs, tolerance = 1e-6)
#'
#' # Edge cases
#' expect_equal(surv_quantile(jdist, 1), 0)
#' expect_equal(surv_quantile(jdist, 0), Inf)
#'
surv_quantile.surv_join <- function(x, probs, ...) {
    check_probs(probs)

    # Precompute survival at each cut point using the joined distribution
    n_cuts <- length(x$cuts)
    surv_at_cuts <- surv_prob(x, x$cuts, ...)

    # Precompute scale factors for each segment beyond first
    scalars <- numeric(n_cuts)
    surv_at_cut_prev <- surv_prob(x$dists[[1]], x$cuts[1], ...)
    for (i in seq_len(n_cuts)) {
        surv_at_cut_this <- surv_prob(x$dists[[i + 1]], x$cuts[i], ...)
        scalars[i] <- surv_at_cuts[i] / surv_at_cut_this
        if (i < n_cuts) {
            surv_at_cut_prev <- scalars[i] * surv_prob(x$dists[[i + 1]], x$cuts[i + 1], ...)
        }
    }

    vapply(probs, function(p) {
        if (p == 1) return(0)
        if (p == 0) return(Inf)

        # Find which segment p falls in
        # Segment 1: p >= surv_at_cuts[1]
        # Segment k: surv_at_cuts[k-1] > p >= surv_at_cuts[k]
        # Last segment: p < surv_at_cuts[n_cuts]
        seg <- 1
        for (i in seq_len(n_cuts)) {
            if (p < surv_at_cuts[i]) {
                seg <- i + 1
            } else {
                break
            }
        }

        if (seg == 1) {
            return(surv_quantile(x$dists[[1]], p, ...))
        }

        # Rescale p for the segment's distribution
        # S_joined(t) = scalar * S_seg(t) in this segment
        # So S_seg(t) = p / scalar
        p_seg <- p / scalars[seg - 1]
        surv_quantile(x$dists[[seg]], p_seg, ...)
    }, numeric(1))
}

#' @export
#' @tests
#'
#' dist1 <- define_surv_param('exp', rate = 0.05)
#' dist2 <- define_surv_param('exp', rate = 0.1)
#' expect_output(
#'  print(join(dist1, 5.12345, dist2)),
#'  'A joined survival distribution:
#'   * Segment 1 (t = 0 - 5.12): An exponential distribution (rate = 0.05).
#'   * Segment 2 (t = 5.12 - \U221E): An exponential distribution (rate = 0.1).',
#'  fixed = TRUE
#' )
print.surv_join <- function(x, ...) {
    args <- list(...)
    if (is.null(args$digits)) {
        digits <- 3
    } else {
        digits <- args$digits
    }
    formatted_cuts <- format(x$cuts, digits = digits)
    items_lines <- map_chr(seq_along(x$dists), function(i) {
        n_dists <- length(x$dists)
        dist_output <- to_list_item_output(x$dists[[i]])

        start <- '0'
        if (i != 1) {
            start <- as.character(formatted_cuts[i - 1])
        }

        end <- '\U221E'
        if (i != n_dists) {
            end <- as.character(formatted_cuts[i])
        }

        glue('    * Segment {i} (t = {start} - {end}): {dist_output}')
    })
    output <- paste0(c('A joined survival distribution:', items_lines), collapse = '\n')

    cat(output)
}