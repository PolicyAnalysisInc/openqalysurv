#' Apply Fixed Shift in Time
#' 
#' Apply a fixed shift in time to move the hazards of a survival
#' distribution forwards or backwards.
#' 
#' @name apply_shift
#' @rdname apply_shift
#' @export
#' 
#' @param dist a survival distribution
#' @param shift amount of time to shift the distribution. A positive
#' value delays hazards (shifts the survival curve to the right), and a
#' negative value accelerates hazards (shifts the survival curve to the left).
#'   
#' @return A `surv_shift` object.
#'
#' @details
#' The shifted survival function is defined as:
#'
#' \deqn{S_{\mathrm{shifted}}(t) = \begin{cases} 1 & t \le \delta \\
#' S(t - \delta) & t > \delta \end{cases}}{S_shifted(t) = 1 for
#' t <= delta; S_shifted(t) = S(t - delta) for t > delta}
#'
#' where \eqn{\delta}{delta} is the shift amount and \eqn{S} is the
#' baseline survival function.
#'
#' A positive shift delays all hazards by \eqn{\delta}{delta} time units
#' (the curve moves to the right); a negative shift accelerates them
#' (the curve moves to the left).
#'
#' A common use of \code{apply_shift} is to realign the timescale of a
#' distribution before passing it to \code{\link{join}}. When a
#' distribution is estimated from data that begins at a cut point
#' \eqn{c} — so that its own time zero corresponds to \eqn{c} — applying
#' \code{apply_shift(dist, c)} maps its internal time to the absolute
#' study timescale, which \code{join} requires.
#'
#' @examples
#' 
#' shifted_dist <- apply_shift(
#'  define_surv_param("exp", rate = 0.25),
#'  3
#' )
#' 
#' @tests
#' dist1 <- define_surv_param("exp", rate = 0.25)
#' expect_equal(
#'  apply_shift(dist1, 2),
#'  create_list_object(c('surv_shift', 'surv_dist'), dist = dist1, shift = 2)
#' )
#' expect_equal(
#'  apply_shift(dist1, 4),
#'  apply_shift(apply_shift(dist1, 2), 2)
#' )
#' expect_error(
#'  apply_shift('foo', 0.5),
#'  'Error applying shift, invalid survival distribution provided.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_shift(dist1, 'foo'),
#'  'Error applying shift, "shift" must be numeric.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  apply_shift(dist1, NA_real_),
#'  'Error applying shift, "shift" cannot be NA.',
#'  fixed = TRUE
#' )
apply_shift <- function(dist, shift) {

    # Check that dist is a valid type
    is_surv_dist <- is_surv_dist(dist)
    if (!is_surv_dist) {
        err <- get_and_populate_message('apply_shift_wrong_type_dist')
        stop(err, call. = show_call_error())
    }

    # Check that shift is numeric
    is_numeric <- any(c('integer', 'numeric') %in% class(shift))
    if (!is_numeric) {
        err <- get_and_populate_message('apply_shift_wrong_type_shift')
        stop(err, call. = show_call_error())
    }

    # Check that shift isn't missing
    missing_shift <- any(is.na(shift))
    if (missing_shift) {
        err <- get_and_populate_message('apply_shift_missing_shift')
        stop(err, call. = show_call_error())
    }

    shift <- truncate_param('shift', shift)

    # If shift equals zero then noop
    if (shift == 0) {
        return(dist)
    }

    # If the baseline distribution is of type surv_shift
    # then we can just add the shifts.
    if (inherits(dist, 'surv_shift')) {
        dist$shift <- dist$shift + shift
        return(dist)
    }

    create_list_object(
        c('surv_shift', 'surv_dist'),
        dist = dist,
        shift = shift
    )
  
}

#' @export
#' 
#' @tests
#' 
#' dist1 <- define_surv_param("exp", rate = 0.50)
#' dist2 <- apply_shift(dist1, 2)
#' expect_equal(
#'  surv_prob(dist1, seq_len(100)),
#'  surv_prob(dist2, seq_len(100) + 2)
#' )
surv_prob.surv_shift <- function(x, time, ...) {

    # Check that times are valid
    check_times(time, 'calculating survival probabilities', 'time')

    # Create a vector to store results
    ret <- vector('numeric', length(time))

    # Shift the times and figure out which times are lte zero
    shifted_times <- time - x$shift
    zero_times <- shifted_times <= 0

    # Survival for t <= 0 is 1 by definition
    ret[zero_times] <- 1

    # Survival for t > 0 given by surv_prob
    ret[!zero_times] <- surv_prob(x$dist, shifted_times[!zero_times])

    ret
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.1)
#'
#' # Shift algebra: Q_shift(p) = shift + Q_inner(p)
#' expect_equal(
#'  surv_quantile(apply_shift(dist1, 5), 0.5),
#'  5 + surv_quantile(dist1, 0.5),
#'  tolerance = 1e-10
#' )
#'
#' # Roundtrip
#' shifted <- apply_shift(dist1, 3)
#' probs <- c(0.9, 0.5, 0.1)
#' times <- surv_quantile(shifted, probs)
#' expect_equal(surv_prob(shifted, times), probs, tolerance = 1e-6)
#'
#' # Negative shift: clamp to 0
#' shifted_neg <- apply_shift(dist1, -100)
#' expect_equal(surv_quantile(shifted_neg, 0.999), 0)
#'
surv_quantile.surv_shift <- function(x, probs, ...) {
    check_probs(probs)
    vapply(probs, function(p) {
        if (p == 1) return(0)
        if (p == 0) return(Inf)
        max(0, x$shift + surv_quantile(x$dist, p, ...))
    }, numeric(1))
}

#' @export
#'
#' @tests
#' dist1 <- apply_shift(define_surv_param('exp', rate = 0.025), 2.5)
#' expect_output(
#'  print(dist1),
#'  'A shifted survival distribution:
#'   * Shift: 2.5
#'   * Baseline Distribution: An exponential distribution (rate = 0.025).',
#'  fixed = T
#' )
#' 
print.surv_shift <- function(x, ...) {
    bl_dist_output <- to_list_item_output(x$dist)
    output <- paste0(
        c(
            'A shifted survival distribution:',
            glue('    * Shift: {x$shift}'),
            glue('    * Baseline Distribution: {bl_dist_output}')
        ),
        collapse = '\n'
    )
    cat(output)
}