#' Evaluate Survival Quantiles
#'
#' Find the time t at which the survival probability equals a given
#' value.  That is, find t such that S(t) = p.
#'
#' @name surv_quantile
#' @rdname surv_quantile
#' @export
#'
#' @param x A `surv_dist` object
#' @param probs A numeric vector of survival probabilities in \[0, 1\]
#' @param ... additional arguments passed to methods
#'
#' @return A numeric vector of times
#'
#' @examples
#' dist1 <- define_surv_param('exp', rate = 0.12)
#' surv_quantile(dist1, c(0.9, 0.5, 0.1))
#'
surv_quantile <- function(x, probs, ...) {
    UseMethod("surv_quantile", x)
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.05)
#'
#' # Roundtrip: surv_prob(dist, surv_quantile(dist, p)) ~ p
#' probs <- c(0.9, 0.5, 0.1)
#' times <- surv_quantile(dist1, probs)
#' expect_equal(surv_prob(dist1, times), probs, tolerance = 1e-6)
#'
#' # Edge cases
#' expect_equal(surv_quantile(dist1, 1), 0)
#' expect_equal(surv_quantile(dist1, 0), Inf)
#'
#' # Validation
#' expect_error(
#'  surv_quantile(dist1, 1.5),
#'  'Error computing survival quantiles, "probs" must be numeric values in the interval [0, 1].',
#'  fixed = TRUE
#' )
#' expect_error(
#'  surv_quantile(dist1, -0.1),
#'  'Error computing survival quantiles, "probs" must be numeric values in the interval [0, 1].',
#'  fixed = TRUE
#' )
#' expect_error(
#'  surv_quantile(dist1, NA_real_),
#'  'Error computing survival quantiles, "probs" cannot contain NA values.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  surv_quantile(dist1, 'foo'),
#'  'Error computing survival quantiles, "probs" must be numeric values in the interval [0, 1].',
#'  fixed = TRUE
#' )
#'
surv_quantile.default <- function(x, probs, ...) {
    check_probs(probs)
    vapply(probs, function(p) {
        if (p == 1) return(0)
        if (p == 0) return(Inf)
        upper <- find_search_upper(x, p)
        if (is.infinite(upper)) return(Inf)
        uniroot(function(t) surv_prob(x, t) - p, c(0, upper),
                tol = .Machine$double.eps^0.5)$root
    }, numeric(1))
}

#' @tests
#' expect_silent(check_probs(c(0, 0.5, 1)))
#' expect_error(check_probs(1.5), 'must be numeric values', fixed = TRUE)
#' expect_error(check_probs(-0.1), 'must be numeric values', fixed = TRUE)
#' expect_error(check_probs(NA_real_), 'cannot contain NA', fixed = TRUE)
#' expect_error(check_probs('foo'), 'must be numeric values', fixed = TRUE)
check_probs <- function(probs) {
    if (!is.numeric(probs) || any(probs < 0, na.rm = TRUE) || any(probs > 1, na.rm = TRUE)) {
        err <- get_and_populate_message('surv_quantile_invalid_probs')
        stop(err, call. = show_call_error())
    }
    if (any(is.na(probs))) {
        err <- get_and_populate_message('surv_quantile_missing_probs')
        stop(err, call. = show_call_error())
    }
}

#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.05)
#' expect_true(is.finite(find_search_upper(dist1, 0.5)))
#' expect_true(surv_prob(dist1, find_search_upper(dist1, 0.01)) < 0.01)
#'
#' # A Kaplan-Meier quantile is the first time S <= p, and Inf when S never
#' # reaches p within follow-up
#' sf_bc <- survival::survfit(survival::Surv(recyrs, censrec) ~ 1, data = flexsurv::bc)
#' q_median <- surv_quantile(sf_bc, 0.5)
#' expect_true(is.finite(q_median))
#' expect_equal(q_median, min(sf_bc$time[sf_bc$surv <= 0.5]), tolerance = 1e-4)
#' expect_equal(surv_quantile(sf_bc, 0.01), Inf)
#'
#' # Stratified fits search only up to the last time every stratum covers
#' sf_strata <- survival::survfit(survival::Surv(recyrs, censrec) ~ group, data = flexsurv::bc)
#' expect_true(is.finite(surv_quantile(sf_strata, 0.5)))
find_search_upper <- function(x, p) {
    cap <- search_time_cap(x)
    upper <- min(1, cap)
    max_iter <- 50
    for (i in seq_len(max_iter)) {
        s <- surv_prob(x, upper)
        if (!is.na(s) && s <= p) return(upper)
        if (upper >= cap) return(Inf)
        upper <- min(upper * 2, cap)
    }
    Inf
}

# Largest time at which survival can be evaluated: for a survfit object, the
# last time with a non-missing prediction (survival is NA beyond follow-up
# unless the curve reaches 0); Inf for other distributions.
search_time_cap <- function(x) {
    if (!inherits(x, "survfit")) return(Inf)
    cand <- sort(unique(c(0, x$time)))
    ok <- !is.na(surv_prob(x, cand))
    if (!any(ok)) return(0)
    max(cand[ok])
}
