#' Generate Random Survival Times
#'
#' Sample random event times from any survival distribution using
#' inverse-transform sampling via `surv_quantile`.
#'
#' @name surv_random
#' @rdname surv_random
#' @export
#'
#' @param x A `surv_dist` object
#' @param n Number of random samples to generate
#' @param ... additional arguments passed to `surv_quantile`
#'
#' @return A numeric vector of random event times
#'
#' @examples
#' dist1 <- define_surv_param('exp', rate = 0.05)
#' samples <- surv_random(dist1, 100)
#'
#' @tests
#'
#' # Exponential: large sample mean ~ 1/rate
#' set.seed(42)
#' dist1 <- define_surv_param('exp', rate = 0.1)
#' samples <- surv_random(dist1, 5000)
#' expect_equal(mean(samples), 10, tolerance = 1)
#'
#' # Cure model: some samples are Inf
#' set.seed(42)
#' cure_dist <- define_surv_cure('exp', theta = 0.3, rate = 0.1, mixture = TRUE)
#' cure_samples <- surv_random(cure_dist, 1000)
#' prop_inf <- mean(is.infinite(cure_samples))
#' expect_equal(prop_inf, 0.3, tolerance = 0.05)
#'
#' # Validation
#' expect_error(
#'  surv_random(dist1, -1),
#'  'Error generating random survival times, "n" must be a single positive integer.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  surv_random(dist1, 0),
#'  'Error generating random survival times, "n" must be a single positive integer.',
#'  fixed = TRUE
#' )
#' expect_error(
#'  surv_random(dist1, 'foo'),
#'  'Error generating random survival times, "n" must be a single positive integer.',
#'  fixed = TRUE
#' )
#'
surv_random <- function(x, n, ...) {
    check_n(n)
    u <- runif(n)
    surv_quantile(x, probs = u, ...)
}

#' @tests
#' expect_silent(check_n(10))
#' expect_error(check_n(-1), 'must be a single positive integer', fixed = TRUE)
#' expect_error(check_n(0), 'must be a single positive integer', fixed = TRUE)
#' expect_error(check_n('foo'), 'must be a single positive integer', fixed = TRUE)
#' expect_error(check_n(c(1, 2)), 'must be a single positive integer', fixed = TRUE)
check_n <- function(n) {
    if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
        err <- get_and_populate_message('surv_random_invalid_n')
        stop(err, call. = show_call_error())
    }
}
