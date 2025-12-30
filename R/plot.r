#' @export
#'
#' @tests
#' # Test with parametric exponential distribution
#' dist_exp <- define_surv_param(distribution = "exp", rate = 0.05)
#' p1 <- plot(dist_exp, max_time = 20)
#' expect_s3_class(p1, "ggplot")
#'
#' # Test with parametric Weibull distribution
#' dist_weibull <- define_surv_param(distribution = "weibull", shape = 1.5, scale = 10)
#' p2 <- plot(dist_weibull, max_time = 30)
#' expect_s3_class(p2, "ggplot")
#'
#' # Test with custom steps parameter
#' p3 <- plot(dist_exp, max_time = 20, steps = 100)
#' expect_s3_class(p3, "ggplot")
#'
#' # Test with cure distribution (tests ... parameter passing)
#' dist_cure <- define_surv_cure(
#'     distribution = "exp",
#'     theta = 0.3,
#'     rate = 0.1,
#'     mixture = TRUE
#' )
#' p4 <- plot(dist_cure, max_time = 50)
#' expect_s3_class(p4, "ggplot")
#'
#' # Test with custom function distribution
#' dist_func <- define_surv_func(pweibull, lower.tail = FALSE, shape = 1.2, scale = 15)
#' p5 <- plot(dist_func, max_time = 30)
#' expect_s3_class(p5, "ggplot")
#'
#' # Verify plot has correct labels (public API, stable)
#' expect_equal(p1$labels$x, "Time")
#' expect_equal(p1$labels$y, "Survival")
plot.surv_dist <- function(x, max_time, steps = 1000, ...) {
    times <- seq(from = 0, to = max_time * 1.2, length.out = steps)
    df <- tibble(
        time = times,
        survival = surv_prob(x, times, ...)
    )

    df[!is.na(df$survival), ]
    ggplot(aes(x = time, y = survival), data = df) +
        geom_line() +
        coord_cartesian(xlim = c(0, max_time), ylim = c(0, 1)) +
        labs(x = 'Time', y = 'Survival')
}