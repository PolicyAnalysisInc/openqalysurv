#' Define Parametric Distribution from survreg Parameters
#'
#' Define a parametric survival distribution using parameter estimates
#' from R's [survival::survreg()] output. Parameters are converted
#' from the AFT parameterization to flexsurv's native parameterization.
#'
#' @name define_surv_survreg
#' @export
#' @rdname define_surv_survreg
#'
#' @param distribution one of `"exponential"`, `"weibull"`,
#'   `"lognormal"`, or `"loglogistic"`.
#' @param intercept the intercept from the survreg model
#'   (`coef(model)["(Intercept)"]`).
#' @param scale the scale parameter from the survreg model
#'   (`model$scale`). Defaults to 1 (used for exponential).
#'
#' @return a `surv_parametric` object.
#'
#' @details
#' The following conversions are applied:
#'
#' | **survreg distribution** | **flexsurv distribution** | **Conversion** |
#' | --- | --- | --- |
#' | exponential | exp | `rate = exp(-intercept)` |
#' | weibull | weibull | `shape = 1/scale`, `scale = exp(intercept)` |
#' | lognormal | lnorm | `meanlog = intercept`, `sdlog = scale` |
#' | loglogistic | llogis | `shape = 1/scale`, `scale = exp(intercept)` |
#'
#' ## Survival function formulations
#'
#' In terms of survreg parameters:
#'
#' - **Exponential**: `S(t) = exp(-exp(-intercept) * t)`
#' - **Weibull**: `S(t) = exp(-(t / exp(intercept))^(1/scale))`
#' - **Lognormal**: `S(t) = 1 - Phi((log(t) - intercept) / scale)` where
#'   `Phi()` is the standard normal CDF
#' - **Log-logistic**: `S(t) = 1 / (1 + (t / exp(intercept))^(1/scale))`
#'
#' @tests
#'
#' # Exponential
#' dist1 <- define_surv_survreg("exponential", intercept = 3)
#' ref1 <- define_surv_param("exp", rate = exp(-3))
#' expect_equal(class(dist1), c('surv_parametric', 'surv_dist'))
#' expect_equal(
#'  surv_prob(dist1, c(0, 5, 10, 20)),
#'  surv_prob(ref1, c(0, 5, 10, 20))
#' )
#'
#' # Weibull
#' dist2 <- define_surv_survreg("weibull", intercept = 3, scale = 0.8)
#' ref2 <- define_surv_param("weibull", shape = 1/0.8, scale = exp(3))
#' expect_equal(
#'  surv_prob(dist2, c(0, 5, 10, 20)),
#'  surv_prob(ref2, c(0, 5, 10, 20))
#' )
#'
#' # Lognormal
#' dist3 <- define_surv_survreg("lognormal", intercept = 2.5, scale = 0.7)
#' ref3 <- define_surv_param("lnorm", meanlog = 2.5, sdlog = 0.7)
#' expect_equal(
#'  surv_prob(dist3, c(0, 5, 10, 20)),
#'  surv_prob(ref3, c(0, 5, 10, 20))
#' )
#'
#' # Log-logistic
#' dist4 <- define_surv_survreg("loglogistic", intercept = 3, scale = 0.6)
#' ref4 <- define_surv_param("llogis", shape = 1/0.6, scale = exp(3))
#' expect_equal(
#'  surv_prob(dist4, c(0, 5, 10, 20)),
#'  surv_prob(ref4, c(0, 5, 10, 20))
#' )
#'
#' # Case insensitive
#' dist5 <- define_surv_survreg("Weibull", intercept = 3, scale = 0.8)
#' expect_equal(
#'  surv_prob(dist5, c(0, 5, 10)),
#'  surv_prob(dist2, c(0, 5, 10))
#' )
#'
#' # Invalid distribution
#' expect_error(
#'  define_surv_survreg("gamma", intercept = 3),
#'  'is not a supported distribution',
#'  fixed = TRUE
#' )
#'
#' # Invalid parameter value (NA)
#' expect_error(
#'  define_surv_survreg("weibull", intercept = NA, scale = 0.8),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
#'
#' # Invalid parameter value (non-numeric)
#' expect_error(
#'  define_surv_survreg("weibull", intercept = "a", scale = 0.8),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
#'
#' @examples
#'
#' define_surv_survreg("exponential", intercept = 3)
#' define_surv_survreg("weibull", intercept = 3, scale = 0.8)
#'
define_surv_survreg <- function(distribution, intercept, scale = 1) {
    source <- "survreg"

    # Match distribution name
    dist <- match_dist_name(distribution, survreg_dists, source)

    # Build args list and validate
    args <- list(intercept = intercept, scale = scale)
    required <- survreg_dist_params[[dist]]
    check_source_params(args, required, source, dist)
    for (p in required) {
        check_numeric_param(args[[p]], p, source, dist)
    }

    # Convert to flexsurv parameters
    flexsurv_dist <- survreg_dists[[dist]]
    flexsurv_params <- switch(dist,
        exponential = list(rate = exp(-intercept)),
        weibull = list(shape = 1 / scale, scale = exp(intercept)),
        lognormal = list(meanlog = intercept, sdlog = scale),
        loglogistic = list(shape = 1 / scale, scale = exp(intercept))
    )

    do.call(define_surv_param, c(list(distribution = flexsurv_dist), flexsurv_params))
}
