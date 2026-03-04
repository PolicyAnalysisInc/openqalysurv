#' Define Parametric Distribution from SAS PROC LIFEREG Parameters
#'
#' Define a parametric survival distribution using parameter estimates
#' from SAS PROC LIFEREG output. Parameters are converted from the
#' AFT parameterization to flexsurv's native parameterization.
#'
#' @name define_surv_lifereg
#' @export
#' @rdname define_surv_lifereg
#'
#' @param distribution one of `"exponential"`, `"weibull"`,
#'   `"lnormal"`, `"llogistic"`, or `"gamma"`.
#' @param intercept the Intercept from the LIFEREG Parameter
#'   Estimates table.
#' @param scale the Scale parameter from the LIFEREG output.
#'   Defaults to 1 (used for exponential).
#' @param shape the Shape parameter from the LIFEREG output.
#'   Only used for `distribution = "gamma"`.
#'
#' @return a `surv_parametric` object.
#'
#' @details
#' The following conversions are applied:
#'
#' | **SAS DIST=** | **flexsurv distribution** | **Conversion** |
#' | --- | --- | --- |
#' | EXPONENTIAL | exp | `rate = exp(-intercept)` |
#' | WEIBULL | weibull | `shape = 1/scale`, `scale = exp(intercept)` |
#' | LNORMAL | lnorm | `meanlog = intercept`, `sdlog = scale` |
#' | LLOGISTIC | llogis | `shape = 1/scale`, `scale = exp(intercept)` |
#' | GAMMA | gengamma | `mu = intercept`, `sigma = scale`, `Q = shape` |
#'
#' ## Survival function formulations
#'
#' In terms of LIFEREG parameters:
#'
#' - **Exponential**: `S(t) = exp(-exp(-intercept) * t)`
#' - **Weibull**: `S(t) = exp(-(t / exp(intercept))^(1/scale))`
#' - **Lognormal**: `S(t) = 1 - Phi((log(t) - intercept) / scale)` where
#'   `Phi()` is the standard normal CDF
#' - **Log-logistic**: `S(t) = 1 / (1 + (t / exp(intercept))^(1/scale))`
#' - **Generalized gamma**: No simple closed form;
#'   `mu = intercept`, `sigma = scale`, `Q = shape`;
#'   see Prentice (1974)
#'
#' @references Prentice, R. L. (1974). A log gamma model and its maximum
#' likelihood estimation. Biometrika 61(3):539-544.
#'
#' @tests
#'
#' # Exponential
#' dist1 <- define_surv_lifereg("exponential", intercept = 3)
#' ref1 <- define_surv_param("exp", rate = exp(-3))
#' expect_equal(class(dist1), c('surv_parametric', 'surv_dist'))
#' expect_equal(
#'  surv_prob(dist1, c(0, 5, 10, 20)),
#'  surv_prob(ref1, c(0, 5, 10, 20))
#' )
#'
#' # Weibull
#' dist2 <- define_surv_lifereg("weibull", intercept = 3, scale = 0.8)
#' ref2 <- define_surv_param("weibull", shape = 1/0.8, scale = exp(3))
#' expect_equal(
#'  surv_prob(dist2, c(0, 5, 10, 20)),
#'  surv_prob(ref2, c(0, 5, 10, 20))
#' )
#'
#' # Lognormal
#' dist3 <- define_surv_lifereg("lnormal", intercept = 2.5, scale = 0.7)
#' ref3 <- define_surv_param("lnorm", meanlog = 2.5, sdlog = 0.7)
#' expect_equal(
#'  surv_prob(dist3, c(0, 5, 10, 20)),
#'  surv_prob(ref3, c(0, 5, 10, 20))
#' )
#'
#' # Log-logistic
#' dist4 <- define_surv_lifereg("llogistic", intercept = 3, scale = 0.6)
#' ref4 <- define_surv_param("llogis", shape = 1/0.6, scale = exp(3))
#' expect_equal(
#'  surv_prob(dist4, c(0, 5, 10, 20)),
#'  surv_prob(ref4, c(0, 5, 10, 20))
#' )
#'
#' # Generalized gamma
#' dist5 <- define_surv_lifereg("gamma", intercept = 2.3, scale = 0.4, shape = -0.03)
#' ref5 <- define_surv_param("gengamma", mu = 2.3, sigma = 0.4, Q = -0.03)
#' expect_equal(
#'  surv_prob(dist5, c(0, 5, 10, 20)),
#'  surv_prob(ref5, c(0, 5, 10, 20))
#' )
#'
#' # Case insensitive
#' dist6 <- define_surv_lifereg("WEIBULL", intercept = 3, scale = 0.8)
#' expect_equal(
#'  surv_prob(dist6, c(0, 5, 10)),
#'  surv_prob(dist2, c(0, 5, 10))
#' )
#'
#' # Invalid distribution
#' expect_error(
#'  define_surv_lifereg("gompertz", intercept = 3),
#'  'is not a supported distribution',
#'  fixed = TRUE
#' )
#'
#' # Missing shape for gamma
#' expect_error(
#'  define_surv_lifereg("gamma", intercept = 2.3, scale = 0.4),
#'  'required parameters missing',
#'  fixed = TRUE
#' )
#'
#' # Invalid parameter value
#' expect_error(
#'  define_surv_lifereg("weibull", intercept = "a", scale = 0.8),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
#'
#' @examples
#'
#' define_surv_lifereg("exponential", intercept = 3)
#' define_surv_lifereg("weibull", intercept = 3, scale = 0.8)
#' define_surv_lifereg("gamma", intercept = 2.3, scale = 0.4, shape = -0.03)
#'
define_surv_lifereg <- function(distribution, intercept, scale = 1, shape = NULL) {
    source <- "PROC LIFEREG"

    # Match distribution name
    dist <- match_dist_name(distribution, lifereg_dists, source)

    # Build args list and validate
    args <- list(intercept = intercept, scale = scale)
    if (!is.null(shape)) args$shape <- shape
    required <- lifereg_dist_params[[dist]]
    check_source_params(args, required, source, dist)
    for (p in required) {
        check_numeric_param(args[[p]], p, source, dist)
    }

    # Convert to flexsurv parameters
    flexsurv_dist <- lifereg_dists[[dist]]
    flexsurv_params <- switch(dist,
        exponential = list(rate = exp(-intercept)),
        weibull = list(shape = 1 / scale, scale = exp(intercept)),
        lnormal = list(meanlog = intercept, sdlog = scale),
        llogistic = list(shape = 1 / scale, scale = exp(intercept)),
        gamma = list(mu = intercept, sigma = scale, Q = shape)
    )

    do.call(define_surv_param, c(list(distribution = flexsurv_dist), flexsurv_params))
}
