#' Define Parametric Distribution (flexsurv Parameterization)
#'
#' Define a parametric survival distribution using flexsurv's native
#' parameterization. A complete listing of supported distributions
#' is provided in the details section.
#'
#' @name define_surv_flexsurv
#' @export
#' @rdname define_surv_flexsurv
#'
#' @param distribution a parametric survival distribution.
#' @param ... additional distribution parameters
#' (see details section below)
#'
#' @return a `surv_parametric` object.
#' @details
#' Supported distributions are listed in the table below.
#'
#' | **Distribution** | **Description** | **Parameters** | **Notes** |
#' | --- | --- |  --- | --- |
#' | "exp" | Exponential | rate | |
#' | "lnorm" | Lognormal | meanlog, sdlog | |
#' | "llogis" | Log-Logistic | shape, scale | |
#' | "weibull" | Weibull (AFT) | shape, scale | |
#' | "weibullPH" | Weibull (PH) | shape, scale | |
#' | "gompertz" | Gompertz | shape, rate | |
#' | "gamma" | Gamma | shape, scale | |
#' | "gengamma" | Generalized Gamma (stable) | mu, sigma, Q | Described in Prentice (1974) |
#' | "gengamma.orig" | Generalized Gamma (original) | shape, scale, k | Described in Stacy (1962) |
#' | "genf" | Generalized F (stable) | mu, sigma, Q, P | Described in Prentice (1975) |
#' | "genf.org" | Generalized F (original) | mu, sigma, s1, s2 | Described in Prentice (1975) |
#'
#' ## Survival function formulations
#'
#' - **Exponential**: `S(t) = exp(-rate * t)`
#' - **Weibull (AFT)**: `S(t) = exp(-(t / scale)^shape)`
#' - **Weibull (PH)**: `S(t) = exp(-scale * t^shape)`
#' - **Lognormal**: `S(t) = 1 - Phi((log(t) - meanlog) / sdlog)` where
#'   `Phi()` is the standard normal CDF
#' - **Log-logistic**: `S(t) = 1 / (1 + (t / scale)^shape)`
#' - **Gompertz**: hazard `h(t) = rate * exp(shape * t)`,
#'   survival `S(t) = exp(-(rate / shape) * (exp(shape * t) - 1))`
#' - **Gamma**: No closed-form survival function; computed numerically
#'   via `pgamma()`
#' - **Generalized gamma (stable/original)**: No simple closed form;
#'   see Prentice (1974) and Stacy (1962)
#' - **Generalized F (stable/original)**: No simple closed form;
#'   see Prentice (1975)
#'
#' @references Stacy, E. W. (1962). A generalization of the gamma
#' distribution.  Annals of Mathematical Statistics 33:1187-92.
#'
#' Prentice, R. L. (1974). A log gamma model and its maximum likelihood
#' estimation. Biometrika 61(3):539-544.
#'
#' R. L. Prentice (1975). Discrimination among some parametric
#' models. Biometrika 62(3):607-614.
#'
#' @tests
#'
#' dist1 <- define_surv_flexsurv(distribution = "exp", rate = 0.05)
#' expect_equal(class(dist1), c('surv_parametric', 'surv_dist'))
#' expect_equal(dist1$distribution, 'exp')
#' expect_equal(dist1$parameters, list(rate = 0.05))
#'
#' expect_error(
#'  define_surv_flexsurv(distribution = "weibull", shape = 1.2),
#'  'Error defining Weibull (AFT) distribution, parameters missing from function call: "scale".',
#'  fixed = TRUE
#' )
#'
#' dist2 <- define_surv_flexsurv("weibull", shape = 1.5, scale = 20)
#' ref2 <- define_surv_param("weibull", shape = 1.5, scale = 20)
#' expect_equal(
#'  surv_prob(dist2, c(0, 5, 10, 20)),
#'  surv_prob(ref2, c(0, 5, 10, 20))
#' )
#'
#' @examples
#'
#' define_surv_flexsurv(distribution = "exp", rate = .5)
#' define_surv_flexsurv(distribution = "gompertz", rate = .5, shape = 1)
#'
define_surv_flexsurv <- function(distribution, ...) {

  args <- list(...)

  # Match distribution against list
  dist_string <- match.arg(distribution, choices = flexsurv_dists)

  # Extract params from arguments
  params <- get_dist_params_from_args(dist_string, args)

  # Return object
  create_list_object(
    c("surv_parametric", "surv_dist"),
    distribution = dist_string,
    parameters = params
  )
}
