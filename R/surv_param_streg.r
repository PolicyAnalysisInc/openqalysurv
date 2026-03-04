#' Define Parametric Distribution from Stata streg Parameters
#'
#' Define a parametric survival distribution using parameter estimates
#' from Stata's `streg` output. Parameters are converted from the
#' AFT or PH parameterization to flexsurv's native parameterization.
#'
#' @name define_surv_streg
#' @export
#' @rdname define_surv_streg
#'
#' @param distribution one of `"exponential"`, `"weibull"`,
#'   `"lognormal"`, `"loglogistic"`, `"ggamma"` (AFT metric), or
#'   `"gompertz"` (PH metric only).
#' @param ... distribution parameters from streg output.
#'   The `_cons` coefficient should be passed as `cons`.
#'   See details for parameter names by distribution.
#' @param metric either `"aft"` or `"ph"`, indicating which
#'   metric was used to fit the model.
#'
#' @return a `surv_parametric` object.
#'
#' @details
#'
#' ## AFT metric parameters
#'
#' | **Distribution** | **Parameters** | **Conversion** |
#' | --- | --- | --- |
#' | exponential | `cons` | `rate = exp(-cons)` |
#' | weibull | `cons`, `ln_p` | `shape = exp(ln_p)`, `scale = exp(cons)` |
#' | lognormal | `cons`, `sigma` | `meanlog = cons`, `sdlog = sigma` |
#' | loglogistic | `cons`, `gamma` | `shape = gamma`, `scale = exp(cons)` |
#' | ggamma | `cons`, `sigma`, `kappa` | `mu = cons`, `sigma = sigma`, `Q = kappa` |
#'
#' ## PH metric parameters
#'
#' | **Distribution** | **Parameters** | **Conversion** |
#' | --- | --- | --- |
#' | exponential | `cons` | `rate = exp(cons)` |
#' | weibull | `cons`, `ln_p` | `shape = exp(ln_p)`, `scale = exp(cons)` |
#' | gompertz | `cons`, `gamma` | `shape = gamma`, `rate = exp(cons)` |
#'
#' ## Survival function formulations
#'
#' ### AFT metric
#'
#' - **Exponential**: `S(t) = exp(-exp(-cons) * t)`
#' - **Weibull**: `S(t) = exp(-(t / exp(cons))^exp(ln_p))`
#' - **Lognormal**: `S(t) = 1 - Phi((log(t) - cons) / sigma)` where
#'   `Phi()` is the standard normal CDF
#' - **Log-logistic**: `S(t) = 1 / (1 + (t / exp(cons))^gamma)`
#' - **Generalized gamma**: No simple closed form;
#'   `mu = cons`, `sigma = sigma`, `Q = kappa`
#'
#' ### PH metric
#'
#' - **Exponential**: `S(t) = exp(-exp(cons) * t)`
#' - **Weibull**: `S(t) = exp(-exp(cons) * t^exp(ln_p))`
#' - **Gompertz**: hazard `h(t) = exp(cons) * exp(gamma * t)`,
#'   survival `S(t) = exp(-(exp(cons) / gamma) * (exp(gamma * t) - 1))`
#'
#' @tests
#'
#' # --- AFT metric ---
#'
#' # Exponential AFT
#' dist1 <- define_surv_streg("exponential", cons = 3, metric = "aft")
#' ref1 <- define_surv_param("exp", rate = exp(-3))
#' expect_equal(class(dist1), c('surv_parametric', 'surv_dist'))
#' expect_equal(
#'  surv_prob(dist1, c(0, 5, 10, 20)),
#'  surv_prob(ref1, c(0, 5, 10, 20))
#' )
#'
#' # Weibull AFT
#' dist2 <- define_surv_streg("weibull", cons = 3, ln_p = 0.405, metric = "aft")
#' ref2 <- define_surv_param("weibull", shape = exp(0.405), scale = exp(3))
#' expect_equal(
#'  surv_prob(dist2, c(0, 5, 10, 20)),
#'  surv_prob(ref2, c(0, 5, 10, 20))
#' )
#'
#' # Lognormal AFT
#' dist3 <- define_surv_streg("lognormal", cons = 2.5, sigma = 0.7, metric = "aft")
#' ref3 <- define_surv_param("lnorm", meanlog = 2.5, sdlog = 0.7)
#' expect_equal(
#'  surv_prob(dist3, c(0, 5, 10, 20)),
#'  surv_prob(ref3, c(0, 5, 10, 20))
#' )
#'
#' # Log-logistic AFT
#' dist4 <- define_surv_streg("loglogistic", cons = 3, gamma = 1.667, metric = "aft")
#' ref4 <- define_surv_param("llogis", shape = 1.667, scale = exp(3))
#' expect_equal(
#'  surv_prob(dist4, c(0, 5, 10, 20)),
#'  surv_prob(ref4, c(0, 5, 10, 20))
#' )
#'
#' # Generalized gamma AFT
#' dist5 <- define_surv_streg("ggamma", cons = 2.3, sigma = 0.4, kappa = -0.03, metric = "aft")
#' ref5 <- define_surv_param("gengamma", mu = 2.3, sigma = 0.4, Q = -0.03)
#' expect_equal(
#'  surv_prob(dist5, c(0, 5, 10, 20)),
#'  surv_prob(ref5, c(0, 5, 10, 20))
#' )
#'
#' # --- PH metric ---
#'
#' # Exponential PH
#' dist6 <- define_surv_streg("exponential", cons = -3, metric = "ph")
#' ref6 <- define_surv_param("exp", rate = exp(-3))
#' expect_equal(
#'  surv_prob(dist6, c(0, 5, 10, 20)),
#'  surv_prob(ref6, c(0, 5, 10, 20))
#' )
#'
#' # Weibull PH
#' dist7 <- define_surv_streg("weibull", cons = -2, ln_p = 0.405, metric = "ph")
#' ref7 <- define_surv_param("weibullPH", shape = exp(0.405), scale = exp(-2))
#' expect_equal(
#'  surv_prob(dist7, c(0, 5, 10, 20)),
#'  surv_prob(ref7, c(0, 5, 10, 20))
#' )
#'
#' # Gompertz PH
#' dist8 <- define_surv_streg("gompertz", cons = -3, gamma = 0.05, metric = "ph")
#' ref8 <- define_surv_param("gompertz", shape = 0.05, rate = exp(-3))
#' expect_equal(
#'  surv_prob(dist8, c(0, 5, 10, 20)),
#'  surv_prob(ref8, c(0, 5, 10, 20))
#' )
#'
#' # Case insensitive
#' dist9 <- define_surv_streg("Weibull", cons = 3, ln_p = 0.405, metric = "aft")
#' expect_equal(
#'  surv_prob(dist9, c(0, 5, 10)),
#'  surv_prob(dist2, c(0, 5, 10))
#' )
#'
#' # Invalid distribution for metric
#' expect_error(
#'  define_surv_streg("gompertz", cons = 3, gamma = 0.05, metric = "aft"),
#'  'is not a supported distribution',
#'  fixed = TRUE
#' )
#'
#' # Invalid distribution
#' expect_error(
#'  define_surv_streg("gamma", cons = 3, metric = "aft"),
#'  'is not a supported distribution',
#'  fixed = TRUE
#' )
#'
#' # Missing parameter
#' expect_error(
#'  define_surv_streg("weibull", cons = 3, metric = "aft"),
#'  'required parameters missing',
#'  fixed = TRUE
#' )
#'
#' # Invalid parameter value
#' expect_error(
#'  define_surv_streg("weibull", cons = NA, ln_p = 0.5, metric = "aft"),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
#'
#' @examples
#'
#' define_surv_streg("exponential", cons = 3, metric = "aft")
#' define_surv_streg("weibull", cons = 3, ln_p = 0.405, metric = "aft")
#' define_surv_streg("gompertz", cons = -3, gamma = 0.05, metric = "ph")
#'
define_surv_streg <- function(distribution, ..., metric = c("aft", "ph")) {
    source <- "streg"
    metric <- match.arg(metric)

    # Select mapping and param lists based on metric
    if (metric == "aft") {
        valid_dists <- streg_aft_dists
        param_map <- streg_aft_dist_params
    } else {
        valid_dists <- streg_ph_dists
        param_map <- streg_ph_dist_params
    }

    # Match distribution name
    dist <- match_dist_name(distribution, valid_dists, source)

    # Collect and validate parameters
    args <- list(...)
    required <- param_map[[dist]]
    check_source_params(args, required, source, dist)
    for (p in required) {
        check_numeric_param(args[[p]], p, source, dist)
    }

    # Convert to flexsurv parameters
    flexsurv_dist <- valid_dists[[dist]]

    if (metric == "aft") {
        flexsurv_params <- switch(dist,
            exponential = list(rate = exp(-args$cons)),
            weibull = list(shape = exp(args$ln_p), scale = exp(args$cons)),
            lognormal = list(meanlog = args$cons, sdlog = args$sigma),
            loglogistic = list(shape = args$gamma, scale = exp(args$cons)),
            ggamma = list(mu = args$cons, sigma = args$sigma, Q = args$kappa)
        )
    } else {
        flexsurv_params <- switch(dist,
            exponential = list(rate = exp(args$cons)),
            weibull = list(shape = exp(args$ln_p), scale = exp(args$cons)),
            gompertz = list(shape = args$gamma, rate = exp(args$cons))
        )
    }

    do.call(define_surv_param, c(list(distribution = flexsurv_dist), flexsurv_params))
}
