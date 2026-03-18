#' Define Parametric Distribution
#'
#' Alias for [define_surv_flexsurv()]. Define a parametric survival
#' distribution using flexsurv's native parameterization.
#'
#' @name define_surv_param
#' @export
#' @rdname define_surv_param
#'
#' @param distribution a parametric survival distribution
#'   (see [define_surv_flexsurv()] for supported distributions).
#' @param ... distribution parameters
#'   (see [define_surv_flexsurv()] for details).
#'
#' @return a `surv_parametric` object.
#'
#' @tests
#'
#' dist1 <- define_surv_param(distribution = "exp", rate = 0.05)
#' expect_equal(class(dist1), c('surv_parametric', 'surv_dist'))
#' expect_equal(dist1$distribution, 'exp')
#' expect_equal(dist1$parameters, list(rate = 0.05))
#'
#' expect_error(
#'  define_surv_param(distribution = "weibull", shape = 1.2),
#'  'Error defining Weibull (AFT) distribution, parameters missing from function call: "scale".',
#'  fixed = TRUE
#' )
#'
#'
#' @examples
#'
#' define_surv_param(distribution = "exp", rate = .5)
#' define_surv_param(distribution = "gompertz", rate = .5, shape = 1)
#'
define_surv_param <- function(distribution, ...) {
    define_surv_flexsurv(distribution, ...)
}

#' @export
#' @tests
#' 
#' surv_dist1 <- define_surv_param('weibull', shape = 1.2438, scale = 20.3984)
#' expect_output(
#'  print(surv_dist1),
#'  "A Weibull (AFT) distribution (shape = 1.24, scale = 20.40).",
#'  fixed = TRUE
#' )
#' 
#' surv_dist2 <- define_surv_param('exp', rate = 0.34)
#' expect_output(
#'  print(surv_dist2),
#'  "An exponential distribution (rate = 0.34).",
#'  fixed = TRUE
#' )
print.surv_parametric <- function(x, ...) {
    formatter <- create_param_formatter(...)
    dist_name <- get_dist_display_name(x$dist)
    indef_article <- str_to_title(get_indefinite_article(dist_name))
    param_string <- paste(
        paste(names(x$parameters), '=', formatter(as.numeric(x$parameters))),
        collapse = ', '
    )
    output <- glue('{indef_article} {dist_name} distribution ({param_string}).')

    cat(output)
}

#' @export
#' 
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.12)
#' expect_equal(
#'  surv_prob(dist1, c(0, 1, 2, 3)),
#'  c(1.0000000, 0.8869204, 0.7866279, 0.6976763),
#'  tolerance = 0.00001
#' )
#' 
#' dist1 <- define_surv_param('gengamma', mu = 2.321, sigma = 0.434, Q = -0.034)
#' expect_equal(
#'  surv_prob(dist1, c(0, 1, 2, 3)),
#'  c(1.0000000, 1.0000000, 0.9999393, 0.9979701),
#'  tolerance = 0.00001
#' )
surv_prob.surv_parametric <- function(x, time, ...) {
    
    check_times(time, 'calculating survival probabilities', 'time')

    # Collect extra arguments
    args <- list(...)

    # Get survival distribution function
    surv_dist <- get_flexsurv_dist(x$distribution)

    # Put together arguments for call to survival distribution
    args_for_surv_dist <- append(
        list(q = time, lower.tail = FALSE),
        x$parameters
    )

    # Call survival distribution function with arguments
    ret <- do.call(surv_dist, args_for_surv_dist)

    ret
}

#' @export
#'
#' @tests
#' # Exponential matches R builtin
#' dist_exp <- define_surv_param('exp', rate = 0.12)
#' expect_equal(
#'  surv_quantile(dist_exp, 0.5),
#'  qexp(0.5, rate = 0.12),
#'  tolerance = 1e-10
#' )
#'
#' # Weibull matches R builtin
#' dist_weib <- define_surv_param('weibull', shape = 1.5, scale = 20)
#' expect_equal(
#'  surv_quantile(dist_weib, 0.5),
#'  qweibull(0.5, shape = 1.5, scale = 20),
#'  tolerance = 1e-10
#' )
#'
#' # Roundtrip
#' probs <- c(0.9, 0.5, 0.1)
#' times <- surv_quantile(dist_weib, probs)
#' expect_equal(surv_prob(dist_weib, times), probs, tolerance = 1e-10)
#'
#' # Edge cases
#' expect_equal(surv_quantile(dist_exp, 1), 0)
#' expect_equal(surv_quantile(dist_exp, 0), Inf)
#'
surv_quantile.surv_parametric <- function(x, probs, ...) {
    check_probs(probs)
    q_func <- get_flexsurv_quantile(x$distribution)
    args_for_q <- append(list(p = probs, lower.tail = FALSE), x$parameters)
    result <- do.call(q_func, args_for_q)
    result[probs == 1] <- 0
    result[probs == 0] <- Inf
    result
}

#' @tests
#' expect_equal(get_flexsurv_dist('weibull'), pweibull)
#' expect_equal(get_flexsurv_dist('genf'), pgenf)
#' expect_equal(get_flexsurv_dist('llogis'), pllogis)
get_flexsurv_dist <- function(dist_name) {
    get(paste0("p", dist_name))
}

#' @tests
#' expect_equal(get_flexsurv_quantile('weibull'), qweibull)
#' expect_equal(get_flexsurv_quantile('exp'), qexp)
get_flexsurv_quantile <- function(dist_name) {
    get(paste0("q", dist_name))
}

#' @tests
#' expect_equal(
#'  get_flexsurv_dist_params('weibull'), c('shape', 'scale')
#' )
#' expect_equal(
#'  get_flexsurv_dist_params('gengamma'),
#'  c('mu', 'sigma', 'Q')
#' )
#' expect_equal(
#'  get_flexsurv_dist_params('genf'),
#'  c('mu', 'sigma', 'Q', 'P')
#' )
get_flexsurv_dist_params <- function(dist_name) {
    dist <- get_flexsurv_dist(dist_name)
    all_param_names <- names(formals(dist))
    dist_param_names <- setdiff(all_param_names, c('q', 'lower.tail', 'log.p'))

    dist_param_names
}

#' @tests
#' expect_equal(
#'  get_dist_params_from_args(
#'      'weibull',
#'      list(foo=1,shape=2,scale=c(3,3,3),bar=4)
#' ),
#'  list(shape=2,scale=3)
#' )
#' 
get_dist_params_from_args <- function(distribution, args) {

  # Extract parameter names
  param_names <- get_flexsurv_dist_params(distribution)

  if (!is.null(names(args)) || length(args) != length(param_names)) {
    # Run checks
    check_param_names(args, distribution)
  } else {
      names(args) <- param_names
  }

  # Create named list with parameters
  params <- map(param_names, function(name) get_dist_param_from_args(name, args))
  names(params) <- param_names

  params
}

#' @tests
#' expect_equal(
#'  get_dist_param_from_args(
#'      'scale',
#'      list(foo=1,shape=2,scale=c(3,3,3),bar=4)
#' ),
#'  3
#' )
get_dist_param_from_args <- function(name, args) {

    values <- args[[name]]
    truncate_param(name, values)

}

#' @tests
#' 
#' expect_equal(
#'  get_dist_display_name('foo'),
#'  'foo'
#' )
#' 
#' expect_equal(
#'  get_dist_display_name('exp'),
#'  'exponential'
#' )
get_dist_display_name <- function(name) {
    if (!name %in% names(flexsurv_dist_aliases)) {
        return(name)
    }

    flexsurv_dist_aliases[[name]]
}

#' @tests
#' expect_error(
#'  check_param_names(list(shape=1,foo=2), 'weibullPH'), 
#'  'Error defining Weibull (PH) distribution, parameters missing from function call: "scale".',
#'  fixed = T
#' )
#' 
check_param_names <- function(params, dist) {
    surv_func_params <- get_flexsurv_dist_params(dist)
    missing_params <- surv_func_params[!surv_func_params %in% names(params)]
    if (length(missing_params) > 0) {
        dist_str <- get_dist_display_name(dist)
        param_str <- quoted_list_string(missing_params)
        err <- get_and_populate_message(
            'missing_parameters',
            dist = dist_str,
            params = param_str
        )
        stop(err, call. = show_call_error())
    }
}