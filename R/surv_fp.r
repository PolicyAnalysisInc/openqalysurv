

#' Define Fractional Polynomial Survival Distribution
#'
#' Define a fractional polynomial (FP) survival distribution on the
#' log-hazard scale, as used in network meta-analysis of survival data.
#' Supports FP1 (one power) and FP2 (two powers) models.
#'
#' The log-hazard is modeled as:
#' \itemize{
#'   \item FP1: \code{log h(t) = beta0 + beta1 * t^(p1)}
#'   \item FP2: \code{log h(t) = beta0 + beta1 * t^(p1) + beta2 * t^(p2)}
#' }
#'
#' where \code{t^0 = log(t)} by convention, and repeated powers use
#' the Royston-Altman convention of multiplying by \code{log(t)}.
#'
#' Survival probabilities are computed as \code{S(t) = exp(-H(t))} where
#' \code{H(t)} is the cumulative hazard obtained by numerical integration
#' of the hazard function.
#'
#' @name define_surv_fp
#' @rdname define_surv_fp
#' @export
#'
#' @param betas numeric vector of coefficients. Length must equal
#'   \code{length(powers) + 1} (intercept plus one per power).
#' @param powers numeric vector of fractional polynomial powers
#'   (length 1 for FP1, length 2 for FP2).
#'
#' @return a \code{surv_fp} object.
#'
#' @details
#' The standard power set is \{-2, -1, -0.5, 0, 0.5, 1, 2, 3\}, but
#' any real-valued powers are accepted.
#'
#' @references
#' Jansen, J. P. (2011). Network meta-analysis of survival data with
#' fractional polynomials. BMC Medical Research Methodology, 11, 61.
#'
#' Royston, P. and Altman, D.G. (1994). Regression using fractional
#' polynomials of continuous covariates: parsimonious parametric modelling.
#' Journal of the Royal Statistical Society Series C, 43, 429-467.
#'
#' @examples
#'
#' # FP1 model
#' define_surv_fp(betas = c(-3, 0.5), powers = c(1))
#'
#' # FP2 model
#' define_surv_fp(betas = c(-3, 0.5, -0.2), powers = c(0, 1))
#'
#' @tests
#'
#' # Basic construction
#' fp1 <- define_surv_fp(betas = c(-3, 0.5), powers = c(1))
#' expect_equal(class(fp1), c('surv_fp', 'surv_dist'))
#' expect_equal(fp1$betas, c(beta0 = -3, beta1 = 0.5))
#' expect_equal(fp1$powers, c(1))
#'
#' # FP2 construction
#' fp2 <- define_surv_fp(betas = c(-3, 0.5, -0.2), powers = c(0, 1))
#' expect_equal(class(fp2), c('surv_fp', 'surv_dist'))
#' expect_equal(fp2$betas, c(beta0 = -3, beta1 = 0.5, beta2 = -0.2))
#'
#' # Error: non-numeric betas
#' expect_error(
#'  define_surv_fp(betas = "foo", powers = c(1)),
#'  '"betas" must be numeric',
#'  fixed = TRUE
#' )
#'
#' # Error: non-numeric powers
#' expect_error(
#'  define_surv_fp(betas = c(-3, 0.5), powers = "bar"),
#'  '"powers" must be numeric',
#'  fixed = TRUE
#' )
#'
#' # Error: NA in betas
#' expect_error(
#'  define_surv_fp(betas = c(-3, NA), powers = c(1)),
#'  '"betas" cannot contain NA',
#'  fixed = TRUE
#' )
#'
#' # Error: NA in powers
#' expect_error(
#'  define_surv_fp(betas = c(-3, 0.5), powers = c(NA_real_)),
#'  '"powers" cannot contain NA',
#'  fixed = TRUE
#' )
#'
#' # Error: length mismatch
#' expect_error(
#'  define_surv_fp(betas = c(-3, 0.5, 0.1), powers = c(1)),
#'  '"betas" must have length equal to length of "powers" plus one',
#'  fixed = TRUE
#' )
#'
#' # Error: empty powers
#' expect_error(
#'  define_surv_fp(betas = c(-3), powers = numeric(0)),
#'  'must specify at least one power',
#'  fixed = TRUE
#' )
#'
define_surv_fp <- function(betas, powers) {

    # Validate betas type
    if (!inherits(betas, c('numeric', 'integer'))) {
        err <- get_and_populate_message('fp_betas_wrong_type')
        stop(err, call. = show_call_error())
    }

    # Validate powers type
    if (!inherits(powers, c('numeric', 'integer'))) {
        err <- get_and_populate_message('fp_powers_wrong_type')
        stop(err, call. = show_call_error())
    }

    # Validate no NA in betas
    if (any(is.na(betas))) {
        err <- get_and_populate_message('fp_betas_missing')
        stop(err, call. = show_call_error())
    }

    # Validate no NA in powers
    if (any(is.na(powers))) {
        err <- get_and_populate_message('fp_powers_missing')
        stop(err, call. = show_call_error())
    }

    # Validate at least one power
    if (length(powers) < 1) {
        err <- get_and_populate_message('fp_min_powers')
        stop(err, call. = show_call_error())
    }

    # Validate length of betas = length of powers + 1
    if (length(betas) != length(powers) + 1) {
        err <- get_and_populate_message('fp_length_mismatch')
        stop(err, call. = show_call_error())
    }

    # Name betas
    beta_names <- paste0('beta', seq_along(betas) - 1)
    names(betas) <- beta_names

    create_list_object(
        c("surv_fp", "surv_dist"),
        betas = betas,
        powers = powers
    )
}


#' @tests
#'
#' # p=0 means log(t)
#' expect_equal(fp_transform_time(exp(1), c(0), 1), 1)
#'
#' # p=1 means t^1 = t
#' expect_equal(fp_transform_time(5, c(1), 1), 5)
#'
#' # p=2 means t^2
#' expect_equal(fp_transform_time(3, c(2), 1), 9)
#'
#' # p=-1 means t^(-1) = 1/t
#' expect_equal(fp_transform_time(4, c(-1), 1), 0.25)
#'
#' # Repeated powers: second index multiplied by log(t)
#' expect_equal(
#'  fp_transform_time(5, c(1, 1), 2),
#'  5 * log(5)
#' )
#'
#' # Repeated powers with p=0: t^0 * log(t) = log(t) * log(t)
#' expect_equal(
#'  fp_transform_time(exp(2), c(0, 0), 2),
#'  2 * 2
#' )
#'
fp_transform_time <- function(t, powers, index) {
    p <- powers[index]
    if (p == 0) {
        result <- log(t)
    } else {
        result <- t^p
    }

    # Repeated power convention (Royston-Altman 1994)
    if (index > 1 && powers[index] == powers[index - 1]) {
        result <- result * log(t)
    }

    result
}


#' @tests
#'
#' # FP1: log h(t) = beta0 + beta1 * t^p1
#' expect_equal(
#'  fp_log_hazard(c(-3, 0.5), c(1), 2),
#'  -3 + 0.5 * 2
#' )
#'
#' # FP2: log h(t) = beta0 + beta1 * t^p1 + beta2 * t^p2
#' expect_equal(
#'  fp_log_hazard(c(-3, 0.5, -0.2), c(0, 1), 2),
#'  -3 + 0.5 * log(2) + -0.2 * 2
#' )
#'
#' # Vectorized over t
#' expect_equal(
#'  fp_log_hazard(c(-3, 0.5), c(1), c(1, 2, 3)),
#'  c(-3 + 0.5 * 1, -3 + 0.5 * 2, -3 + 0.5 * 3)
#' )
#'
fp_log_hazard <- function(betas, powers, t) {
    result <- rep(betas[1], length(t))
    for (i in seq_along(powers)) {
        result <- result + betas[i + 1] * fp_transform_time(t, powers, i)
    }
    result
}


#' @tests
#'
#' # Cumulative hazard at t=0 should be 0
#' expect_equal(
#'  fp_cumulative_hazard(c(-3, 0.5), c(1), 0),
#'  0
#' )
#'
#' # Cumulative hazard should be positive for t > 0
#' ch <- fp_cumulative_hazard(c(-3, 0.5), c(1), 5)
#' expect_true(ch > 0)
#'
#' # Vectorized: multiple time points
#' ch_vec <- fp_cumulative_hazard(c(-3, 0.5), c(1), c(0, 1, 5))
#' expect_equal(ch_vec[1], 0)
#' expect_true(all(diff(ch_vec) >= 0))
#'
fp_cumulative_hazard <- function(betas, powers, t) {
    vapply(t, function(ti) {
        if (ti == 0) return(0)
        tryCatch(
            integrate(
                function(u) exp(fp_log_hazard(betas, powers, u)),
                lower = 0,
                upper = ti,
                rel.tol = 1e-8,
                subdivisions = 200
            )$value,
            error = function(e) {
                warning(
                    "Numerical integration failed at t=", ti, ": ",
                    e$message
                )
                NA_real_
            }
        )
    }, numeric(1))
}


#' @export
#' @tests
#'
#' fp1 <- define_surv_fp(betas = c(-3, 0.5), powers = c(1))
#' expect_output(
#'  print(fp1),
#'  "An FP1 survival distribution (beta0 = -3.0, beta1 = 0.5, powers = [1]).",
#'  fixed = TRUE
#' )
#'
#' fp2 <- define_surv_fp(betas = c(-3, 0.5, -0.2), powers = c(0, 1))
#' expect_output(
#'  print(fp2),
#'  "An FP2 survival distribution (beta0 = -3.0, beta1 = 0.5, beta2 = -0.2, powers = [0, 1]).",
#'  fixed = TRUE
#' )
#'
print.surv_fp <- function(x, ...) {
    formatter <- create_param_formatter(...)
    order <- length(x$powers)
    fp_label <- paste0('FP', order)
    beta_str <- paste(
        names(x$betas), '=', formatter(as.numeric(x$betas)),
        collapse = ', '
    )
    power_str <- paste0(formatter(x$powers), collapse = ', ')
    output <- glue(
        'An {fp_label} survival distribution ({beta_str}, powers = [{power_str}]).'
    )
    cat(output)
}


#' @export
#'
#' @tests
#'
#' # ---- Equivalence: FP1 p=1 is Gompertz ----
#' # Gompertz hazard: h(t) = rate * exp(shape * t)
#' # FP1 p=1: h(t) = exp(beta0 + beta1 * t)
#' # => rate = exp(beta0), shape = beta1
#' fp_gomp <- define_surv_fp(
#'  betas = c(log(0.1), 0.05),
#'  powers = c(1)
#' )
#' gomp <- define_surv_param('gompertz', shape = 0.05, rate = 0.1)
#' expect_equal(
#'  surv_prob(fp_gomp, c(0, 1, 5, 10, 20)),
#'  surv_prob(gomp, c(0, 1, 5, 10, 20)),
#'  tolerance = 1e-5
#' )
#'
#' # ---- Equivalence: FP1 p=0 is Weibull PH ----
#' # Weibull PH hazard: h(t) = shape * scale * t^(shape-1)
#' # FP1 p=0: h(t) = exp(beta0 + beta1 * log(t))
#' #        = exp(beta0) * t^beta1
#' # => shape = beta1 + 1, scale = exp(beta0) / shape
#' fp_weib <- define_surv_fp(
#'  betas = c(log(0.15), 0.5),
#'  powers = c(0)
#' )
#' weib_shape <- 0.5 + 1
#' weib_scale <- exp(log(0.15)) / weib_shape
#' weib <- define_surv_param('weibullPH', shape = weib_shape, scale = weib_scale)
#' expect_equal(
#'  surv_prob(fp_weib, c(0, 1, 5, 10, 20)),
#'  surv_prob(weib, c(0, 1, 5, 10, 20)),
#'  tolerance = 1e-5
#' )
#'
#' # ---- Equivalence: FP1 beta1=0 is Exponential ----
#' # h(t) = exp(beta0) = constant => exponential with rate = exp(beta0)
#' fp_exp <- define_surv_fp(
#'  betas = c(log(0.05), 0),
#'  powers = c(1)
#' )
#' exp_dist <- define_surv_param('exp', rate = 0.05)
#' expect_equal(
#'  surv_prob(fp_exp, c(0, 1, 5, 10, 20)),
#'  surv_prob(exp_dist, c(0, 1, 5, 10, 20)),
#'  tolerance = 1e-5
#' )
#'
#' # ---- FP2: valid probabilities ----
#' fp2 <- define_surv_fp(
#'  betas = c(-3, 0.5, -0.2),
#'  powers = c(0, 1)
#' )
#' probs2 <- surv_prob(fp2, c(0, 1, 5, 10, 20))
#' expect_true(all(probs2 >= 0 & probs2 <= 1))
#' expect_true(all(diff(probs2) <= 0))
#'
#' # ---- FP2 repeated powers: valid probabilities ----
#' fp2r <- define_surv_fp(
#'  betas = c(-3, 0.5, -0.2),
#'  powers = c(1, 1)
#' )
#' probs2r <- surv_prob(fp2r, c(0, 1, 5, 10, 20))
#' expect_true(all(probs2r >= 0 & probs2r <= 1))
#' expect_true(all(diff(probs2r) <= 0))
#'
#' # ---- FP1 negative power ----
#' fp_neg <- define_surv_fp(
#'  betas = c(-2, -0.5),
#'  powers = c(-1)
#' )
#' probs_neg <- surv_prob(fp_neg, c(1, 5, 10, 20))
#' expect_true(all(is.finite(probs_neg)))
#' expect_true(all(probs_neg >= 0 & probs_neg <= 1))
#'
#' # ---- Large time value ----
#' fp_large <- define_surv_fp(
#'  betas = c(-3, 0.5),
#'  powers = c(1)
#' )
#' prob_large <- surv_prob(fp_large, 1000)
#' expect_true(!is.nan(prob_large))
#' expect_true(prob_large >= 0)
#'
#' # ---- Scalar time input ----
#' fp_scalar <- define_surv_fp(
#'  betas = c(-3, 0.5),
#'  powers = c(1)
#' )
#' prob_scalar <- surv_prob(fp_scalar, 5)
#' expect_true(length(prob_scalar) == 1)
#' expect_true(prob_scalar >= 0 & prob_scalar <= 1)
#'
#' # ---- S(0) = 1 ----
#' expect_equal(surv_prob(fp_scalar, 0), 1)
#'
#' # ---- Compatibility with apply_hr ----
#' fp_hr <- apply_hr(fp_scalar, 0.5)
#' prob_hr <- surv_prob(fp_hr, c(0, 1, 5))
#' expect_true(all(prob_hr >= 0 & prob_hr <= 1))
#' expect_equal(prob_hr[1], 1)
#'
#' # ---- Compatibility with apply_af ----
#' fp_af <- apply_af(fp_scalar, 2)
#' prob_af <- surv_prob(fp_af, c(0, 1, 5))
#' expect_true(all(prob_af >= 0 & prob_af <= 1))
#' expect_equal(prob_af[1], 1)
#'
#' # ---- Compatibility with apply_or ----
#' fp_or <- apply_or(fp_scalar, 0.5)
#' prob_or <- surv_prob(fp_or, c(0, 1, 5))
#' expect_true(all(prob_or >= 0 & prob_or <= 1))
#' expect_equal(prob_or[1], 1)
#'
#' # ---- Compatibility with apply_shift ----
#' fp_sh <- apply_shift(fp_scalar, 2)
#' prob_sh <- surv_prob(fp_sh, c(0, 1, 5))
#' expect_true(all(prob_sh >= 0 & prob_sh <= 1))
#' expect_equal(prob_sh[1], 1)
#'
#' # ---- Compatibility with apply_td_hr ----
#' fp_td <- apply_td_hr(fp_scalar, time = c(0, 5, 10), hr = c(0.8, 1.0, 1.2))
#' prob_td <- surv_prob(fp_td, c(0, 1, 3))
#' expect_true(all(prob_td >= 0 & prob_td <= 1))
#' expect_equal(prob_td[1], 1)
#'
#' # ---- Compatibility with add_hazards ----
#' exp_ref <- define_surv_param('exp', rate = 0.05)
#' fp_ah <- add_hazards(fp_scalar, exp_ref)
#' prob_ah <- surv_prob(fp_ah, c(0, 1, 5))
#' expect_true(all(prob_ah >= 0 & prob_ah <= 1))
#' expect_equal(prob_ah[1], 1)
#'
#' # ---- Compatibility with mix ----
#' fp_mx <- mix(fp_scalar, 0.5, exp_ref, 0.5)
#' prob_mx <- surv_prob(fp_mx, c(0, 1, 5))
#' expect_true(all(prob_mx >= 0 & prob_mx <= 1))
#' expect_equal(prob_mx[1], 1)
#'
#' # ---- Compatibility with join ----
#' fp_jn <- join(fp_scalar, 5, exp_ref)
#' prob_jn <- surv_prob(fp_jn, c(0, 1, 5))
#' expect_true(all(prob_jn >= 0 & prob_jn <= 1))
#' expect_equal(prob_jn[1], 1)
#'
#' # ---- Compatibility with event_prob ----
#' ep <- event_prob(fp_scalar, 0, 5)
#' expect_true(ep >= 0 & ep <= 1)
#'
#' # ---- Compatibility with plot ----
#' p <- plot(fp_scalar, max_time = 10)
#' expect_true(inherits(p, 'gg'))
#'
surv_prob.surv_fp <- function(x, time, ...) {
    check_times(time, 'calculating survival probabilities', 'time')
    exp(-fp_cumulative_hazard(x$betas, x$powers, time))
}
