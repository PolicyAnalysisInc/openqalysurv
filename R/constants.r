# Error/Warning messages
messages <- list(
    missing_parameters = 'Error defining {dist} distribution, parameters missing from function call: {params}.',
    invalid_theta = 'Error defining cure model, cure fraction (theta) must be in range [0-1].',
    truncated_vector = 'Parameter {param} was length > 1 and only the first element will be used.',
    n_spline_params = 'Error defining restricted cubic spline distribution, must provide at least two parameter values followed by a matching number of knot times.',
    spline_param_type = 'Error defining restricted cubic spline distribution, parameter was of type "{class}" instead of "numeric".',
    spline_param_names = 'Error defining restricted cubic spline distribution, incorrect argument names were provided.',
    model_no_covariates = 'Generating prediction from model with covariates but no covariates were provided. Predictions will reflect weighted average of predictions for subjects used to fit model.',
    life_table_missing_gender_mix = 'Error defining life-table, must provide either "percent_male", or "percent_female", but not both.',
    life_table_missing_columns = 'Error defining life-table, the following columns were expected but not found: {missing_cols}.',
    life_table_dupe_age = 'Error defining life-table, column "{age_col}" contained duplicate values.',
    life_table_varying_bands = 'Error defining life-table, life-table must use constant age bands.',
    km_missing_columns = 'Error defining KM, the following columns were expected but not found: {missing_cols}.',
    km_dupe_time = 'Error defining KM, column "{time_col}" contained duplicate values.',
    km_invalid_start = 'Error defining KM, column "{time_col}" must start with a value 0 and "{surv_col}" must start with a value 1.',
    km_increasing_surv = 'Error defining KM, column "{surv_col}" may not be increasing with respect to "{time_col}".',
    km_invalid_prob = 'Error defining KM, values in column "{surv_col}" must be within the interval [0-1].',
    km_invalid_types = 'Error defining KM, the following columns were of invalid type: {invalid_type_cols}.',
    km_missing_values = 'Error defining KM, the following columns contained missing values: {missing_value_cols}.',
    surv_prob_wrong_type = 'Error calculating survival probabilities, invalid survival distribution provided.',
    event_prob_wrong_type = 'Error calculating event probabilities, invalid survival distribution provided.',
    apply_hr_wrong_type_dist = 'Error applying hazard ratio, invalid survival distribution provided.',
    apply_hr_wrong_type_hr = 'Error applying hazard ratio, "hr" must be numeric.',
    apply_hr_missing_hr = 'Error applying hazard ratio, "hr" cannot be NA.',
    apply_hr_invalid_hr = 'Error applying hazard ratio, "hr" cannot be negative.',
    apply_td_hr_wrong_type_dist = 'Error applying time-dependent hazard ratio, invalid survival distribution provided.',
    apply_td_hr_wrong_type_hr = 'Error applying time-dependent hazard ratio, "hr" must be numeric.',
    apply_td_hr_missing_hr = 'Error applying time-dependent hazard ratio, "hr" cannot contain NA values.',
    apply_td_hr_invalid_hr = 'Error applying time-dependent hazard ratio, "hr" cannot contain negative values.',
    apply_td_hr_length_mismatch = 'Error applying time-dependent hazard ratio, "time" and "hr" must have the same length.',
    apply_td_hr_conflicting_hr = 'Error applying time-dependent hazard ratio, conflicting hazard ratios found for the same time value.',
    apply_td_hr_time_exceeds_stored = 'Error calculating survival probabilities, query time exceeds maximum stored time.',
    apply_af_wrong_type_dist = 'Error applying acceleration factor, invalid survival distribution provided.',
    apply_af_wrong_type_af = 'Error applying acceleration factor, "af" must be numeric.',
    apply_af_missing_af = 'Error applying acceleration factor, "af" cannot be NA.',
    apply_af_invalid_af = 'Error applying acceleration factor, "af" cannot be negative.',
    apply_or_wrong_type_dist = 'Error applying odds ratio, invalid survival distribution provided.',
    apply_or_wrong_type_or = 'Error applying odds ratio, "or" must be numeric.',
    apply_or_missing_or = 'Error applying odds ratio, "or" cannot be NA.',
    apply_or_invalid_or = 'Error applying odds ratio, "or" cannot be negative.',
    apply_shift_wrong_type_dist = 'Error applying shift, invalid survival distribution provided.',
    apply_shift_wrong_type_shift = 'Error applying shift, "shift" must be numeric.',
    apply_shift_missing_shift = 'Error applying shift, "shift" cannot be NA.',
    join_wrong_type_cut = 'Error joining distributions, cuts times must be numeric.',
    join_missing_cut = 'Error joining distributions, cuts times cannot be NA.',
    join_invalid_cut = 'Error joining distributions, cut times cannot be negative.',
    join_wrong_n_args = 'Error joining distributions, must provide an odd number of arguments corresponding to n distributions and n - 1 cut points.',
    join_wrong_type_dist = 'Error joining distributions, invalid survival distribution provided.',
    join_cuts_order = 'Error joining distributions, distributions and cutpoints must be provided in order.',
    set_covariates_wrong_type_dist = 'Error setting covariates, only survfit and flexsurv models are supported.',
    set_covariates_wrong_type_data = 'Error setting covariates, "data" must be provided as a data frame.',
    set_covariates_missing_data = 'Error setting covariates, must provide either "data" or named arguments for covariate values.',
    mix_wrong_type_weight = 'Error mixing distributions, weights must be numeric.',
    mix_missing_weight = 'Error mixing distributions, weights cannot be NA.',
    mix_invalid_weight = 'Error mixing distributions, weights cannot be negative.',
    mix_wrong_n_args = 'Error mixing distributions, must provide an even number of arguments corresponding to n distributions and weights.',
    mix_wrong_type_dist = 'Error mixing distributions, invalid survival distribution provided.',
    mix_weights_wrong_sum = 'Error mixing distributions, weights must sum to 1.',
    add_hazards_wrong_type_dist = 'Error adding hazards, invalid survival distribution provided.',
    check_time_wrong_class = 'Error {context}, "{time_name}" must be numeric.',
    check_time_negative = 'Error {context}, "{time_name}" cannot be negative.',
    check_time_missing = 'Error {context}, "{time_name}" cannot be NA.',
    lifetable_mortality_wrong_type = 'Error defining life-table, death rates must be numeric.',
    fp_betas_wrong_type = 'Error defining fractional polynomial distribution, "betas" must be numeric.',
    fp_powers_wrong_type = 'Error defining fractional polynomial distribution, "powers" must be numeric.',
    fp_betas_missing = 'Error defining fractional polynomial distribution, "betas" cannot contain NA values.',
    fp_powers_missing = 'Error defining fractional polynomial distribution, "powers" cannot contain NA values.',
    fp_length_mismatch = 'Error defining fractional polynomial distribution, "betas" must have length equal to length of "powers" plus one.',
    fp_min_powers = 'Error defining fractional polynomial distribution, must specify at least one power.',
    invalid_source_dist = 'Error defining {source} {dist} distribution, "{dist}" is not a supported distribution for {source}.',
    missing_source_params = 'Error defining {source} {dist} distribution, required parameters missing: {params}.',
    invalid_source_param = 'Error defining {source} {dist} distribution, parameter "{param}" must be a single non-NA numeric value.',
    surv_quantile_invalid_probs = 'Error computing survival quantiles, "probs" must be numeric values in the interval [0, 1].',
    surv_quantile_missing_probs = 'Error computing survival quantiles, "probs" cannot contain NA values.',
    surv_random_invalid_n = 'Error generating random survival times, "n" must be a single positive integer.',
    max_hazards_wrong_type_dist = 'Error computing maximum hazards, invalid survival distribution provided.',
    max_hazards_invalid_cycle_length = 'Error computing maximum hazards, "cycle_length" must be a single positive numeric value.'
)

# Possible values for distribution argument to flexsurvreg
flexsurv_dists <- c("exp", "weibull", "weibullPH", "lnorm", "llogis", "gamma", "gompertz", "gengamma", "gengamma.orig", "genf", "genf.orig")

# Display names for flexsurv distributions
flexsurv_dist_aliases <- list(
    exp = 'exponential',
    weibull = 'Weibull (AFT)',
    weibullPH = 'Weibull (PH)',
    llogis = 'log-logistic',
    lnorm = 'lognormal',
    gamma = 'gamma',
    gengamma = 'generalized gamma (stable)',
    gengamma.org = 'generalized gamma (original)',
    genf = 'generalized F (stable)',
    genf.orif = 'generalized F (original)',
    gompertz = 'Gompertz'
)

# Possible values for scale argument to flexsurvspline
flexsurv_spline_scales <- c("hazard", "odds", "normal")

flexsurv_spline_scale_aliases <- list(
    hazard = 'log cumulative hazard',
    odds = 'log cumulative odds',
    normal = 'inverse normal CDF'
)

# Default values for options
default_options <- list(
    openqalysurv.show_call_signature_in_errors = FALSE,
    openqalysurv.show_call_signature_in_warnings = FALSE
)

# Use 'an' instead of 'a' when a word starts with one of
# these letters
word_start_vowels <- c('a','e','i','o','u')

# ============================================================================
# Software-specific distribution mappings
# ============================================================================

# survreg (R survival package) -> flexsurv distribution names
survreg_dists <- c(
    exponential = "exp",
    weibull = "weibull",
    lognormal = "lnorm",
    loglogistic = "llogis"
)

# SAS PROC LIFEREG -> flexsurv distribution names
lifereg_dists <- c(
    exponential = "exp",
    weibull = "weibull",
    lnormal = "lnorm",
    llogistic = "llogis",
    gamma = "gengamma"
)

# Stata streg AFT metric -> flexsurv distribution names
streg_aft_dists <- c(
    exponential = "exp",
    weibull = "weibull",
    lognormal = "lnorm",
    loglogistic = "llogis",
    ggamma = "gengamma"
)

# Stata streg PH metric -> flexsurv distribution names
streg_ph_dists <- c(
    exponential = "exp",
    weibull = "weibullPH",
    gompertz = "gompertz"
)

# ============================================================================
# Expected parameter names per software/distribution
# ============================================================================

survreg_dist_params <- list(
    exponential = "intercept",
    weibull = c("intercept", "scale"),
    lognormal = c("intercept", "scale"),
    loglogistic = c("intercept", "scale")
)

lifereg_dist_params <- list(
    exponential = "intercept",
    weibull = c("intercept", "scale"),
    lnormal = c("intercept", "scale"),
    llogistic = c("intercept", "scale"),
    gamma = c("intercept", "scale", "shape")
)

streg_aft_dist_params <- list(
    exponential = "cons",
    weibull = c("cons", "ln_p"),
    lognormal = c("cons", "sigma"),
    loglogistic = c("cons", "gamma"),
    ggamma = c("cons", "sigma", "kappa")
)

streg_ph_dist_params <- list(
    exponential = "cons",
    weibull = c("cons", "ln_p"),
    gompertz = c("cons", "gamma")
)

# ============================================================================
# Shared helpers for software-specific distribution functions
# ============================================================================

#' @tests
#' expect_equal(match_dist_name("Weibull", survreg_dists, "survreg"), "weibull")
#' expect_equal(match_dist_name("EXPONENTIAL", survreg_dists, "survreg"), "exponential")
#' expect_error(
#'  match_dist_name("gamma", survreg_dists, "survreg"),
#'  'is not a supported distribution',
#'  fixed = TRUE
#' )
match_dist_name <- function(distribution, valid_dists, source) {
    dist_lower <- tolower(distribution)
    matched <- tryCatch(
        match.arg(dist_lower, choices = names(valid_dists)),
        error = function(e) {
            err <- get_and_populate_message(
                'invalid_source_dist',
                source = source,
                dist = distribution
            )
            stop(err, call. = show_call_error())
        }
    )
    matched
}

#' @tests
#' expect_silent(
#'  check_source_params(
#'      list(intercept = 3, scale = 0.8),
#'      c("intercept", "scale"),
#'      "survreg",
#'      "weibull"
#'  )
#' )
#' expect_error(
#'  check_source_params(
#'      list(intercept = 3),
#'      c("intercept", "scale"),
#'      "survreg",
#'      "weibull"
#'  ),
#'  'required parameters missing',
#'  fixed = TRUE
#' )
check_source_params <- function(args, required, source, dist) {
    missing_params <- required[!required %in% names(args)]
    if (length(missing_params) > 0) {
        err <- get_and_populate_message(
            'missing_source_params',
            source = source,
            dist = dist,
            params = quoted_list_string(missing_params)
        )
        stop(err, call. = show_call_error())
    }
}

#' @tests
#' expect_silent(check_numeric_param(3.5, "intercept", "survreg", "weibull"))
#' expect_error(
#'  check_numeric_param("a", "intercept", "survreg", "weibull"),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
#' expect_error(
#'  check_numeric_param(NA, "intercept", "survreg", "weibull"),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
#' expect_error(
#'  check_numeric_param(c(1, 2), "intercept", "survreg", "weibull"),
#'  'must be a single non-NA numeric value',
#'  fixed = TRUE
#' )
check_numeric_param <- function(value, name, source, dist) {
    if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
        err <- get_and_populate_message(
            'invalid_source_param',
            source = source,
            dist = dist,
            param = name
        )
        stop(err, call. = show_call_error())
    }
}