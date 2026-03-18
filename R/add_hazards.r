#' Add Hazards
#' 
#' Combine two or more survival distributions as independent risks
#' by adding their hazards.
#' 
#' @name add_hazards
#' @rdname add_hazards
#' @export
#' 
#' @param dist1 survival distribution to add
#' @param dist2 second survival distribution to add
#' @param ... additional survival distributions to add
#'   
#' @return a `surv_add_haz` object
#' 
#' @examples
#' 
#' dist1 <- define_surv_param("exp", rate = .125)
#' dist2 <- define_surv_param("weibull", shape = 1.2, scale = 50)
#' combined_dist <- add_hazards(dist1, dist2)
#' 
#' @tests
#' 
#' dist1 <- define_surv_param("exp", rate = .125)
#' dist2 <- define_surv_param("weibull", shape = 1.2, scale = 50)
#' dist3 <- define_surv_param("weibull", shape = 1.1, scale = 30)
#' expect_equal(
#'  add_hazards(dist1, dist2, dist3),
#'  create_list_object(
#'      c('surv_add_haz', 'surv_combined', 'surv_dist'),
#'      dists = list(dist1, dist2, dist3)
#'  )
#' )
#' expect_error(
#'  add_hazards(dist1, dist2, 'foo'),
#'  'Error adding hazards, invalid survival distribution provided.',
#'  fixed = TRUE
#' )
#' 
add_hazards <- function(dist1, dist2, ...) {
  
    # Compile and check distributions
    dists <- append(list(dist1, dist2), list(...))
    walk(dists, function(x) {
        valid <- is_surv_dist(x)
        if (!valid) {
            err <- get_and_populate_message('add_hazards_wrong_type_dist')
            stop(err, call. = show_call_error())
        }
    })

    create_list_object(
        c('surv_add_haz', 'surv_combined', 'surv_dist'),
        dists = dists
    )
}

#' @export
#' 
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.025)
#' dist2 <- define_surv_param('exp', rate = 0.010)
#' dist3 <- define_surv_param('exp', rate = 0.002)
#' dist4 <- define_surv_param('exp', rate = 0.035)
#' dist5 <- define_surv_param('exp', rate = 0.037)
#' expect_equal(
#'  surv_prob(dist4, seq_len(100)),
#'  surv_prob(add_hazards(dist1, dist2), seq_len(100))
#' )
#' expect_equal(
#'  surv_prob(dist5, seq_len(100)),
#'  surv_prob(add_hazards(dist1, dist2, dist3), seq_len(100))
#' )
surv_prob.surv_add_haz <- function(x, time, ...) {
    check_times(time, 'calculating survival probabilities', 'time')
    Reduce(`*`, map(x$dists, function(dist) surv_prob(dist, time, ...)))
}

#' @export
#'
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.05)
#' dist2 <- define_surv_param('exp', rate = 0.10)
#' combined <- add_hazards(dist1, dist2)
#'
#' # Roundtrip
#' probs <- c(0.9, 0.5, 0.1)
#' times <- surv_quantile(combined, probs)
#' expect_equal(surv_prob(combined, times), probs, tolerance = 1e-6)
#'
#' # Edge cases
#' expect_equal(surv_quantile(combined, 1), 0)
#' expect_equal(surv_quantile(combined, 0), Inf)
#'
surv_quantile.surv_add_haz <- function(x, probs, ...) {
    check_probs(probs)
    n <- length(x$dists)
    vapply(probs, function(p) {
        if (p == 1) return(0)
        if (p == 0) return(Inf)

        # Bracket: [min Q_i(p^(1/n)), max Q_i(p^(1/n))]
        p_adj <- p^(1/n)
        component_quantiles <- map_dbl(x$dists, function(d) surv_quantile(d, p_adj, ...))
        lo <- min(component_quantiles)
        hi <- max(component_quantiles)

        if (abs(hi - lo) < 1e-12) return(lo)

        finite_qs <- component_quantiles[is.finite(component_quantiles)]
        if (length(finite_qs) == 0) return(Inf)
        if (any(is.infinite(component_quantiles))) {
            hi <- max(finite_qs) * 2
            while (surv_prob(x, hi) > p) hi <- hi * 2
        }

        uniroot(function(t) surv_prob(x, t) - p, c(lo, hi),
                tol = .Machine$double.eps^0.5)$root
    }, numeric(1))
}

#' @export 
#' 
#' @tests
#' dist1 <- define_surv_param('exp', rate = 0.12)
#' dist2 <- define_surv_param('exp', rate = 0.18)
#' expect_output(
#'  print(add_hazards(dist1, dist2)),
#'  'A survival distribution combining the hazards of:
#'   * Distribution 1: An exponential distribution (rate = 0.12).
#'   * Distribution 2: An exponential distribution (rate = 0.18).',
#'  fixed = TRUE
#' )
print.surv_add_haz <- function(x, ...) {
    args <- list(...)
    items_lines <- map_chr(seq_along(x$dists), function(i) {
        dist_output <- to_list_item_output(x$dists[[i]])
        glue('    * Distribution {i}: {dist_output}')
    })
    output <- paste0(c('A survival distribution combining the hazards of:', items_lines), collapse = '\n')

    cat(output)
}