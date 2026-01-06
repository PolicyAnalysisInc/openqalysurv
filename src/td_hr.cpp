#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_td_surv_prob(
    NumericVector query_time,      // Query times
    NumericVector stored_time,     // Stored time points (sorted)
    NumericVector hr,              // HR for each stored time point
    NumericVector log_surv_stored, // log(S0) at stored times
    NumericVector log_surv_query   // log(S0) at query times
) {
    int n_query = query_time.size();
    int n_stored = stored_time.size();
    NumericVector result(n_query);

    // Precompute log differences between consecutive stored times
    // log_diff[i] = log(S0(stored_time[i+1])) - log(S0(stored_time[i]))
    NumericVector log_diff(n_stored - 1);
    for (int i = 0; i < n_stored - 1; i++) {
        log_diff[i] = log_surv_stored[i + 1] - log_surv_stored[i];
    }

    // Precompute cumulative weighted log-survival at each stored time point
    // cum_weighted[i] = sum_{k=0}^{i-1} hr[k] * log_diff[k]
    NumericVector cum_weighted(n_stored);
    cum_weighted[0] = 0.0;
    for (int i = 0; i < n_stored - 1; i++) {
        cum_weighted[i + 1] = cum_weighted[i] + hr[i] * log_diff[i];
    }

    // Compute survival for each query time
    for (int j = 0; j < n_query; j++) {
        double t = query_time[j];

        // Handle t = 0 specially
        if (t == 0.0) {
            result[j] = 1.0;
            continue;
        }

        // Find which interval t falls into using binary search
        // Find largest i such that stored_time[i] <= t
        int lo = 0, hi = n_stored - 1;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (stored_time[mid] <= t) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        int interval_idx = lo;

        // Contribution from complete intervals up to interval_idx
        double cum_log_surv = cum_weighted[interval_idx];

        // Partial interval contribution (if t is between stored points)
        if (t > stored_time[interval_idx] && interval_idx < n_stored - 1) {
            double partial_log_diff = log_surv_query[j] - log_surv_stored[interval_idx];
            cum_log_surv += hr[interval_idx] * partial_log_diff;
        }

        result[j] = exp(cum_log_surv);
    }

    return result;
}
