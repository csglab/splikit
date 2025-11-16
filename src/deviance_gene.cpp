#include <RcppArmadillo.h>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Helper function to compute deviance for a single row
// Extracted to avoid code duplication between single and multi-threaded paths
inline double compute_row_deviance(const arma::sp_mat& gene_expression, int i, int n_cols) {
  double sum_y  = 0.0;
  double sum_y2 = 0.0;

  // First pass: compute sums
  for (int j = 0; j < n_cols; j++){
    double y = gene_expression(i, j);
    sum_y  += y;
    sum_y2 += y * y;
  }

  // If there are no counts, return 0
  if (sum_y <= 0) {
    return 0.0;
  }

  double mu_hat   = sum_y / n_cols;                            // intercept-only mean
  double variance = (sum_y2 / n_cols) - (mu_hat * mu_hat);     // empirical variance

  // Estimate NB dispersion: theta = mu^2/(var − mu) if var > mu, else (approx. infinite)
  double theta_est = (variance > mu_hat)
    ? (mu_hat * mu_hat) / (variance - mu_hat)
    : 1e12;

  // Second pass: compute deviance
  double dev_row = 0.0;
  for (int j = 0; j < n_cols; j++){
    double y = gene_expression(i, j);
    double term = 0.0;

    if (y > 0)
      term = y * std::log(y / mu_hat);

    term -= (y + theta_est)
      * std::log((y + theta_est) / (mu_hat + theta_est));

    dev_row += 2.0 * term;
  }

  return dev_row;
}

// [[Rcpp::export]]
arma::vec calcNBDeviancesWithThetaEstimation(const arma::sp_mat& gene_expression, int num_threads = 1) {
  try {
    int n_rows = gene_expression.n_rows;
    int n_cols = gene_expression.n_cols;
    arma::vec dev(n_rows, arma::fill::zeros);

    // Notify about OpenMP availability only once (thread-safe)
    static std::atomic<bool> has_printed(false);
    bool expected = false;
    if (has_printed.compare_exchange_strong(expected, true)) {
      #ifdef _OPENMP
        if (num_threads > 1) {
          Rcpp::Rcout << "OpenMP is available. Using " << num_threads << " threads for NB deviance calculation.\n";
        }
      #else
        if (num_threads > 1) {
          Rcpp::Rcout << "OpenMP is not available. Running in single-threaded mode.\n";
        }
      #endif
    }

    #ifdef _OPENMP
    if (num_threads > 1) {
      omp_set_num_threads(num_threads);
      #pragma omp parallel for schedule(dynamic)
      for (int i = 0; i < n_rows; i++){
        dev[i] = compute_row_deviance(gene_expression, i, n_cols);
      }
    } else {
    #endif
      // Single-threaded version
      for (int i = 0; i < n_rows; i++){
        dev[i] = compute_row_deviance(gene_expression, i, n_cols);
      }
    #ifdef _OPENMP
    }
    #endif

    return dev;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("C++ exception in calcNBDeviancesWithThetaEstimation (unknown reason)");
  }
  return arma::vec(); // never reached
}