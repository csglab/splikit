#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec calcNBDeviancesWithThetaEstimation(const arma::sp_mat& gene_expression, int num_threads = 1) {

  int n_rows = gene_expression.n_rows;
  int n_cols = gene_expression.n_cols;  // number of cells
  arma::vec dev(n_rows, arma::fill::zeros);

  // Notify about OpenMP availability only once
  static bool has_printed = false;
  if (!has_printed) {
    #ifdef _OPENMP
      if (num_threads > 1) {
        Rcpp::Rcout << "OpenMP is available. Using " << num_threads << " threads for NB deviance calculation.\n";
      }
    #else
      if (num_threads > 1) {
        Rcpp::Rcout << "OpenMP is not available. Running in single-threaded mode.\n";
      }
    #endif
    has_printed = true;
  }

  #ifdef _OPENMP
  if (num_threads > 1) {
    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_rows; i++){
      double sum_y  = 0.0;
      double sum_y2 = 0.0;

      for (int j = 0; j < n_cols; j++){
        double y = gene_expression(i, j);
        sum_y  += y;        // compute total counts
        sum_y2 += y * y;    // and total of squares for gene i.
      }

      // If there are no counts, set deviance to 0.
      if (sum_y <= 0) {
        dev[i] = 0;
        continue;
      }

      double mu_hat   = sum_y / n_cols;                            // intercept-only mean
      double variance = (sum_y2 / n_cols) - (mu_hat * mu_hat);     // empirical variance

      // Estimate NB dispersion: theta = mu^2/(var − mu) if var > mu, else (approx. infinite)
      double theta_est = (variance > mu_hat)
        ? (mu_hat * mu_hat) / (variance - mu_hat)
        : 1e12;

      double dev_row = 0.0;
      // Second pass: deviance contribution from each cell
      for (int j = 0; j < n_cols; j++){
        double y = gene_expression(i, j);
        double term = 0.0;

        if (y > 0)
          term = y * std::log(y / mu_hat);

        term -= (y + theta_est)
          * std::log((y + theta_est) / (mu_hat + theta_est));

        dev_row += 2.0 * term;
      }

      dev[i] = dev_row;
    }
  } else {
  #endif
    // Single-threaded version
    for (int i = 0; i < n_rows; i++){
      double sum_y  = 0.0;
      double sum_y2 = 0.0;

      for (int j = 0; j < n_cols; j++){
        double y = gene_expression(i, j);
        sum_y  += y;        // compute total counts
        sum_y2 += y * y;    // and total of squares for gene i.
      }

      // If there are no counts, set deviance to 0.
      if (sum_y <= 0) {
        dev[i] = 0;
        continue;
      }

      double mu_hat   = sum_y / n_cols;                            // intercept-only mean
      double variance = (sum_y2 / n_cols) - (mu_hat * mu_hat);     // empirical variance

      // Estimate NB dispersion: theta = mu^2/(var − mu) if var > mu, else (approx. infinite)
      double theta_est = (variance > mu_hat)
        ? (mu_hat * mu_hat) / (variance - mu_hat)
        : 1e12;

      double dev_row = 0.0;
      // Second pass: deviance contribution from each cell
      for (int j = 0; j < n_cols; j++){
        double y = gene_expression(i, j);
        double term = 0.0;

        if (y > 0)
          term = y * std::log(y / mu_hat);

        term -= (y + theta_est)
          * std::log((y + theta_est) / (mu_hat + theta_est));

        dev_row += 2.0 * term;
      }

      dev[i] = dev_row;
    }
  #ifdef _OPENMP
  }
  #endif

  return dev;
}