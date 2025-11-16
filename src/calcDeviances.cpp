// calcDeviances_ratio.cpp
#include <RcppArmadillo.h>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec calcDeviances_ratio(const arma::sp_mat& M1,
                              const arma::sp_mat& M2,
                              int num_threads = 1) {
  try {
    int n_rows = M1.n_rows;
    int n_cols = M1.n_cols;
    arma::vec dev(n_rows, arma::fill::zeros);

    // Notify about OpenMP availability only once (thread-safe)
    // Improved to reduce message spam
    static std::atomic<bool> has_printed(false);
    bool expected = false;
    if (has_printed.compare_exchange_strong(expected, true)) {
#ifdef _OPENMP
      if (num_threads > 1) {
        Rcpp::Rcout << "OpenMP detected: using " << num_threads << " threads for deviance calculation.\n";
      }
#else
      // Only warn if user requested multi-threading but it's not available
      if (num_threads > 1) {
        Rcpp::Rcout << "OpenMP not available: running in single-threaded mode despite num_threads="
                    << num_threads << " request.\n";
      }
#endif
    }

    // Set the number of threads if OpenMP is available
#ifdef _OPENMP
    if (num_threads > 1) {
      omp_set_num_threads(num_threads);
    }
#endif

    // Parallelize the outer loop if OpenMP is available and num_threads > 1
#ifdef _OPENMP
#pragma omp parallel for if(num_threads > 1)
#endif
    for (int i = 0; i < n_rows; i++) {
      double sum_y = 0.0, sum_n = 0.0;
      for (int k = 0; k < n_cols; k++) {
        double y = M1(i,k), f = M2(i,k);
        sum_y += y;
        sum_n += (y + f);
      }
      if (sum_n <= 0) { dev[i] = 0; continue; }
      double p_hat = sum_y / sum_n;
      if (p_hat <= 0 || p_hat >= 1) { dev[i] = 0; continue; }

      double dev_row = 0.0;
      for (int k = 0; k < n_cols; k++) {
        double y = M1(i,k), f = M2(i,k);
        double n_i = y + f;
        if (n_i <= 0) continue;
        if (y > 0)
          dev_row += 2 * y * std::log(y / (n_i * p_hat));
        if (n_i - y > 0)
          dev_row += 2 * (n_i - y) * std::log((n_i - y) / (n_i * (1 - p_hat)));
      }
      dev[i] = dev_row;
    }

    return dev;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("C++ exception in calcDeviances_ratio (unknown reason)");
  }
  return arma::vec(); // never reached
}
