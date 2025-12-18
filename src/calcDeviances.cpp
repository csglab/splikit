// calcDeviances_ratio.cpp
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec calcDeviances_ratio(const arma::sp_mat& M1,
                              const arma::sp_mat& M2,
                              int num_threads = 1) {
  int n_rows = M1.n_rows;
  int n_cols = M1.n_cols;
  arma::vec dev(n_rows, arma::fill::zeros);
  
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
}
