#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Fast Computation of Row Variances for a Sparse Matrix
//' 
//' This function computes the variance for each row of a sparse matrix using a
//' single pass over the nonzero elements by iterating in column-major order.
//'
//' @param M A sparse matrix (dgCMatrix) of gene expression counts.
//' @return A numeric vector with the variance for each row.
//' @examples
//' library(Matrix)
//' set.seed(42)
//' M <- rsparsematrix(1000, 500, density = 0.05) * 10
//' row_variances <- rowVarsSparseFast(M)
// [[Rcpp::export]]
arma::vec get_row_variance(const arma::sp_mat& M) {
  int n_rows = M.n_rows;
  int n_cols = M.n_cols;
  
  // Preallocate vectors for row sums and row squared sums.
  arma::vec row_sum(n_rows, arma::fill::zeros);
  arma::vec row_sum2(n_rows, arma::fill::zeros);
  
  // Loop over columns (efficient for CSC storage).
  for (arma::sp_mat::const_iterator it = M.begin();
       it != M.end(); ++it) {
    int i = it.row();  // row index of the nonzero element
    double value = *it;
    row_sum[i] += value;
    row_sum2[i] += value * value;
  }
  
  // Compute the variance for each row.
  arma::vec variances(n_rows, arma::fill::zeros);
  for (int i = 0; i < n_rows; i++) {
    double mean = row_sum[i] / n_cols;
    variances[i] = row_sum2[i] / n_cols - mean * mean;
  }
  
  return variances;
}
