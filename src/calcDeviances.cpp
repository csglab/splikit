// calcDeviances.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Calculate Deviances for an Intercept-Only Binomial GLM on Sparse Matrices
 //'
 //' @param M1 A sparse matrix (dgCMatrix) of “success” counts.
 //' @param M2 A sparse matrix (dgCMatrix) of “failure” counts.
 //' @return A numeric vector of deviances, one per row (junction event).
 //' @details For each row, the function computes the MLE \(\hat{p}\) over all cells
 //' and then sums the binomial deviance contributions over all cells.
 //' If the row has zero trials or if \(\hat{p}\) is 0 or 1, the deviance is set to 0.
 //' @examples
 //' # Assuming M1_sub and M2_sub have been subset for a given library:
 //' dev_vals <- calcDeviances(M1_sub, M2_sub)
 //' 
 // [[Rcpp::export]]
 arma::vec calcDeviances(const arma::sp_mat& M1, const arma::sp_mat& M2) {
   int n_rows = M1.n_rows;
   int n_cols = M1.n_cols;  // number of cells in the library subset
   arma::vec dev(n_rows, arma::fill::zeros);
   
   for (int i = 0; i < n_rows; i++){
     double sum_y = 0.0;
     double sum_n = 0.0;
     
     // First pass: compute total successes and total trials for row i.
     for (int k = 0; k < n_cols; k++){
       double y_val = M1(i, k); // Access element (returns 0 if not stored)
       double f_val = M2(i, k);
       double n_i = y_val + f_val;
       sum_y += y_val;
       sum_n += n_i;
     }
     
     // If there are no trials, return 0 deviance for this row.
     if (sum_n <= 0) {
       dev[i] = 0;
       continue;
     }
     
     double p_hat = sum_y / sum_n;
     // If p_hat is on the boundary, set deviance to 0.
     if (p_hat <= 0 || p_hat >= 1) {
       dev[i] = 0;
       continue;
     }
     
     // Second pass: compute the deviance contribution from each cell.
     double dev_row = 0.0;
     for (int k = 0; k < n_cols; k++){
       double y_val = M1(i, k);
       double f_val = M2(i, k);
       double n_i = y_val + f_val;
       if (n_i <= 0) continue; // no trials in this cell
       
       double term = 0.0;
       if (y_val > 0)
         term += y_val * std::log(y_val / (n_i * p_hat));
       if ((n_i - y_val) > 0)
         term += (n_i - y_val) * std::log((n_i - y_val) / (n_i * (1 - p_hat)));
       
       dev_row += 2 * term;
     }
     
     dev[i] = dev_row;
   }
   
   return dev;
 }
 