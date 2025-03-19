// calcNBDeviancesWithThetaEstimation.cpp
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Calculate Deviances for an Intercept-Only Negative Binomial GLM on Sparse Matrices,
//' estimating theta from the relationship: var = mean + mean^2/theta.
//'
//' @param gene_expression A sparse matrix (dgCMatrix) of gene expression counts.
//' @return A numeric vector of deviances, one per row (gene).
//' @details For each gene (row), the function computes the fitted mean (total counts divided by the number of cells)
//' and estimates the dispersion parameter \(\theta\) using
//' \deqn{ \theta = \frac{\mu^2}{\text{Var}(Y) - \mu}, }
//' where the variance is estimated as the average of the squared counts minus the square of the mean.
//' If \(\text{Var}(Y) \le \mu\) (i.e. no overdispersion), theta is set to a large number to approximate a Poisson model.
//' The per-cell negative binomial deviance is then computed as
//' \deqn{ d(y, \mu) = 2\left[y \log\left(\frac{y}{\mu}\right) - (y+\theta) \log\left(\frac{y+\theta}{\mu+\theta}\right) \right]. }
//' @examples
//' dev_vals <- calcNBDeviancesWithThetaEstimation(gene_expr)
// [[Rcpp::export]]
arma::vec calcNBDeviancesWithThetaEstimation(const arma::sp_mat& gene_expression) {
  int n_rows = gene_expression.n_rows;
  int n_cols = gene_expression.n_cols;  // number of cells
  arma::vec dev(n_rows, arma::fill::zeros);
  
  for (int i = 0; i < n_rows; i++){
    double sum_y = 0.0;
    double sum_y2 = 0.0;
    
   
    for (int j = 0; j < n_cols; j++){
      double y = gene_expression(i, j);
      sum_y += y;  // compute total counts 
      sum_y2 += y * y; // and total of squares for gene i.
    }
    
    // If there are no counts, set deviance to 0.
    if (sum_y <= 0) {
      dev[i] = 0;
      continue;
    }
    
   
    double mu_hat = sum_y / n_cols;  // Fitted mean under the intercept-only model (MLE).
    
    // Estimate the variance as (mean of squares) - (square of the mean).
    double variance = (sum_y2 / n_cols) - (mu_hat * mu_hat);
    
    // Estimate theta using: theta = mu_hat^2 / (variance - mu_hat) if variance > mu_hat.
    double theta_est;
    if (variance > mu_hat) {
      theta_est = (mu_hat * mu_hat) / (variance - mu_hat);
    } else {
      // If variance is not greater than mean, approximate Poisson behavior.
      theta_est = 1e12;  // large number to approximate theta -> infinity
    }
    
    double dev_row = 0.0;
    // Second pass: compute deviance contribution from each cell.
    for (int j = 0; j < n_cols; j++){
      double y = gene_expression(i, j);
      double term = 0.0;
      
      // Add y * log(y/mu_hat); define 0*log(0) as 0.
      if (y > 0)
        term = y * std::log(y / mu_hat);
      
      // Subtract (y + theta_est)*log((y + theta_est)/(mu_hat + theta_est)).
      term -= (y + theta_est) * std::log((y + theta_est) / (mu_hat + theta_est));
      
      // Multiply by 2 as per the deviance definition.
      dev_row += 2 * term;
    }
    
    dev[i] = dev_row;
  }
  
  return dev;
}
