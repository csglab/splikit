#include <RcppEigen.h>
#include <cmath>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::NumericVector standardizeSparse_variance_loess(const Eigen::SparseMatrix<double>& X,
                                                       bool display_progress = false) {
  // Convert the input matrix to row-major order for efficient row iteration.
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
  SpMat X_row = X;
  
  int nrows = X_row.rows();  // number of features (rows)
  int ncols = X_row.cols();  // number of observations (columns)
  double vmax = std::sqrt(ncols);  // clipping threshold
  
  // Vectors to store each row's mean and variance (computed including zeros)
  Rcpp::NumericVector rowMeans(nrows);
  Rcpp::NumericVector rowVars(nrows);
  
  // Calculate mean and variance for each row.
  for (int i = 0; i < nrows; i++){
    double sum = 0.0;
    int countNonZero = 0;
    for (SpMat::InnerIterator it(X_row, i); it; ++it) {
      sum += it.value();
      countNonZero++;
    }
    int total = ncols;  // total entries (nonzeros + implicit zeros)
    int nZero = total - countNonZero;
    double mean = sum / total;
    rowMeans[i] = mean;
    
    double sqDiffSum = 0.0;
    for (SpMat::InnerIterator it(X_row, i); it; ++it) {
      double diff = it.value() - mean;
      sqDiffSum += diff * diff;
    }
    // Zeros contribute as (0 - mean)^2
    sqDiffSum += nZero * (mean * mean);
    
    double var = (total > 1) ? sqDiffSum / (total - 1) : 0.0;
    rowVars[i] = var;
  }
  
  // Prepare vectors to hold the expected variance and standard deviation.
  // We'll compute these by fitting a loess model to log10(variance) ~ log10(mean)
  // for rows where variance > 0 and mean > 0.
  Rcpp::NumericVector expectedVar(nrows, 0.0);
  Rcpp::NumericVector sd(nrows, 1.0); // default sd
  
  // Identify indices for non-constant rows (variance > 0 and mean > 0 for log transform)
  std::vector<int> nonConstIndices;
  for (int i = 0; i < nrows; i++) {
    if (rowVars[i] > 0 && rowMeans[i] > 0) {
      nonConstIndices.push_back(i);
    }
  }
  
  if(nonConstIndices.size() > 0) {
    int nNonConst = nonConstIndices.size();
    Rcpp::NumericVector loessMean(nNonConst);
    Rcpp::NumericVector loessVar(nNonConst);
    
    for (int j = 0; j < nNonConst; j++){
      int idx = nonConstIndices[j];
      loessMean[j] = rowMeans[idx];
      loessVar[j] = rowVars[idx];
    }
    
    // Create a data frame with columns "mean" and "variance" for the loess fit.
    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["mean"] = loessMean,
                                                   Rcpp::_["variance"] = loessVar);
    
    // Optionally display progress.
    if (display_progress) {
      Rcpp::Rcerr << "Fitting loess to compute expected variance" << std::endl;
    }
    
    // Call R's loess function: fit <- loess(log10(variance) ~ log10(mean), data = df, span = 0.3)
    Rcpp::Function loess("loess");
    Rcpp::List fit = loess(Rcpp::Named("formula") = Rcpp::CharacterVector::create("log10(variance) ~ log10(mean)"),
                           Rcpp::Named("data") = df,
                           Rcpp::Named("span") = 0.3);
    
    // Get the fitted values from the loess model using R's predict function.
    Rcpp::Function predict("predict");
    Rcpp::NumericVector fitted = predict(fit, Rcpp::Named("newdata") = df);
    
    // For each non-constant row, set expected variance as 10^(fitted value)
    // and compute the standard deviation as sqrt(expected variance).
    for (int j = 0; j < nNonConst; j++){
      int idx = nonConstIndices[j];
      expectedVar[idx] = std::pow(10.0, fitted[j]);
      sd[idx] = std::sqrt(expectedVar[idx]);
      if(sd[idx] == 0) sd[idx] = 1.0; // safeguard
    }
  }
  
  // Now, for each row, standardize the values using the computed mean and the loess-derived SD,
  // clip the standardized values at vmax, and compute the variance of these standardized values.
  Rcpp::NumericVector result(nrows);
  for (int i = 0; i < nrows; i++){
    double mean = rowMeans[i];
    double curr_sd = sd[i];
    int total = ncols;
    
    double sumSquares = 0.0;
    int countNonZero = 0;
    for (SpMat::InnerIterator it(X_row, i); it; ++it) {
      double standardized_val = (it.value() - mean) / curr_sd;
      double clipped_val = std::min(vmax, standardized_val);
      sumSquares += clipped_val * clipped_val;
      countNonZero++;
    }
    int nZero = total - countNonZero;
    double zero_std = (0 - mean) / curr_sd;
    sumSquares += nZero * (zero_std * zero_std);
    
    result[i] = (total > 1) ? sumSquares / (total - 1) : 0.0;
  }
  
  return result;
}
