// Implements the variance-stabilizing transformation (VST) step of Hafemeister & Satija (2019)
// without requiring Seurat. Operates on sparse dgCMatrix input using Rcpp + S4.
//
// Hafemeister, C. & Satija, R. (2019).
// "Normalization and variance stabilization of single-cell RNA-seq data using
// regularized negative binomial regression." Cell, 176(4), 1375–1389.e4.
// DOI: https://doi.org/10.1016/j.cell.2019.05.031
//
// This file provides a lean Rcpp implementation that:
//  - Extracts slots (i, p, x, Dim) from a dgCMatrix S4 object
//  - Computes per-row means and variances (including zero entries)
//  - Fits a loess curve to log10(var) ~ log10(mean)
//  - Returns standardized variances (squared z-scores) per row
//  - Logs progress via Rcpp::Rcout
//
// Only 32-bit indexing is supported, consistent with R's native dgCMatrix.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector standardizeSparse_variance_vst(SEXP matSEXP,
                                               bool display_progress = false) {
  Rcout << "[SSVL] entering\n";

  // 1) validate & extract
  if (!Rf_isS4(matSEXP) || !Rf_inherits(matSEXP, "dgCMatrix"))
    stop("`mat` must be a dgCMatrix");
  S4 M(matSEXP);
  IntegerVector dim = M.slot("Dim");
  IntegerVector i   = M.slot("i");
  IntegerVector p   = M.slot("p");
  NumericVector x   = M.slot("x");
  int nrow = dim[0], ncol = dim[1];
  if (display_progress) Rcout << "[SSVL] matrix is " << nrow << "×" << ncol << "\n";

  // 2) first pass: sums, sums2, nnzCount per row
  std::vector<double> sum(nrow, 0.0), sum2(nrow, 0.0);
  std::vector<int>    nnz(nrow, 0);
  for (int col = 0; col < ncol; ++col) {
    for (int idx = p[col]; idx < p[col+1]; ++idx) {
      int r = i[idx];
      double v = x[idx];
      sum[r]  += v;
      sum2[r] += v * v;
      nnz[r]  += 1;
    }
  }

  // 3) compute raw means and variances (including zeros)
  NumericVector rowMean(nrow), rowVar(nrow);
  for (int r = 0; r < nrow; ++r) {
    double  mu = sum[r] / ncol;
    // zeros contribute 0 to sum2
    double mean_of_squares = sum2[r] / ncol;
    rowMean[r] =  mu;
    rowVar [r] = mean_of_squares -  mu* mu;
  }
  if (display_progress) Rcout << "[SSVL] computed raw means & vars\n";

  // 4) fit loess on log10(var) ~ log10(mean) for those with var>0
  std::vector<int> good;
  good.reserve(nrow);
  for (int r = 0; r < nrow; ++r)
    if (rowVar[r] > 0 && rowMean[r] > 0)
      good.push_back(r);
    NumericVector expectedVar(nrow, 0.0), sd(nrow, 1.0);

    if (!good.empty()) {
      int k = good.size();
      NumericVector lm(k), lv(k);
      for (int j = 0; j < k; ++j) {
        lm[j] = rowMean[ good[j] ];
        lv[j] = rowVar [ good[j] ];
      }
      DataFrame df = DataFrame::create(
        Named("mean")     = lm,
        Named("variance") = lv
      );
      if (display_progress) Rcout << "[SSVL] fitting loess on " << k << " points\n";

      Function loess("loess"), predict("predict");
      List fit = loess(_["formula"] = "log10(variance) ~ log10(mean)",
                       _["data"]    = df,
                       _["span"]    = 0.3);
      NumericVector fitted = predict(fit, _["newdata"] = df);

      for (int j = 0; j < k; ++j) {
        int r = good[j];
        double ev = std::pow(10.0, fitted[j]);
        expectedVar[r] = ev;
        sd[r] = (ev > 0 ? std::sqrt(ev) : 1.0);
      }
    }
    if (display_progress) Rcout << "[SSVL] computed expectedVar & sd\n";

    // 5) final pass: standardize each row (one loop over non-zeros + O(1) for zeros)
    NumericVector result(nrow, 0.0);
    double vmax = std::sqrt((double)ncol);
    for (int col = 0; col < ncol; ++col) {
      for (int idx = p[col]; idx < p[col+1]; ++idx) {
        int r = i[idx];
        double z = (x[idx] - rowMean[r]) / sd[r];
        // clamp
        if      (z >  vmax) z =  vmax;
        else if (z < -vmax) z = -vmax;
        result[r] += z * z;
      }
    }
    // add zero contributions
    for (int r = 0; r < nrow; ++r) {
      int nz = nnz[r];
      int zr = ncol - nz;
      double z0 = (0.0 - rowMean[r]) / sd[r];
      // clamp z0 as well
      if      (z0 >  vmax) z0 =  vmax;
      else if (z0 < -vmax) z0 = -vmax;
      result[r] += zr * (z0 * z0);

      // normalize by (ncol - 1) or ncol, depending on your definition
      result[r] = (ncol > 1) ? ( result[r] / (ncol - 1) ) : 0.0;
    }

    Rcout << "[SSVL] exiting\n";
    return result;
}

