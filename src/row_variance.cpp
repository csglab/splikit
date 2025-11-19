#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rowVariance_cpp(SEXP mat) {

  // Dense numeric matrix path
  if (Rf_isMatrix(mat) && TYPEOF(mat) == REALSXP) {
    NumericMatrix M(mat);
    int nrow = M.nrow(), ncol = M.ncol();

    NumericVector out(nrow);
    for (int i = 0; i < nrow; ++i) {
      double sum = 0.0, sum2 = 0.0;
      for (int j = 0; j < ncol; ++j) {
        double v = M(i, j);
        sum  += v;
        sum2 += v * v;
      }
      double mean = sum / ncol;
      out[i] = sum2 / ncol - mean * mean;
    }

    return out;
  }

  // Dense integer matrix path (Issue #16 from deep analysis)
  if (Rf_isMatrix(mat) && TYPEOF(mat) == INTSXP) {
    IntegerMatrix M(mat);
    int nrow = M.nrow(), ncol = M.ncol();

    NumericVector out(nrow);
    for (int i = 0; i < nrow; ++i) {
      double sum = 0.0, sum2 = 0.0;
      for (int j = 0; j < ncol; ++j) {
        double v = static_cast<double>(M(i, j));
        sum  += v;
        sum2 += v * v;
      }
      double mean = sum / ncol;
      out[i] = sum2 / ncol - mean * mean;
    }

    return out;
  }

  // Sparse matrix path (dgCMatrix)
  if (Rf_isS4(mat) && Rf_inherits(mat, "dgCMatrix")) {
    S4 M(mat);
    IntegerVector dims = M.slot("Dim");
    IntegerVector i     = M.slot("i");
    IntegerVector p     = M.slot("p");
    NumericVector x     = M.slot("x");

    int nrow = dims[0], ncol = dims[1];

    NumericVector rowSum(nrow, 0.0), rowSum2(nrow, 0.0);

    // compressed-column iteration
    for (int col = 0; col < ncol; ++col) {
      for (int idx = p[col]; idx < p[col + 1]; ++idx) {
        int row = i[idx];
        double v = x[idx];
        rowSum[row]  += v;
        rowSum2[row] += v * v;
      }
    }

    NumericVector out(nrow);
    for (int row = 0; row < nrow; ++row) {
      double mean = rowSum[row] / ncol;
      out[row] = rowSum2[row] / ncol - mean * mean;
    }

    return out;
  }

  // Unsupported type
  stop("`mat` must be a numeric matrix or a dgCMatrix.");
  return NumericVector(0); // never reached
}
