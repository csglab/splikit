// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// A small epsilon used to clamp probabilities away from 0 and 1.
const double EPS = 1e-8;

// Logistic (sigmoid) function.
inline double logistic(double eta) {
  return 1.0 / (1.0 + std::exp(-eta));
}

// Helper: Compute deviance for one observation in binomial regression.
// n = m1 + m2 and mu = n * p (fitted successes).
inline double deviance_obs(double m1, double n, double mu) {
  double dev = 0.0;
  if (m1 > 0)
    dev += 2.0 * m1 * std::log(m1 / (mu + EPS));
  if ((n - m1) > 0)
    dev += 2.0 * (n - m1) * std::log((n - m1) / ((n - mu) + EPS));
  return dev;
}

// Function: Fit a binomial GLM with one covariate (plus intercept) using IRLS.
// Parameters:
//   x: predictor vector
//   m1: successes vector
//   m2: failures vector
// Returns a List containing:
//   beta: the coefficient vector [intercept, slope]
//   deviance: the model deviance (sum over the observations)
//   n: number of observations used.
List fit_logistic_regression(const vec& x, const vec& m1, const vec& m2) {
  int nobs = x.n_elem;
  // Design matrix: column 0 is intercept, column 1 is predictor x.
  mat X(nobs, 2);
  X.col(0).ones();
  X.col(1) = x;
  
  // Response is successes (m1) with number of trials = m1 + m2.
  vec y = m1;
  vec n_trials = m1 + m2;
  
  // Initialize beta to zeros.
  vec beta = zeros<vec>(2);
  
  int max_iter = 100;
  double tol = 1e-6;
  
  for (int iter = 0; iter < max_iter; iter++) {
    // Compute the linear predictor: eta = X * beta.
    vec eta = X * beta;
    // Compute fitted probabilities with clamping.
    vec p = eta;
    for (uword i = 0; i < p.n_elem; i++) {
      p(i) = std::min(1.0 - EPS, std::max(EPS, logistic(eta(i))));
    }
    // Fitted mean successes: mu = n_trials * p.
    vec mu = n_trials % p;
    
    // Weight vector: w = n_trials * p * (1-p)
    vec w = n_trials % (p % (1 - p));
    
    // Compute the working variable z.
    vec z = eta + (y - mu) / (w + EPS);
    
    // IRLS update: beta_new = inv(X' W X) * (X' W z)
    mat W = diagmat(w);
    mat XtWX = X.t() * W * X;
    vec XtWz = X.t() * W * z;
    
    vec beta_new;
    // Use no_approx option to disable approximation for ill-conditioned matrices
    // Warnings are handled at R level via suppressWarnings()
    bool solved = solve(beta_new, XtWX, XtWz, solve_opts::no_approx);
    if (!solved) {
      return List::create(Named("beta") = NumericVector::create(NA_REAL, NA_REAL),
                          Named("deviance") = NA_REAL,
                          Named("nobs") = nobs);
    }
    
    if (norm(beta_new - beta, "inf") < tol) {
      beta = beta_new;
      break;
    }
    
    beta = beta_new;
  }
  
  // After convergence, compute deviance.
  vec eta = X * beta;
  vec p = eta;
  for (uword i = 0; i < p.n_elem; i++) {
    p(i) = std::min(1.0 - EPS, std::max(EPS, logistic(eta(i))));
  }
  vec mu = n_trials % p;
  
  double dev = 0.0;
  for (int i = 0; i < nobs; i++) {
    dev += deviance_obs(y(i), n_trials(i), mu(i));
  }
  
  return List::create(Named("beta") = beta,
                      Named("deviance") = dev,
                      Named("nobs") = nobs);
}

// Function: Fit an intercept-only binomial GLM.
// Returns a List containing beta (intercept), deviance, and number of observations.
List fit_logistic_regression_null(const vec& m1, const vec& m2) {
  int nobs = m1.n_elem;
  vec n_trials = m1 + m2;
  
  double sum_y = sum(m1);
  double sum_n = sum(n_trials);
  
  // Empirical probability clamped by EPS
  double p_hat = sum_y / sum_n;
  p_hat = std::min(1.0 - EPS, std::max(EPS, p_hat));
  
  vec mu = n_trials * p_hat;
  
  double dev = 0.0;
  for (int i = 0; i < nobs; i++) {
    dev += deviance_obs(m1(i), n_trials(i), mu(i));
  }
  
  double beta0 = std::log(p_hat / (1.0 - p_hat));
  
  return List::create(Named("beta") = NumericVector::create(beta0),
                      Named("deviance") = dev,
                      Named("nobs") = nobs);
}

// Template function to handle both dense and sparse matrices
template<typename T1, typename T2>
NumericVector cppBetabinPseudoR2_impl(const arma::mat& Z,
                                      const T1& m1,
                                      const T2& m2,
                                      const std::string& metric) {
  int nrows = Z.n_rows;
  int p = Z.n_cols; // number of columns / observations per row
  NumericVector result(nrows, NA_REAL); // to store pseudo R² for each row
  
  // Validate metric parameter
  bool use_nagelkerke = (metric == "nagelkerke" || metric == "Nagelkerke");
  bool use_coxsnell = (metric == "coxsnell" || metric == "CoxSnell" || metric == "cox-snell" || metric == "Cox-Snell");
  
  if (!use_nagelkerke && !use_coxsnell) {
    Rcpp::stop("Invalid metric. Must be 'CoxSnell' or 'Nagelkerke'");
  }
  
  for (int i = 0; i < nrows; i++) {
    // Extract predictor vector for row i.
    rowvec x_full = Z.row(i);
    
    // Extract row i from m1 and m2 element-wise.
    vec m1_row(p), m2_row(p);
    for (int j = 0; j < p; j++) {
      m1_row(j) = m1(i, j);
      m2_row(j) = m2(i, j);
    }
    
    // Determine valid indices where total trials (m1 + m2) are greater than 0.
    std::vector<int> idx;
    for (int j = 0; j < p; j++) {
      if ((m1_row(j) + m2_row(j)) > 0) {
        idx.push_back(j);
      }
    }
    
    int n_valid = idx.size();
    if (n_valid < 2) {
      result[i] = NA_REAL;
      continue;
    }
    
    // Create filtered vectors.
    vec x(n_valid);
    vec m1_sub(n_valid);
    vec m2_sub(n_valid);
    for (int j = 0; j < n_valid; j++) {
      x(j)      = x_full(idx[j]);
      m1_sub(j) = m1_row(idx[j]);
      m2_sub(j) = m2_row(idx[j]);
    }
    
    // If all successes or all failures, we cannot fit the model reliably.
    if (sum(m1_sub) == 0 || sum(m2_sub) == 0) {
      result[i] = NA_REAL;
      continue;
    }
    
    // Fit the model. Check if the binomial (logistic) model is adequate.
    double max_m1 = m1_sub.max();
    double max_m2 = m2_sub.max();
    
    List fit_full, fit_reduced;
    if ((max_m1 <= 1.0) && (max_m2 <= 1.0)) {
      // Full model: predictor + intercept.
      fit_full = fit_logistic_regression(x, m1_sub, m2_sub);
      // Reduced model: intercept only.
      fit_reduced = fit_logistic_regression_null(m1_sub, m2_sub);
    } else {
      // For now, use the same logistic regression fit as an approximation.
      fit_full = fit_logistic_regression(x, m1_sub, m2_sub);
      fit_reduced = fit_logistic_regression_null(m1_sub, m2_sub);
    }
    
    // Check if the regression converged.
    vec beta_full = fit_full["beta"];
    if (beta_full.has_nan() || beta_full.n_elem < 2) {
      result[i] = NA_REAL;
      continue;
    }
    
    double dev_full = as<double>(fit_full["deviance"]);
    double dev_reduced = as<double>(fit_reduced["deviance"]);
    int n_obs = as<int>(fit_full["nobs"]);
    
    // Compute Cox-Snell R².
    double CoxSnell_R2 = 1.0 - std::exp((dev_full - dev_reduced) / double(n_obs));
    
    // Return the appropriate metric
    if (use_nagelkerke) {
      // Compute Nagelkerke R²
      double denom = 1.0 - std::exp(-dev_reduced / double(n_obs));
      double Nagelkerke_R2 = (denom != 0.0) ? CoxSnell_R2 / denom : NA_REAL;
      
      // Return the signed square-root of Nagelkerke_R2, using the sign of the predictor coefficient.
      double slope = beta_full(1);
      double Nagelkerke_r = (Nagelkerke_R2 >= 0.0) ? std::sqrt(Nagelkerke_R2) * ((slope >= 0) ? 1.0 : -1.0) : NA_REAL;
      result[i] = Nagelkerke_r;
    } else {
      // Return the signed square-root of CoxSnell_R2 (default), using the sign of the predictor coefficient.
      double slope = beta_full(1);
      double CoxSnell_r = (CoxSnell_R2 >= 0.0) ? std::sqrt(CoxSnell_R2) * ((slope >= 0) ? 1.0 : -1.0) : NA_REAL;
      result[i] = CoxSnell_r;
    }
  }
  
  return result;
}

// Main exported function for dense matrices (backward compatibility)
// [[Rcpp::export]]
NumericVector cppBetabinPseudoR2(const arma::mat& Z,
                                 const arma::mat& m1,
                                 const arma::mat& m2,
                                 std::string metric = "CoxSnell") {
  return cppBetabinPseudoR2_impl(Z, m1, m2, metric);
}

// Exported function for sparse m1 and m2 matrices
// [[Rcpp::export]]
NumericVector cppBetabinPseudoR2_sparse(const arma::mat& Z,
                                        const arma::sp_mat& m1,
                                        const arma::sp_mat& m2,
                                        std::string metric = "CoxSnell") {
  return cppBetabinPseudoR2_impl(Z, m1, m2, metric);
}

// Exported function for mixed types: sparse m1, dense m2
// [[Rcpp::export]]
NumericVector cppBetabinPseudoR2_mixed1(const arma::mat& Z,
                                        const arma::sp_mat& m1,
                                        const arma::mat& m2,
                                        std::string metric = "CoxSnell") {
  return cppBetabinPseudoR2_impl(Z, m1, m2, metric);
}

// Exported function for mixed types: dense m1, sparse m2
// [[Rcpp::export]]
NumericVector cppBetabinPseudoR2_mixed2(const arma::mat& Z,
                                        const arma::mat& m1,
                                        const arma::sp_mat& m2,
                                        std::string metric = "CoxSnell") {
  return cppBetabinPseudoR2_impl(Z, m1, m2, metric);
}