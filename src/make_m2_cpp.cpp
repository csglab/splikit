// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Compute M2 exclusion matrix from M1 inclusion matrix (C++ implementation)
//'
//' For each event in a group, M2 = group_sum - M1.
//' This is a fast C++ implementation with optional OpenMP parallelization.
//'
//' @param M1 Sparse matrix (dgCMatrix) of inclusion counts (events x cells)
//' @param group_ids Integer vector of group IDs for each event (0-indexed internally)
//' @param n_threads Number of threads for OpenMP (default 1)
//'
//' @return Sparse matrix M2 with same dimensions as M1
//'
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat make_m2_cpp(const arma::sp_mat& M1,
                          const IntegerVector& group_ids,
                          int n_threads = 1) {

  int n_events = M1.n_rows;
  int n_cells = M1.n_cols;

  // Validate inputs

  if (group_ids.size() != n_events) {
    stop("group_ids length must match number of rows in M1");
  }

  // Find number of unique groups

  int max_group = 0;
  for (int i = 0; i < n_events; i++) {
    if (group_ids[i] > max_group) {
      max_group = group_ids[i];
    }
  }
  int n_groups = max_group + 1;

  Rcout << "Computing M2 matrix (C++ implementation)\n";
  Rcout << "  Events: " << n_events << ", Cells: " << n_cells << ", Groups: " << n_groups << "\n";

#ifdef _OPENMP
  if (n_threads > 1) {
    omp_set_num_threads(n_threads);
    Rcout << "  Using " << n_threads << " OpenMP threads\n";
  }
#else
  if (n_threads > 1) {
    Rcout << "  OpenMP not available, running single-threaded\n";
  }
#endif

  // Step 1: Pre-compute which events belong to each group
  // group_members[g] = vector of event indices in group g
  std::vector<std::vector<int>> group_members(n_groups);
  for (int i = 0; i < n_events; i++) {
    int g = group_ids[i];
    group_members[g].push_back(i);
  }

  // Step 2: Compute group sums for each cell
  // group_sums[g * n_cells + j] = sum of M1 for all events in group g, cell j
  // Using dense matrix for group sums since we need random access
  arma::mat group_sums(n_groups, n_cells, arma::fill::zeros);

  // Iterate through sparse matrix to compute group sums
  // Sparse matrices are stored in CSC format (column-major)
  for (sp_mat::const_iterator it = M1.begin(); it != M1.end(); ++it) {
    int row = it.row();
    int col = it.col();
    double val = *it;
    int g = group_ids[row];
    group_sums(g, col) += val;
  }

  Rcout << "  Group sums computed\n";

  // Step 3: Build M2 matrix
  // For efficiency, we'll build triplets and construct sparse matrix
  // M2[i, j] = group_sums[group_ids[i], j] - M1[i, j]

  // Count non-zeros: M2 is non-zero where either:
  // - M1 is non-zero, or
  // - group_sum is non-zero and M1 is zero

  // We'll use a different approach: iterate by groups
  // For each group, for each cell, distribute the exclusion counts

  std::vector<uword> row_indices;
  std::vector<uword> col_indices;
  std::vector<double> values;

  // Reserve space (estimate based on M1 density * group expansion)
  size_t estimated_nnz = M1.n_nonzero * 2;
  row_indices.reserve(estimated_nnz);
  col_indices.reserve(estimated_nnz);
  values.reserve(estimated_nnz);

  // Process each group
  // Note: We process sequentially here to avoid thread-safety issues with vectors
  // But the inner loops could be parallelized

  for (int g = 0; g < n_groups; g++) {
    const std::vector<int>& members = group_members[g];
    int group_size = members.size();

    if (group_size <= 1) {
      // Single-member groups have M2 = 0 everywhere
      continue;
    }

    // For each cell
    for (int j = 0; j < n_cells; j++) {
      double total = group_sums(g, j);

      if (total == 0) {
        // No counts in this group for this cell
        continue;
      }

      // For each event in this group
      for (int k = 0; k < group_size; k++) {
        int i = members[k];
        double m1_val = M1(i, j);
        double m2_val = total - m1_val;

        if (m2_val != 0) {
          row_indices.push_back(i);
          col_indices.push_back(j);
          values.push_back(m2_val);
        }
      }
    }
  }

  Rcout << "  Triplets computed: " << values.size() << " non-zeros\n";

  // Construct sparse matrix from triplets
  arma::umat locations(2, values.size());
  arma::vec vals(values.size());

  for (size_t k = 0; k < values.size(); k++) {
    locations(0, k) = row_indices[k];
    locations(1, k) = col_indices[k];
    vals(k) = values[k];
  }

  arma::sp_mat M2(locations, vals, n_events, n_cells);

  Rcout << "  M2 matrix constructed\n";

  return M2;
}


//' Faster M2 computation using parallel column processing
//'
//' This version processes columns in parallel, which is more efficient
//' for the CSC sparse matrix format.
//'
//' @param M1 Sparse matrix (dgCMatrix) of inclusion counts (events x cells)
//' @param group_ids Integer vector of group IDs for each event
//' @param n_threads Number of threads for OpenMP (default 1)
//'
//' @return Sparse matrix M2 with same dimensions as M1
//'
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat make_m2_cpp_parallel(const arma::sp_mat& M1,
                                   const IntegerVector& group_ids,
                                   int n_threads = 1) {

  int n_events = M1.n_rows;
  int n_cells = M1.n_cols;

  // Validate inputs
  if (group_ids.size() != n_events) {
    stop("group_ids length must match number of rows in M1");
  }

  // Find number of unique groups
  int max_group = 0;
  for (int i = 0; i < n_events; i++) {
    if (group_ids[i] > max_group) {
      max_group = group_ids[i];
    }
  }
  int n_groups = max_group + 1;

  Rcout << "Computing M2 matrix (C++ parallel implementation)\n";
  Rcout << "  Events: " << n_events << ", Cells: " << n_cells << ", Groups: " << n_groups << "\n";

#ifdef _OPENMP
  if (n_threads > 1) {
    omp_set_num_threads(n_threads);
    Rcout << "  Using " << n_threads << " OpenMP threads\n";
  }
#else
  if (n_threads > 1) {
    Rcout << "  OpenMP not available, running single-threaded\n";
  }
#endif

  // Pre-compute group members
  std::vector<std::vector<int>> group_members(n_groups);
  for (int i = 0; i < n_events; i++) {
    group_members[group_ids[i]].push_back(i);
  }

  // Process each column in parallel
  // Each thread builds its own triplet list, then we merge

  std::vector<std::vector<uword>> all_rows(n_cells);
  std::vector<std::vector<uword>> all_cols(n_cells);
  std::vector<std::vector<double>> all_vals(n_cells);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(n_threads > 1)
#endif
  for (int j = 0; j < n_cells; j++) {
    // Compute group sums for this column
    std::vector<double> group_sum(n_groups, 0.0);

    // Get column slice from sparse matrix
    for (sp_mat::const_col_iterator it = M1.begin_col(j); it != M1.end_col(j); ++it) {
      int row = it.row();
      double val = *it;
      group_sum[group_ids[row]] += val;
    }

    // Compute M2 values for this column
    for (int g = 0; g < n_groups; g++) {
      double total = group_sum[g];
      if (total == 0) continue;

      const std::vector<int>& members = group_members[g];
      for (int i : members) {
        double m1_val = M1(i, j);
        double m2_val = total - m1_val;

        if (m2_val != 0) {
          all_rows[j].push_back(i);
          all_cols[j].push_back(j);
          all_vals[j].push_back(m2_val);
        }
      }
    }
  }

  // Count total non-zeros
  size_t total_nnz = 0;
  for (int j = 0; j < n_cells; j++) {
    total_nnz += all_vals[j].size();
  }

  Rcout << "  Triplets computed: " << total_nnz << " non-zeros\n";

  // Merge all triplets
  arma::umat locations(2, total_nnz);
  arma::vec vals(total_nnz);

  size_t idx = 0;
  for (int j = 0; j < n_cells; j++) {
    for (size_t k = 0; k < all_vals[j].size(); k++) {
      locations(0, idx) = all_rows[j][k];
      locations(1, idx) = all_cols[j][k];
      vals(idx) = all_vals[j][k];
      idx++;
    }
  }

  arma::sp_mat M2(locations, vals, n_events, n_cells);

  Rcout << "  M2 matrix constructed\n";

  return M2;
}
