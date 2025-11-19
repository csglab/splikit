// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Memory-efficient M2 computation building CSC format directly
//'
//' This version minimizes memory usage by:
//' - Building CSC format directly (no triplet intermediate)
//' - Using O(n_groups) workspace per column instead of dense group_sums matrix
//' - Two-pass algorithm: count then fill
//'
//' Memory usage: O(nnz_output) + O(n_groups) workspace
//' vs previous: O(n_groups * n_cells) + O(6 * nnz_output)
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

#ifdef _OPENMP
  if (n_threads > 1) {
    omp_set_num_threads(n_threads);
  }
#endif

  // Pre-compute group members
  // group_members[g] = vector of event indices in group g
  std::vector<std::vector<int>> group_members(n_groups);
  for (int i = 0; i < n_events; i++) {
    group_members[group_ids[i]].push_back(i);
  }

  // ============================================================================
  // PASS 1: Count non-zeros per column
  // ============================================================================

  std::vector<size_t> col_nnz(n_cells, 0);

#ifdef _OPENMP
#pragma omp parallel if(n_threads > 1)
{
  // Thread-local workspace for group sums
  std::vector<double> group_sum_local(n_groups);

#pragma omp for schedule(dynamic)
  for (int j = 0; j < n_cells; j++) {
    // Reset workspace
    std::fill(group_sum_local.begin(), group_sum_local.end(), 0.0);

    // Compute group sums for this column
    for (sp_mat::const_col_iterator it = M1.begin_col(j); it != M1.end_col(j); ++it) {
      int row = it.row();
      group_sum_local[group_ids[row]] += *it;
    }

    // Count non-zeros for this column
    size_t count = 0;
    for (int g = 0; g < n_groups; g++) {
      double total = group_sum_local[g];
      if (total == 0) continue;

      const std::vector<int>& members = group_members[g];
      for (int i : members) {
        double m1_val = M1(i, j);
        double m2_val = total - m1_val;
        if (m2_val != 0) {
          count++;
        }
      }
    }
    col_nnz[j] = count;
  }
}
#else
  // Single-threaded version
  std::vector<double> group_sum(n_groups);
  for (int j = 0; j < n_cells; j++) {
    std::fill(group_sum.begin(), group_sum.end(), 0.0);

    for (sp_mat::const_col_iterator it = M1.begin_col(j); it != M1.end_col(j); ++it) {
      int row = it.row();
      group_sum[group_ids[row]] += *it;
    }

    size_t count = 0;
    for (int g = 0; g < n_groups; g++) {
      double total = group_sum[g];
      if (total == 0) continue;

      const std::vector<int>& members = group_members[g];
      for (int i : members) {
        double m1_val = M1(i, j);
        double m2_val = total - m1_val;
        if (m2_val != 0) {
          count++;
        }
      }
    }
    col_nnz[j] = count;
  }
#endif

  // ============================================================================
  // Build column pointers (CSC format)
  // ============================================================================

  std::vector<size_t> col_ptr(n_cells + 1);
  col_ptr[0] = 0;
  for (int j = 0; j < n_cells; j++) {
    col_ptr[j + 1] = col_ptr[j] + col_nnz[j];
  }
  size_t total_nnz = col_ptr[n_cells];

  // Allocate output arrays (exactly sized - no waste)
  std::vector<uword> row_indices(total_nnz);
  std::vector<double> values(total_nnz);

  // ============================================================================
  // PASS 2: Fill in values
  // ============================================================================

#ifdef _OPENMP
#pragma omp parallel if(n_threads > 1)
{
  // Thread-local workspace
  std::vector<double> group_sum_local(n_groups);

#pragma omp for schedule(dynamic)
  for (int j = 0; j < n_cells; j++) {
    // Reset workspace
    std::fill(group_sum_local.begin(), group_sum_local.end(), 0.0);

    // Compute group sums for this column
    for (sp_mat::const_col_iterator it = M1.begin_col(j); it != M1.end_col(j); ++it) {
      int row = it.row();
      group_sum_local[group_ids[row]] += *it;
    }

    // Fill in M2 values for this column
    size_t write_idx = col_ptr[j];

    // We need to write rows in sorted order for valid CSC format
    // Collect (row, value) pairs first, then sort
    std::vector<std::pair<uword, double>> col_entries;
    col_entries.reserve(col_nnz[j]);

    for (int g = 0; g < n_groups; g++) {
      double total = group_sum_local[g];
      if (total == 0) continue;

      const std::vector<int>& members = group_members[g];
      for (int i : members) {
        double m1_val = M1(i, j);
        double m2_val = total - m1_val;
        if (m2_val != 0) {
          col_entries.emplace_back(static_cast<uword>(i), m2_val);
        }
      }
    }

    // Sort by row index (required for CSC format)
    std::sort(col_entries.begin(), col_entries.end(),
              [](const std::pair<uword, double>& a, const std::pair<uword, double>& b) {
                return a.first < b.first;
              });

    // Write to output arrays
    for (const auto& entry : col_entries) {
      row_indices[write_idx] = entry.first;
      values[write_idx] = entry.second;
      write_idx++;
    }
  }
}
#else
  // Single-threaded version
  std::vector<double> group_sum2(n_groups);
  for (int j = 0; j < n_cells; j++) {
    std::fill(group_sum2.begin(), group_sum2.end(), 0.0);

    for (sp_mat::const_col_iterator it = M1.begin_col(j); it != M1.end_col(j); ++it) {
      int row = it.row();
      group_sum2[group_ids[row]] += *it;
    }

    size_t write_idx = col_ptr[j];
    std::vector<std::pair<uword, double>> col_entries;
    col_entries.reserve(col_nnz[j]);

    for (int g = 0; g < n_groups; g++) {
      double total = group_sum2[g];
      if (total == 0) continue;

      const std::vector<int>& members = group_members[g];
      for (int i : members) {
        double m1_val = M1(i, j);
        double m2_val = total - m1_val;
        if (m2_val != 0) {
          col_entries.emplace_back(static_cast<uword>(i), m2_val);
        }
      }
    }

    std::sort(col_entries.begin(), col_entries.end(),
              [](const std::pair<uword, double>& a, const std::pair<uword, double>& b) {
                return a.first < b.first;
              });

    for (const auto& entry : col_entries) {
      row_indices[write_idx] = entry.first;
      values[write_idx] = entry.second;
      write_idx++;
    }
  }
#endif

  // ============================================================================
  // Construct sparse matrix directly from CSC components
  // ============================================================================

  // Convert to arma types
  arma::uvec rowind(total_nnz);
  arma::vec vals(total_nnz);
  arma::uvec colptr(n_cells + 1);

  for (size_t k = 0; k < total_nnz; k++) {
    rowind[k] = row_indices[k];
    vals[k] = values[k];
  }
  for (int j = 0; j <= n_cells; j++) {
    colptr[j] = col_ptr[j];
  }

  // Create sparse matrix from CSC components
  arma::sp_mat M2(rowind, colptr, vals, n_events, n_cells);

  return M2;
}


//' Legacy function - redirects to memory-efficient version
//'
//' @param M1 Sparse matrix (dgCMatrix) of inclusion counts (events x cells)
//' @param group_ids Integer vector of group IDs for each event
//' @param n_threads Number of threads for OpenMP (default 1)
//'
//' @return Sparse matrix M2 with same dimensions as M1
//'
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat make_m2_cpp(const arma::sp_mat& M1,
                          const IntegerVector& group_ids,
                          int n_threads = 1) {
  // Redirect to memory-efficient parallel version
  return make_m2_cpp_parallel(M1, group_ids, n_threads);
}
