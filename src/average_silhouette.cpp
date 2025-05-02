#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
double silhouette_avg(const arma::mat& X, const IntegerVector& cluster_assignments, int n_threads = 1) {
  int n = X.n_rows;
  IntegerVector clusters = clone(cluster_assignments);
  IntegerVector unique_clusters = sort_unique(clusters);
  NumericVector silhouette_values(n);

  #pragma omp parallel for num_threads(n_threads)
  for (int i = 0; i < n; ++i) {
    int this_cluster = clusters[i];

    double a = 0.0;
    int a_count = 0;
    double b = R_PosInf;

    for (int c = 0; c < unique_clusters.size(); ++c) {
      int other_cluster = unique_clusters[c];
      double total_dist = 0.0;
      int count = 0;

      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        if (clusters[j] == other_cluster) {
          double dist_ij = arma::norm(X.row(i) - X.row(j), 2);
          total_dist += dist_ij;
          count++;
        }
      }

      if (count > 0) {
        double avg_dist = total_dist / count;
        if (other_cluster == this_cluster) {
          a = avg_dist;
          a_count = count;
        } else if (avg_dist < b) {
          b = avg_dist;
        }
      }
    }

    double s = 0.0;
    if (a_count > 0 && std::max(a, b) > 0) {
      s = (b - a) / std::max(a, b);
    }
    silhouette_values[i] = s;
  }

  return mean(silhouette_values);
}

