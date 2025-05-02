#' @useDynLib splikit, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Compute Pseudo-Correlation Using Beta-Binomial Model
#'
#' This function calculates a pseudo R²-like correlation metric using a beta-binomial model
#' implemented in C++. It takes in a data matrix `ZDB_matrix` and two model matrices
#' for inclusion and exclusion, respectively. Optionally, it can suppress warnings
#' generated during computation.
#'
#' @param ZDB_matrix A numeric dense matrix of shape (events x samples). Should have rownames representing events.
#' @param m1_inclusion A numeric dense matrix of the same number of rows as `ZDB_matrix`, representing inclusion features.
#' @param m2_exclusion A numeric dense matrix of the same number of rows as `ZDB_matrix`, representing exclusion features.
#' @param suppress_warnings Logical. If \code{TRUE} (default), suppresses warnings during any warnings triggered during
#' computation (e.g., due to ill-conditioned inputs)
#'
#' @return A `data.table` with the following columns:
#' \describe{
#'   \item{event}{The event names from `ZDB_matrix` rownames.}
#'   \item{pseudo_correlation}{The computed pseudo R² correlation values.}
#'   \item{null_distribution}{Null correlation values from a permuted version of `ZDB_matrix`.}
#' }
#'
#' @examples
#' set.seed(42)
#' # get the m1 object
#' junction_abundance_object <- load_toy_SJ_object()
#' m1_obj <- make_m1(junction_ab_object = junction_abundance_object)
#'
#' # obtaining the m1 and eventdata
#' m1_inclusion <- m1_obj$m1_inclusion_matrix
#' eventdata <- m1_obj$event_data
#' m2_exclusion <- make_m2(m1_inclusion, eventdata)
#'
#' # creating a dummy ZDB
#' ZDB_matrix <- matrix(rnorm(n = (nrow(m1_inclusion) * ncol(m1_inclusion)), sd = 7),
#' nrow = nrow(m1_inclusion),
#' ncol = ncol(m1_inclusion))
#' rownames(ZDB_matrix) <- rownames(m1_inclusion)
#'
#' # both m1 and m2 need to be a dense matrix
#' m1_inclusion <- as.matrix(m1_inclusion)
#' m2_exclusion <- as.matrix(m2_exclusion)
#'
#' # getting the pseudo-correlation
#' pseudo_r_squre <- get_pseudo_correlation(ZDB_matrix, m1_inclusion, m2_exclusion)
#' print(pseudo_r_squre)
#'
#' @export
get_pseudo_correlation <- function(ZDB_matrix, m1_inclusion, m2_exclusion, suppress_warnings=TRUE) {

  if (!is.matrix(ZDB_matrix)) stop("ZDB_matrix must be a matrix.")
  if (!is.matrix(m1_inclusion)) stop("m1_inclusion must be a matrix.")
  if (!is.matrix(m2_exclusion)) stop("m2_exclusion must be a matrix.")

  if (nrow(m1_inclusion) != nrow(ZDB_matrix)) {
    stop("m1_inclusion must have the same number of rows as ZDB_matrix.")
  }
  if (nrow(m2_exclusion) != nrow(ZDB_matrix)) {
    stop("m2_exclusion must have the same number of rows as ZDB_matrix.")
  }

  cat("Input dimension check passed. Proceeding with computation.\n")
  if(suppress_warnings){suppressWarnings({

      correls <- cppBetabinPseudoR2(Z = ZDB_matrix,
                                    m1 = as.matrix(m1_inclusion),
                                    m2 = as.matrix(m2_exclusion))

      if (is.null(rownames(ZDB_matrix))) {
        warning("ZDB_matrix has no row names; assigning default event names.")
        events <- paste0("Event_", seq_len(nrow(ZDB_matrix)))
      } else {
        events <- rownames(ZDB_matrix)
      }

      null_dist <- cppBetabinPseudoR2(Z = ZDB_matrix[, sample.int(ncol(ZDB_matrix))],
                                      m1 = as.matrix(m1_inclusion),
                                      m2 = as.matrix(m2_exclusion))

      results <- data.table::data.table(
        event = events,
        pseudo_correlation = correls,
        null_distribution = null_dist

      )
    })
  }

  results <- na.omit(results)

  cat("Computation completed successfully.\n")
  return(results)
}

#' Calculate Row-wise Variance for Dense or Sparse Matrices
#'
#' @description
#' Efficiently computes the variance of each row for either a base R dense matrix
#' or a sparse dgCMatrix, via a single Rcpp entry point. Logs progress messages
#' to the R console.
#' @param M A numeric matrix (base R matrix) or a sparse matrix of class \code{"dgCMatrix"}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress and informational messages.
#' @return A numeric vector of length \code{nrow(M)} containing the variance of each row.
#' @details
#' Dispatches in C++ between dense and sparse implementations to avoid unnecessary
#' overhead or external dependencies. Uses compressed-column traversal for sparse inputs.
#' @examples
#' \dontrun{
#'   library(Matrix)
#'   # Dense example
#'   dm <- matrix(rnorm(1000), nrow = 100)
#'   get_rowVar(dm)
#'   # Sparse example
#'   sm <- rsparsematrix(100, 10, density = 0.1)
#'   get_rowVar(sm)
#' }
#' @note
#' Only 32-bit integer indices are supported, due to limitations in R's internal matrix representations.
#' This function will not work with matrices that exceed the 32-bit integer indexing range.
#' @export
get_rowVar <- function(M, verbose=FALSE) {
  if (! (is.matrix(M) || inherits(M, "dgCMatrix")) ) {
    stop("`M` must be either a dense numeric matrix or a dgCMatrix.")
  }
  if(verbose) message("[get_rowVar] Starting computation")
  result <- rowVariance_cpp(M)
  if(verbose) message("[get_rowVar] Computation finished")
  result
}

#' Compute Average Silhouette Width with Logging
#'
#' Computes the average silhouette width for a clustering solution using Euclidean distance.
#'
#' @param X A numeric matrix where rows are observations and columns are features.
#' @param cluster_assignments An integer vector of cluster assignments, which must be the same length as the number of rows in \code{X}.
#' @param n_threads Number of threads to use for parallel processing.
#'
#' @note This process can be very slow for large matrices if single-threaded. Use multiple threads to take advantage of parallel computation for significantly faster results.
#'
#' @return A single numeric value: the average silhouette score.
#'
#' @examples
#' # Preparing the inputs
#' set.seed(42)
#' pc_matrix <- matrix(data = rnorm(n = 10000 * 15, sd = 2), nrow = 10000, ncol = 15)
#' cluster_numbers <- runif(n = 10000, min = 1, max = 10) |> as.integer()
#'
#' # Getting the mean silhouette score
#' n_threads <- parallel::detectCores()
#' score <- get_silhouette_mean(pc_matrix, cluster_numbers, n_threads)
#' print(score)
#'
#' @export
get_silhouette_mean <- function(X, cluster_assignments, n_threads = 1) {
  stopifnot(is.matrix(X), is.numeric(X))
  stopifnot(is.integer(cluster_assignments) || is.numeric(cluster_assignments))
  stopifnot(nrow(X) == length(cluster_assignments))

  cat("[silhouette_avg] Starting computation...\n")
  cat(sprintf("[silhouette_avg] Using %d threads...", n_threads), "\n")

  score <- silhouette_avg(X, as.integer(cluster_assignments), as.integer(n_threads))
  return(score)
}
