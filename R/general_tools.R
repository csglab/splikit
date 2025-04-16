#' Compute Pseudo-Correlation Using Beta-Binomial Model
#'
#' This function calculates a pseudo R²-like correlation metric using a beta-binomial model
#' implemented in C++. It takes in a data matrix `ZDB_matrix` and two model matrices
#' for inclusion and exclusion, respectively.
#'
#' @param ZDB_matrix A numeric dense matrix of shape (events x samples). Should have rownames representing events.
#' @param m1_inclusion A numeric dense matrix of the same number of rows as `ZDB_matrix`, representing inclusion features.
#' @param m2_exclusion A numeric dense matrix of the same number of rows as `ZDB_matrix`, representing exclusion features.
#'
#' @return A `data.table` with the following columns:
#' \describe{
#'   \item{event}{The event names from `ZDB_matrix` rownames.}
#'   \item{pseudo_correlation}{The computed pseudo R² correlation values.}
#'   \item{null_distribution}{Null correlation values from a permuted version of `ZDB_matrix`.}
#' }
#'
#' @examples
#' set.seed(123)
#'
#' # Create a dummy ZDB matrix (5 events x 4 samples)
#' ZDB_matrix <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' rownames(ZDB_matrix) <- paste0("Event", 1:5)
#'
#' # Create dummy inclusion and exclusion matrices
#' m1_inclusion <- matrix(rnorm(15), nrow = 5, ncol = 3)
#' m2_exclusion <- matrix(rnorm(10), nrow = 5, ncol = 2)
#'
#' # Run the function
#' result <- multigedi_get_pseudo_correlation(ZDB_matrix, m1_inclusion, m2_exclusion)
#' print(result)
#'
#' @export
multigedi_get_pseudo_correlation <- function(ZDB_matrix, m1_inclusion, m2_exclusion) {
  
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

  suppressPackageStartupMessages(Rcpp::sourceCpp("./src/cpp_pseudoR2.cpp"))
  
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
  
  results <- na.omit(results)
  
  cat("Computation completed successfully.\n")
  return(results)
}

