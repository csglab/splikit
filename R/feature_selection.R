#' Calculate the Sum Deviance for Inclusion and Exclusion Matrices
#'
#' @param m1_matrix A matrix representing the inclusion matrix. Rows are events, columns are barcodes.
#' @param m2_matrix A matrix representing the exclusion matrix. Rows are events, columns are barcodes.
#' @param n_threads If the module OpenPM is available for your device, the function suggests using multi-thread processing for even faster computation.
#' @param min_row_sum A numeric value specifying the minimum row sum threshold for filtering events. Defaults to 50.
#' @param verbose Logical. If \code{TRUE}, prints progress and informational messages. Default is \code{FALSE}.
#' @param ... Additional arguments to be passed.
#'
#' @return A \code{data.table} containing the events and their corresponding sum deviance values.
#'
#' @examples
#' # loading the toy dataset
#' toy_obj <- load_toy_M1_M2_object()
#'
#' # getting HVE (high variable events)
#'  HVE <- find_variable_events(toy_obj$m1, toy_obj$m2)
#'
#'  # printing the results
#'  print(HVE[order(-sum_deviance)])
#' @export
find_variable_events <- function(m1_matrix, m2_matrix, min_row_sum = 50, n_threads = 1, verbose = FALSE, ...) {

  # Check if matrices are sparse
  if (!(inherits(m1_matrix, "Matrix") && inherits(m2_matrix, "Matrix"))) {
    stop("Both 'm1_matrix' and 'm2_matrix' must be sparse matrices of class 'Matrix'.")
  }
  # Check matrix compatibility
  if (!identical(colnames(m1_matrix), colnames(m2_matrix))) {
    stop("The colnames (barcodes) of inclusion and exclusion matrices are not identical.")
  }
  if (!identical(rownames(m1_matrix), rownames(m2_matrix))) {
    stop("The rownames (junction events) of inclusion and exclusion matrices are not identical.")
  }

  # Filter rows based on minimum row sum criteria
  to_keep_events <- which(rowSums(m1_matrix) > min_row_sum & rowSums(m2_matrix) > min_row_sum)
  m1_matrix <- m1_matrix[to_keep_events, , drop = FALSE]
  m2_matrix <- m2_matrix[to_keep_events, , drop = FALSE]

  # Create metadata table
  temp_current_barcodes <- data.table::data.table(brc = colnames(m1_matrix))
  temp_current_barcodes$ID <- sub("^.{16}-(.*$)", "\\1", temp_current_barcodes$brc)
  meta <- temp_current_barcodes

  libraries <- unique(meta$ID)
  if (verbose) message("There are ", length(libraries), " libraries detected...")

  # Initialize deviance sum vector
  sum_deviances <- numeric(nrow(m1_matrix))
  names(sum_deviances) <- rownames(m1_matrix)

  deviance_results <- lapply(libraries, function(lib) {
    filter <- which(meta[, ID] == lib)
    M1_sub <- m1_matrix[, filter, drop = FALSE]
    M2_sub <- m2_matrix[, filter, drop = FALSE]

    deviance_values <- tryCatch({
      calcDeviances_ratio(M1_sub, M2_sub, n_threads)
    }, error = function(e) {
      stop("Error in calcDeviances_ratio function: ", e$message)
    })
    deviance_values <- c(deviance_values)
    names(deviance_values) <- rownames(M1_sub)
    if (verbose) {
      message("Calculating the deviances for sample ", lib, " has been completed!")
    }
    return(deviance_values)
  })

  # Combine all results
  sum_deviances <- Reduce(`+`, deviance_results)
  rez <- data.table::data.table(events = names(sum_deviances), sum_deviance = as.numeric(sum_deviances))
  return(rez)
}

#' Find Variable Genes Using Variance or Deviance-Based Metrics
#'
#' @description
#' Identifies highly variable genes from a sparse gene expression matrix using one of two methods:
#' variance-stabilizing transformation (VST) or deviance-based modeling. The VST method uses a C++-accelerated
#' approach to compute standardized variance, while the deviance-based method models gene variability
#' across libraries using negative binomial deviances.
#'
#' @param gene_expression_matrix A sparse gene expression matrix (of class \code{Matrix}) with gene names as row names.
#' @param method Character string, either \code{"vst"} or \code{"sum_deviance"}. The default is \code{"sum_deviance"}.
#'   \code{"vst"} uses a variance-stabilizing transformation to identify variable genes.
#'   \code{"sum_deviance"} computes per-library deviances and combines them with a row variance metric.
#' @param n_threads If OpenMP is available for your device, the function suggests using multi-thread processing for even faster computation (only for sum_deviance method).
#' @param verbose Logical. If \code{TRUE}, prints progress and informational messages. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{data.table} containing gene names (column \code{events}) and computed metrics.
#'   For the deviance method, this includes \code{sum_deviance} and \code{variance} columns.
#'
#' @examples
#' library(Matrix)
#' # loading the toy dataset
#' toy_obj <- load_toy_M1_M2_object()
#'
#' # getting high variable genes
#' HVG_VST <- find_variable_genes(toy_obj$gene_expression, method = "vst") # vst method
#' # sum_deviance method
#' HVG_DEV <- find_variable_genes(toy_obj$gene_expression, method = "sum_deviance")
#' 
#' # Using multi-threading for faster computation (sum_deviance method only)
#' HVG_DEV_MT <- find_variable_genes(toy_obj$gene_expression, 
#'                                   method = "sum_deviance", 
#'                                   n_threads = 4) # 4 threads
#'
#' # printing the results
#' print(HVG_VST[order(-standardize_variance)])
#' print(HVG_DEV[order(-sum_deviance)])
#'
#' @useDynLib splikit, .registration = TRUE
#' @import Rcpp
#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix dsCMatrix dgTMatrix dsTMatrix
#' @export
find_variable_genes <- function(gene_expression_matrix, method = "vst", n_threads = 1, verbose = FALSE, ...) {
  # adding the vst method as the default
  method <- match.arg(method, choices = c("vst", "sum_deviance"))

  # Verify that gene_expression_matrix is a sparse Matrix
  if (!inherits(gene_expression_matrix, "Matrix")) {
    stop("The 'gene_expression_matrix' must be a sparse matrix of class 'Matrix'.")
  }

  if (method == "vst") {
    if (verbose) message("The method we are using is vst (Seurat)...")
    if (!exists("standardizeSparse_variance_vst")) {
      stop("The function 'standardizeSparse_variance_vst' is not available. Check your C++ source files.")
    }
    rez_vector <- tryCatch({
      standardizeSparse_variance_vst(matSEXP = gene_expression_matrix)
    }, error = function(e) {
      stop("Error in standardizeSparse_variance_vst: ", e$message)
    })
    rez <- data.table::data.table(events = rownames(gene_expression_matrix),
                                  standardize_variance = rez_vector)
  } else {
    if (verbose) message("The method we are using is like deviance summarion per library...")

    # Filter rows based on minimum row sum criteria
    to_keep_features <- which(rowSums(gene_expression_matrix) > 0)
    if (length(to_keep_features) == 0) {
      stop("No genes with a positive row sum were found.")
    }
    gene_expression_matrix <- gene_expression_matrix[to_keep_features, , drop = FALSE]

    # Create metadata table using column names
    temp_current_barcodes <- data.table::data.table(brc = colnames(gene_expression_matrix))
    temp_current_barcodes$ID <- sub("^.{16}-(.*$)", "\\1", temp_current_barcodes$brc)
    meta <- temp_current_barcodes

    libraries <- unique(meta$ID)
    if (verbose) message("There are ", length(libraries), " libraries detected...")

    # Initialize deviance sum vector with gene names
    sum_deviances <- numeric(nrow(gene_expression_matrix))
    names(sum_deviances) <- rownames(gene_expression_matrix)

    calculate_all_deviances <- function(gene_expression_matrix, meta, ID, n_threads, verbose) {
      # Get unique libraries
      libraries <- unique(meta[, ID])

      # Initialize results list
      deviances_list <- vector("list", length(libraries))
      names(deviances_list) <- libraries

      # Getting the sum deviances for NB model
      deviances_list <- lapply(libraries, function(lib) {
        filter <- which(meta[, ID] == lib)
        gene_expression_matrix_sub <- gene_expression_matrix[, filter, drop = FALSE]
        deviance_values <- tryCatch({
          calcNBDeviancesWithThetaEstimation(as(gene_expression_matrix_sub, "dgCMatrix"), n_threads)
        }, error = function(e) {
          stop("Error in calcNBDeviancesWithThetaEstimation function: ", e$message)
        })
        deviance_values <- c(deviance_values)
        names(deviance_values) <- rownames(gene_expression_matrix_sub)
        if (verbose) {
          message("Calculating the deviances for sample ", lib, " has been completed!")
        }
        return(deviance_values)
      })

      # Sum deviances across all libraries
      sum_deviances <- Reduce(`+`, deviances_list)

      return(sum_deviances)
    }

    # Then use it:
    sum_deviances <- calculate_all_deviances(gene_expression_matrix, meta, ID, n_threads, verbose)

    # Compute row variance using the previously defined function
    row_var <- tryCatch({
      splikit::get_rowVar(M = gene_expression_matrix)
    }, error = function(e) {
      stop("Error in splikit::get_rowVar: ", e$message)
    })

    row_var_cpp_dt <- data.table::data.table(events = rownames(gene_expression_matrix),
                                             variance = row_var)

    rez <- data.table::data.table(events = names(sum_deviances),
                                  sum_deviance = as.numeric(sum_deviances))
    rez <- base::merge(rez, row_var_cpp_dt, by = "events")
    data.table::setkey(x = rez, NULL)
  }

  return(rez)
}
