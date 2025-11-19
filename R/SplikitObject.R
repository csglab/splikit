#' @title SplikitObject
#' @description
#' R6 class for splicing analysis in single-cell RNA-seq data.
#' Provides a modern object-oriented interface to splikit functionality
#' while maintaining backward compatibility with existing functions.
#'
#' @details
#' The SplikitObject encapsulates the core data structures for splicing analysis:
#' \itemize{
#'   \item \code{m1}: Inclusion matrix (sparse dgCMatrix)
#'   \item \code{m2}: Exclusion matrix (sparse dgCMatrix)
#'   \item \code{eventData}: Event metadata (data.table)
#'   \item \code{geneExpression}: Optional gene expression matrix
#' }
#'
#' @examples
#' \dontrun{
#' # Create from junction abundance data
#' junction_ab <- load_toy_SJ_object()
#' obj <- SplikitObject$new(junction_ab = junction_ab)
#'
#' # Compute M2 and find variable events
#' obj$makeM2()
#' hve <- obj$findVariableEvents(min_row_sum = 50)
#'
#' # Or create from existing matrices
#' obj <- SplikitObject$new(m1 = my_m1, m2 = my_m2, eventData = my_eventdata)
#'
#' # Chain operations
#' results <- obj$makeM2()$findVariableEvents()
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom Matrix nnzero rowSums
#' @importFrom data.table is.data.table data.table copy
#' @export
SplikitObject <- R6::R6Class("SplikitObject",

  public = list(

    #' @field m1 Inclusion matrix (dgCMatrix). Rows are events, columns are cells.
    m1 = NULL,

    #' @field m2 Exclusion matrix (dgCMatrix). Same dimensions as m1.
    m2 = NULL,

    #' @field eventData Event metadata (data.table). One row per event.
    eventData = NULL,

    #' @field geneExpression Optional gene expression matrix (dgCMatrix).
    geneExpression = NULL,

    #' @field metadata List containing summary statistics and analysis results.
    metadata = NULL,

    #' @description
    #' Create a new SplikitObject.
    #'
    #' @param junction_ab A junction abundance object from \code{make_junction_ab()}.
    #'   If provided, m1 and eventData are computed automatically.
    #' @param m1 An existing inclusion matrix (dgCMatrix).
    #' @param m2 An existing exclusion matrix (dgCMatrix).
    #' @param eventData A data.table with event metadata.
    #' @param min_counts Minimum count threshold for filtering events (default: 1).
    #' @param verbose Print progress messages (default: FALSE).
    #'
    #' @return A new SplikitObject instance.
    initialize = function(junction_ab = NULL,
                          m1 = NULL, m2 = NULL, eventData = NULL,
                          min_counts = 1, verbose = FALSE) {

      if (!is.null(junction_ab)) {
        # Build from raw junction abundance data
        if (verbose) cat("Building SplikitObject from junction abundance data...\n")

        result <- make_m1(
          junction_ab_object = junction_ab,
          min_counts = min_counts,
          verbose = verbose
        )

        self$m1 <- result$m1_inclusion_matrix
        self$eventData <- result$event_data
        self$metadata <- list(
          creation_method = "junction_ab",
          summary = result$summary
        )

        if (verbose) cat("SplikitObject created with", nrow(self$m1), "events and",
                         ncol(self$m1), "cells.\n")

      } else if (!is.null(m1)) {
        # Build from existing matrices
        private$validateInputs(m1, m2, eventData)

        self$m1 <- private$ensureSparse(m1)
        if (!is.null(m2)) {
          self$m2 <- private$ensureSparse(m2)
        }
        self$eventData <- eventData
        self$metadata <- list(creation_method = "matrices")

      } else {
        private$error("Must provide either 'junction_ab' or 'm1' matrix")
      }

      invisible(self)
    },

    #' @description
    #' Compute the M2 exclusion matrix from M1 and eventData.
    #'
    #' @param batch_size Number of groups per batch for memory management (default: 5000).
    #' @param memory_threshold Maximum rows before switching to batched processing.
    #' @param force_fast Force fast processing regardless of size (default: FALSE).
    #' @param multi_thread Use parallel processing for batched operations (default: FALSE).
    #' @param verbose Print progress messages (default: FALSE).
    #'
    #' @return Self (invisibly), for method chaining.
    makeM2 = function(batch_size = 5000, memory_threshold = 2e9,
                      force_fast = FALSE, multi_thread = FALSE, verbose = FALSE) {

      if (is.null(self$m1)) {
        private$error("M1 matrix not initialized")
      }
      if (is.null(self$eventData)) {
        private$error("eventData not initialized. Cannot compute M2 without event metadata.")
      }

      self$m2 <- make_m2(
        m1_inclusion_matrix = self$m1,
        eventdata = self$eventData,
        batch_size = batch_size,
        memory_threshold = memory_threshold,
        force_fast = force_fast,
        multi_thread = multi_thread,
        verbose = verbose
      )

      invisible(self)
    },

    #' @description
    #' Find highly variable splicing events.
    #'
    #' @param min_row_sum Minimum row sum threshold for filtering (default: 50).
    #' @param n_threads Number of threads for parallel computation (default: 1).
    #' @param verbose Print progress messages (default: FALSE).
    #'
    #' @return A data.table with event names and sum_deviance scores.
    findVariableEvents = function(min_row_sum = 50, n_threads = 1, verbose = FALSE) {

      private$ensureM2Computed()

      # Pre-compute row sums efficiently
      m1_sums <- Matrix::rowSums(self$m1)
      m2_sums <- Matrix::rowSums(self$m2)

      # Check for empty results
      valid_count <- sum(m1_sums > min_row_sum & m2_sums > min_row_sum)
      if (valid_count == 0) {
        private$error(
          "No events pass min_row_sum threshold of ", min_row_sum,
          ". Consider lowering threshold. ",
          "Current range: m1 [", min(m1_sums), "-", max(m1_sums), "], ",
          "m2 [", min(m2_sums), "-", max(m2_sums), "]"
        )
      }

      result <- find_variable_events(
        m1_matrix = self$m1,
        m2_matrix = self$m2,
        min_row_sum = min_row_sum,
        n_threads = n_threads,
        verbose = verbose
      )

      # Store in metadata
      self$metadata$variableEvents <- result
      self$metadata$variableEvents_params <- list(
        min_row_sum = min_row_sum,
        n_events_tested = valid_count
      )

      return(result)
    },

    #' @description
    #' Find highly variable genes from gene expression data.
    #'
    #' @param method Method for variable gene selection: "vst" or "sum_deviance" (default: "vst").
    #' @param n_threads Number of threads for parallel computation (default: 1).
    #' @param verbose Print progress messages (default: FALSE).
    #'
    #' @return A data.table with gene names and variability scores.
    findVariableGenes = function(method = "vst", n_threads = 1, verbose = FALSE) {

      if (is.null(self$geneExpression)) {
        private$error("Gene expression matrix not set. Use $setGeneExpression() first.")
      }

      result <- find_variable_genes(
        gene_expression_matrix = self$geneExpression,
        method = method,
        n_threads = n_threads,
        verbose = verbose
      )

      self$metadata$variableGenes <- result
      self$metadata$variableGenes_params <- list(method = method)

      return(result)
    },

    #' @description
    #' Compute pseudo-correlation between splicing and external data.
    #'
    #' @param ZDB_matrix Dense matrix of external data (e.g., gene expression PCs).
    #'   Must have same dimensions as m1.
    #' @param metric R-squared metric: "CoxSnell" or "Nagelkerke" (default: "CoxSnell").
    #' @param suppress_warnings Suppress computation warnings (default: TRUE).
    #'
    #' @return A data.table with event names, pseudo_correlation, and null_distribution.
    getPseudoCorrelation = function(ZDB_matrix, metric = "CoxSnell",
                                     suppress_warnings = TRUE) {

      private$ensureM2Computed()

      # Dimension validation
      if (!is.matrix(ZDB_matrix)) {
        private$error("ZDB_matrix must be a dense matrix")
      }
      if (nrow(ZDB_matrix) != nrow(self$m1)) {
        private$error(
          "Row mismatch: m1 has ", nrow(self$m1),
          " rows but ZDB_matrix has ", nrow(ZDB_matrix), " rows"
        )
      }
      if (ncol(ZDB_matrix) != ncol(self$m1)) {
        private$error(
          "Column mismatch: m1 has ", ncol(self$m1),
          " columns but ZDB_matrix has ", ncol(ZDB_matrix), " columns"
        )
      }

      result <- get_pseudo_correlation(
        ZDB_matrix = ZDB_matrix,
        m1_inclusion = self$m1,
        m2_exclusion = self$m2,
        metric = metric,
        suppress_warnings = suppress_warnings
      )

      # Warn about NA removal
      n_before <- nrow(self$m1)
      n_after <- nrow(result)
      if (n_before > n_after) {
        n_removed <- n_before - n_after
        message("Removed ", n_removed, " event(s) with NA values (",
                round(100 * n_removed / n_before, 1), "% of total).")
        message("Reasons: insufficient data (n<2), no variation, or convergence failure.")
      }

      return(result)
    },

    #' @description
    #' Subset the object by events and/or cells.
    #'
    #' @param events Event indices or names to keep.
    #' @param cells Cell indices or names to keep.
    #'
    #' @return Self (invisibly), for method chaining.
    subset = function(events = NULL, cells = NULL) {

      if (!is.null(events)) {
        # Convert names to indices if needed
        if (is.character(events)) {
          events <- match(events, rownames(self$m1))
          if (any(is.na(events))) {
            n_missing <- sum(is.na(events))
            warning(n_missing, " event name(s) not found and will be ignored")
            events <- events[!is.na(events)]
          }
        }

        if (length(events) == 0) {
          private$error("Subsetting would remove all events")
        }

        self$m1 <- self$m1[events, , drop = FALSE]
        if (!is.null(self$m2)) {
          self$m2 <- self$m2[events, , drop = FALSE]
        }
        if (!is.null(self$eventData)) {
          self$eventData <- self$eventData[events, ]
        }
      }

      if (!is.null(cells)) {
        # Convert names to indices if needed
        if (is.character(cells)) {
          cells <- match(cells, colnames(self$m1))
          if (any(is.na(cells))) {
            n_missing <- sum(is.na(cells))
            warning(n_missing, " cell name(s) not found and will be ignored")
            cells <- cells[!is.na(cells)]
          }
        }

        if (length(cells) == 0) {
          private$error("Subsetting would remove all cells")
        }

        self$m1 <- self$m1[, cells, drop = FALSE]
        if (!is.null(self$m2)) {
          self$m2 <- self$m2[, cells, drop = FALSE]
        }
        if (!is.null(self$geneExpression)) {
          self$geneExpression <- self$geneExpression[, cells, drop = FALSE]
        }
      }

      invisible(self)
    },

    #' @description
    #' Set the gene expression matrix.
    #'
    #' @param gene_matrix A gene expression matrix (will be converted to dgCMatrix).
    #'
    #' @return Self (invisibly), for method chaining.
    setGeneExpression = function(gene_matrix) {

      # Validate dimensions
      if (ncol(gene_matrix) != ncol(self$m1)) {
        private$error(
          "Gene expression must have same number of cells as m1 (",
          ncol(self$m1), "), got ", ncol(gene_matrix)
        )
      }

      # Ensure sparse
      self$geneExpression <- private$ensureSparse(gene_matrix)

      invisible(self)
    },

    #' @description
    #' Annotate events with gene information from a GTF file.
    #'
    #' @param GTF_file Path to a GTF annotation file.
    #'
    #' @return Self (invisibly), for method chaining.
    annotateEvents = function(GTF_file) {

      # File validation
      if (!file.exists(GTF_file)) {
        private$error("GTF file not found: ", GTF_file)
      }

      if (is.null(self$eventData)) {
        private$error("eventData not initialized. Cannot annotate events.")
      }

      self$eventData <- make_eventdata_plus(
        eventdata = data.table::copy(self$eventData),
        GTF_file_direction = GTF_file
      )

      invisible(self)
    },

    #' @description
    #' Get a summary of the object.
    #'
    #' @return A list with object statistics.
    summary = function() {
      m1_nnz <- Matrix::nnzero(self$m1)
      m1_total <- prod(dim(self$m1))

      list(
        events = nrow(self$m1),
        cells = ncol(self$m1),
        has_m2 = !is.null(self$m2),
        has_gene_expression = !is.null(self$geneExpression),
        n_genes = if (!is.null(self$geneExpression)) nrow(self$geneExpression) else 0,
        sparsity_m1 = round(1 - (m1_nnz / m1_total), 4),
        nnz_m1 = m1_nnz,
        eventData_columns = if (!is.null(self$eventData)) names(self$eventData) else NULL,
        metadata_keys = names(self$metadata)
      )
    },

    #' @description
    #' Print a human-readable summary of the object.
    print = function() {
      cat("SplikitObject\n")
      cat(paste(rep("-", 40), collapse = ""), "\n")
      cat("Events:        ", format(nrow(self$m1), big.mark = ","), "\n")
      cat("Cells:         ", format(ncol(self$m1), big.mark = ","), "\n")
      cat("M2 computed:   ", !is.null(self$m2), "\n")

      if (!is.null(self$geneExpression)) {
        cat("Gene expression: ", format(nrow(self$geneExpression), big.mark = ","), " genes\n")
      } else {
        cat("Gene expression: not set\n")
      }

      m1_nnz <- Matrix::nnzero(self$m1)
      m1_total <- prod(dim(self$m1))
      cat("Sparsity (M1): ", round((1 - m1_nnz / m1_total) * 100, 1), "%\n")
      cat("Memory (M1):   ", format(object.size(self$m1), units = "auto"), "\n")

      if (length(self$metadata) > 0) {
        cat("Metadata:      ", paste(names(self$metadata), collapse = ", "), "\n")
      }

      invisible(self)
    },

    #' @description
    #' Create a deep copy of the object.
    #'
    #' @return A new SplikitObject with copied data.
    clone = function(deep = TRUE) {
      if (deep) {
        new_obj <- SplikitObject$new(
          m1 = self$m1,
          m2 = self$m2,
          eventData = if (!is.null(self$eventData)) data.table::copy(self$eventData) else NULL
        )
        if (!is.null(self$geneExpression)) {
          new_obj$geneExpression <- self$geneExpression
        }
        new_obj$metadata <- self$metadata
        return(new_obj)
      } else {
        # Shallow copy (default R6 behavior)
        super$clone(deep = FALSE)
      }
    }
  ),

  private = list(

    #' Validate input matrices and eventData
    validateInputs = function(m1, m2, eventData) {

      # Check m1
      if (!inherits(m1, "Matrix") && !is.matrix(m1)) {
        private$error("m1 must be a matrix or sparse Matrix object")
      }

      # Check m2 if provided
      if (!is.null(m2)) {
        if (!inherits(m2, "Matrix") && !is.matrix(m2)) {
          private$error("m2 must be a matrix or sparse Matrix object")
        }
        if (!identical(dim(m1), dim(m2))) {
          private$error(
            "m1 and m2 must have identical dimensions. ",
            "m1: ", nrow(m1), "x", ncol(m1), ", ",
            "m2: ", nrow(m2), "x", ncol(m2)
          )
        }
        if (!is.null(rownames(m1)) && !is.null(rownames(m2))) {
          if (!identical(rownames(m1), rownames(m2))) {
            private$error("m1 and m2 must have identical row names")
          }
        }
        if (!is.null(colnames(m1)) && !is.null(colnames(m2))) {
          if (!identical(colnames(m1), colnames(m2))) {
            private$error("m1 and m2 must have identical column names")
          }
        }
      }

      # Check eventData if provided
      if (!is.null(eventData)) {
        if (!data.table::is.data.table(eventData) && !is.data.frame(eventData)) {
          private$error("eventData must be a data.table or data.frame")
        }
        if (nrow(eventData) != nrow(m1)) {
          private$error(
            "eventData must have same number of rows as m1. ",
            "eventData: ", nrow(eventData), " rows, ",
            "m1: ", nrow(m1), " rows"
          )
        }
      }
    },

    #' Ensure M2 is computed
    ensureM2Computed = function() {
      if (is.null(self$m2)) {
        private$error("M2 not computed. Call $makeM2() first.")
      }
    },

    #' Ensure matrix is sparse dgCMatrix
    ensureSparse = function(mat) {
      if (inherits(mat, "dgCMatrix")) {
        return(mat)
      } else if (inherits(mat, "Matrix")) {
        return(methods::as(mat, "dgCMatrix"))
      } else if (is.matrix(mat)) {
        return(methods::as(mat, "dgCMatrix"))
      } else {
        private$error("Cannot convert to sparse matrix: unknown type")
      }
    },

    #' Standardized error reporting
    error = function(...) {
      msg <- paste0(...)
      stop(msg, call. = FALSE)
    }
  )
)


#' Create a SplikitObject
#'
#' Convenience function to create a SplikitObject.
#'
#' @param ... Arguments passed to \code{SplikitObject$new()}.
#'
#' @return A new SplikitObject instance.
#'
#' @examples
#' \dontrun{
#' # From junction abundance
#' junction_ab <- load_toy_SJ_object()
#' obj <- splikit(junction_ab = junction_ab)
#'
#' # From existing matrices
#' obj <- splikit(m1 = my_m1, m2 = my_m2, eventData = my_eventdata)
#' }
#'
#' @export
splikit <- function(...) {
  SplikitObject$new(...)
}
