#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#######################################################################
###################### junction abundance #############################
#######################################################################


#' multigedi_make_junction_ab
#'
#' This function processes STARsolo splicing junction output to create a modality list for splicing data.
#' It supports both single and multiple sample processing, optional barcode filtration, and internal whitelist usage.
#'
#' @param STARsolo_SJ_dirs A character vector or list of strings representing the paths to STARsolo SJ directories.
#'                         Each directory should contain the raw splicing junction output files.
#' @param white_barcode_lists A list of character vectors, each containing barcode whitelist(s) for the corresponding sample.
#'                            If `NULL` (default), the function uses the internal STARsolo whitelist if `use_internal_whitelist` is `TRUE`.
#' @param sample_ids A character vector or list of unique sample IDs corresponding to each directory in `STARsolo_SJ_dirs`.
#' @param use_internal_whitelist A logical flag (default `TRUE`) indicating whether to use the internal STARsolo whitelist
#'                                located at `../Gene/filtered/barcodes.tsv` for each sample when `white_barcode_lists` is `NULL`.
#' @return A list containing processed splicing modality data for each sample. If a single sample is provided,
#'         the function returns the processed data as a list. For multiple samples, a named list is returned.
#'
#'         Each sample's result contains:
#'           - `eventdata`: A `data.table` containing feature metadata.
#'           - `junction_ab`: A sparse matrix of junction abundance.
#'
#' @examples
#' # Single sample processing with external whitelist
#' STARsolo_SJ_dir <- "path/to/sample1"
#' white_barcode_list <- c("barcode1", "barcode2")
#' sample_id <- "sample1"
#' result <- multigedi_make_junction_ab(STARsolo_SJ_dir, white_barcode_list, sample_id)
#'
#' # Multiple samples processing with internal whitelist
#' STARsolo_SJ_dirs <- list("path/to/sample1", "path/to/sample2")
#' sample_ids <- c("sample1", "sample2")
#' result <- multigedi_make_junction_ab(STARsolo_SJ_dirs, NULL, sample_ids)
#'
#' @export
multigedi_make_junction_ab <- function(STARsolo_SJ_dirs, white_barcode_lists = NULL, sample_ids, use_internal_whitelist = TRUE) {
  
  # Handle single sample input by converting to lists
  if (!is.list(STARsolo_SJ_dirs)) STARsolo_SJ_dirs <- list(STARsolo_SJ_dirs)
  if (!is.list(sample_ids)) sample_ids <- list(sample_ids)
  
  # If white_barcode_lists is NULL, create a list of NULL values
  if (is.null(white_barcode_lists)) {
    white_barcode_lists <- vector("list", length(STARsolo_SJ_dirs))
  }
  
  # Ensure inputs have the same length
  if (!all(length(STARsolo_SJ_dirs) == length(white_barcode_lists), length(white_barcode_lists) == length(sample_ids))) {
    stop("All input lists (STARsolo_SJ_dirs, white_barcode_lists, sample_ids) must have the same length.", call. = FALSE)
  }
  
  # Helper function to process one sample
  process_sample <- function(STARsolo_SJ_dir, white_barcode_list, sample_id) {
    # Define paths
    mtx_dir <- paste0(STARsolo_SJ_dir, "/raw/matrix.mtx")
    feature_dir <- paste0(STARsolo_SJ_dir, "/raw/features.tsv")
    barcodes_dir <- paste0(STARsolo_SJ_dir, "/raw/barcodes.tsv")
    internal_whitelist_dir <- paste0(STARsolo_SJ_dir, "/../Gene/filtered/barcodes.tsv")
    
    # Check for required files
    if (!file.exists(mtx_dir)) stop("No abundance matrix in STARsolo SJ direction for sample: ", sample_id, call. = FALSE)
    if (!file.exists(feature_dir)) stop("No feature matrix in STARsolo SJ direction for sample: ", sample_id, call. = FALSE)
    if (!file.exists(barcodes_dir)) stop("No barcode matrix in STARsolo SJ direction for sample: ", sample_id, call. = FALSE)
    
    # Read splicing data
    cat("├── Processing sample: ", sample_id, "\n")
    mtx <- Matrix::readMM(mtx_dir)
    raw_brc <- data.table::fread(barcodes_dir, header = FALSE, showProgress = FALSE)
    feature <- data.table::fread(feature_dir, select = c(1, 2, 3, 4, 5, 6), 
                                 col.names = c('chr', 'start', 'end', 'strand', "intron_motif", 'is_annot'),
                                 showProgress = FALSE)
    
    # Use internal whitelist if enabled and no external whitelist is provided
    if (use_internal_whitelist && is.null(white_barcode_list)) {
      if (file.exists(internal_whitelist_dir)) {
        white_barcode_list <- data.table::fread(internal_whitelist_dir, header = FALSE, showProgress = FALSE)$V1
        cat("│  ├──  Using STARsolo internal whitelist for sample: ", sample_id, "\n")
      } else {
        stop("Internal whitelist not found at ", internal_whitelist_dir, " for sample: ", sample_id, call. = FALSE)
      }
    }
    
    # Add assets to eventdata
    feature[, index := .I]
    feature[, start_cor_id := paste0(chr, "_", start, "|", strand)]
    feature[, end_cor_id := paste0(chr, "_", end, "|", strand)]
    feature[, row_names_mtx := paste0(chr, ":", start, "-", end)]
    feature[, sample_id := sample_id]
    
    # Assign rownames and colnames to mtx
    colnames(mtx) <- raw_brc$V1
    rownames(mtx) <- feature$row_names_mtx
    
    # Check if rownames of mtx match row_names_mtx
    if (!identical(rownames(mtx), feature$row_names_mtx)) {
      stop("Mismatch between rownames of mtx and row_names_mtx in eventdata for sample: ", sample_id, call. = FALSE)
    }
    
    # Trim matrix if barcode list is provided
    if (!is.null(white_barcode_list)) {
      initial_barcodes <- ncol(mtx)
      mtx <- mtx[, colnames(mtx) %in% white_barcode_list, drop = FALSE]
      final_barcodes <- ncol(mtx)
      
      # Validate trimming
      if (final_barcodes == initial_barcodes) {
        warning("Barcode trimming did not reduce the number of barcodes for sample: ", sample_id, call. = FALSE)
      }
      if (final_barcodes == 0) {
        stop("All barcodes were removed after trimming for sample: ", sample_id, call. = FALSE)
      }
      
      cat("│  ├──  Trimmed junction abundance matrix for sample: ", sample_id, " (", final_barcodes, " barcodes remaining)\n")
    } else {
      cat("│  ├──  No barcode filtration applied for sample: ", sample_id, "\n")
    }
    
    # Append sample_id to column names
    colnames(mtx) <- paste0(colnames(mtx), "-", sample_id)
    
    # Save splicing modality
    m1 <- list()
    m1[["eventdata"]] <- feature
    m1[["junction_ab"]] <- as(mtx, "CsparseMatrix")
    
    # Collect summary information
    summary_row <- data.frame(
      Sample = sample_id,
      Barcodes = ncol(mtx),
      Events = nrow(mtx),
      stringsAsFactors = FALSE
    )
    
    cat("│  └── Finished processing sample: ", sample_id, "\n")
    return(list(result = m1, summary = summary_row))
  }
  
  # Use mapply to process all samples
  results <- mapply(
    process_sample,
    STARsolo_SJ_dirs,
    white_barcode_lists,
    sample_ids,
    SIMPLIFY = FALSE
  )
  
  # Extract results and summary
  final_results <- lapply(results, function(x) x$result)
  summary_table <- do.call(rbind, lapply(results, function(x) x$summary))
  
  # Print summary table
  cat("\nSummary of Processed Samples in M1 matrix:\n")
  print(summary_table)
  
  # Return all results
  if (length(final_results) == 1) {
    return(final_results[[1]])  # Return the single result as a list
  } else {
    names(final_results) <- sample_ids
    return(final_results)
  }
}


m1 <- multigedi_make_junction_ab(STARsolo_SJ_dirs  = paste0(list.files("/home/arsham79/scratch/Sandrine-s-project/results/star-output-new-whitelist/", full.names = TRUE), "/Solo.out/SJ/") |> as.list(), 
                        sample_ids = sub(pattern = "(^.*)_MPS.*$", 
                                                 replacement = "\\1", 
                                                 x = list.files("/home/arsham79/scratch/Sandrine-s-project/results/star-output-new-whitelist/")) |> as.list()
)



#######################################################################
######################### make M1 #####################################
#######################################################################

#' multigedi_make_m1
#'
#' This function processes junction abundance data from multiple samples to create a splicing modality inclusion matrix (M1).
#' It merges event data, handles start and end coordinate groups, ensures matrix compatibility, and includes robust error handling.
#'
#' The function requires the following libraries: `dplyr`, `data.table`, and `Matrix`.
#'
#' @param junction_ab_object A named list where each element represents a sample's junction abundance data. 
#'                           Each element must contain `eventdata` and a sparse matrix.
#' @return A list containing:
#'         - `m1_inclusion_matrix`: The processed inclusion matrix for all events across all samples.
#'         - `event_data`: A data table containing the merged and grouped event metadata.
#'
#' @examples
#' # Example usage
#' result <- multigedi_make_m1(junction_ab_object)
#' m1_matrix <- result$m1_inclusion_matrix
#' event_metadata <- result$event_data
#'
#' @export
multigedi_make_m1 <- function(junction_ab_object) {
  # Validate inputs
  if (!is.list(junction_ab_object) || length(junction_ab_object) == 0) {
    stop("`junction_ab_object` must be a non-empty list.")
  }
  
  # Extract and combine all eventdata
  all_in_one_eventdata <- tryCatch({
    do.call(rbind, lapply(junction_ab_object, \(x) x$eventdata))
  }, error = function(e) {
    stop("Error combining `eventdata`: ", e$message)
  })
  
  if (!"row_names_mtx" %in% colnames(all_in_one_eventdata)) {
    stop("`eventdata` must contain a column named `row_names_mtx`.")
  }
  
  # Remove `sample_id` and deduplicate events
  all_junctions <- unique(all_in_one_eventdata[, row_names_mtx])
  all_in_one_eventdata[, sample_id := NULL]
  temp_eventdata <- all_in_one_eventdata[match(all_junctions, row_names_mtx), ]
  
  # Grouping based on start and end coordinates
  temp_eventdata_grouped <- temp_eventdata %>%
    group_by(start_cor_id) %>%
    mutate(start_cor_group_id = cur_group_id(), start_cor_group_count = n()) %>%
    ungroup() %>%
    group_by(end_cor_id) %>%
    mutate(end_cor_group_id = cur_group_id(), end_cor_group_count = n()) %>%
    data.table::as.data.table()
  
  # Process start coordinate groups
  temp_eventdata_grouped_start <- data.table::copy(temp_eventdata_grouped)
  temp_eventdata_grouped_start <- temp_eventdata_grouped_start[start_cor_group_count != 1, ]
  temp_eventdata_grouped_start[, start_cor_group_id := paste0(start_cor_group_id, "_S")]
  temp_eventdata_grouped_start[, row_names_mtx := paste0(row_names_mtx, "_S")]
  setnames(temp_eventdata_grouped_start, old = c("row_names_mtx", "row_names_mtx_new"), 
           new = c("raw_row_names_mtx", "row_names_mtx"))
  
  # Process end coordinate groups
  temp_eventdata_grouped_end <- data.table::copy(temp_eventdata_grouped)
  temp_eventdata_grouped_end <- temp_eventdata_grouped_end[end_cor_group_count != 1, ]
  temp_eventdata_grouped_end[, end_cor_group_id := paste0(end_cor_group_id, "_E")]
  temp_eventdata_grouped_end[, row_names_mtx := paste0(row_names_mtx, "_E")]
  setnames(temp_eventdata_grouped_end, old = c("row_names_mtx", "row_names_mtx_new"), 
           new = c("raw_row_names_mtx", "row_names_mtx"))
  
  # Combine start and end groups
  colnames(temp_eventdata_grouped_end) <- colnames(temp_eventdata_grouped_start)
  eventdata <- rbind(temp_eventdata_grouped_start, temp_eventdata_grouped_end)
  
  # Rename columns for clarity
  eventdata <- eventdata[, .(
    chr, start, end, strand, intron_motif, is_annot, 
    start_cor_id, end_cor_id, raw_row_names_mtx, group_id, 
    group_count, row_names_mtx
  )]
  
  # Create a table for all events
  all_events <- eventdata[, .(raw_row_names_mtx, row_names_mtx)]
  all_events[, group_type := sub("^.*(.$)", "\\1", row_names_mtx)]
  
  # Process and merge all junction abundance matrices
  process_junction_ab_matrices <- function(junction_ab_object, all_events, eventdata) {
    # Initialize the final merged matrix
    m1_raw <- Matrix::sparseMatrix(
      i = integer(0), 
      j = integer(0), 
      dims = c(nrow(eventdata), 0),
      dimnames = list(eventdata$row_names_mtx, NULL)
    )
    
    lapply(seq_along(junction_ab_object), function(j) {
      # Extract current matrix and events
      little_m1 <- junction_ab_object[[j]]
      current_matrix <- little_m1[[2]]
      current_events <- data.table(raw_row_names_mtx = rownames(current_matrix))
      
      # Process start and end events
      start_events <- all_events[group_type == "S"]
      matched_start <- merge(current_events, start_events, by = "raw_row_names_mtx", sort = FALSE)
      start_matrix <- current_matrix[matched_start$raw_row_names_mtx, ]
      rownames(start_matrix) <- matched_start$row_names_mtx
      
      end_events <- all_events[group_type == "E"]
      matched_end <- merge(current_events, end_events, by = "raw_row_names_mtx", sort = FALSE)
      end_matrix <- current_matrix[matched_end$raw_row_names_mtx, ]
      rownames(end_matrix) <- matched_end$row_names_mtx
      
      # Merge start and end matrices
      merged_matrix <- rbind(start_matrix, end_matrix)
      
      # Add missing events
      missing_events <- setdiff(eventdata$row_names_mtx, rownames(merged_matrix))
      if (length(missing_events) > 0) {
        missing_matrix <- Matrix(0, nrow = length(missing_events), ncol = ncol(current_matrix), sparse = TRUE)
        rownames(missing_matrix) <- missing_events
        merged_matrix <- rbind(merged_matrix, missing_matrix)
      }
      
      # Reorder rows
      final_matrix <- merged_matrix[eventdata$row_names_mtx, ]
      
      # Combine into the final matrix
      m1_raw <<- cbind(m1_raw, final_matrix)
      
      # Log progress
      cat("Processed matrix:", names(junction_ab_object)[j], "\n")
    })
    
    return(m1_raw)
  }
  
  # Process all matrices
  m1 <- process_junction_ab_matrices(junction_ab_object, all_events, eventdata)
  
  # Verify matrix row order and return results
  if (!identical(rownames(m1), eventdata$row_names_mtx)) {
    m1 <- m1[eventdata$row_names_mtx, ]
  }
  
  cat("Finished processing M1.\n")
  return(list(
    m1_inclusion_matrix = m1,
    event_data = eventdata
  ))
}




#######################################################################
##################### Gene expression #################################
#######################################################################

#' process_gene_expression
#'
#' This function processes gene expression data from a given directory and creates a sparse matrix
#' for gene expression. It supports barcode filtration using a provided whitelist or the filtered barcodes file.
#'
#' @param expression_dirs A character vector or list of strings representing paths to directories containing
#'                        the gene expression matrix (`matrix.mtx`), barcodes, and features files.
#' @param sample_ids A character vector or list of unique sample IDs corresponding to each directory in `expression_dirs`.
#' @param whitelist_barcodes A list of character vectors, each containing barcode whitelist(s) for the corresponding sample.
#'                            If `NULL` (default), the function uses the filtered barcodes file if available.
#' @param use_filtered A logical flag (default `TRUE`) indicating whether to use the `filtered` data for barcode filtration.
#'                     If `FALSE`, no barcode filtration is applied unless a whitelist is provided.
#' @return A list containing processed gene expression data for each sample. If a single sample is provided,
#'         the function returns the processed data as a sparse matrix. For multiple samples, a named list is returned.
#'
#' @examples
#' # Single sample processing with external whitelist
#' expression_dir <- "path/to/sample1"
#' sample_id <- "sample1"
#' whitelist_barcode <- c("barcode1", "barcode2")
#' result <- process_gene_expression(expression_dir, sample_id, list(whitelist_barcode), use_filtered = FALSE)
#'
#' # Multiple samples processing with default filtered data
#' expression_dirs <- list("path/to/sample1", "path/to/sample2")
#' sample_ids <- c("sample1", "sample2")
#' result <- process_gene_expression(expression_dirs, sample_ids, NULL, use_filtered = TRUE)
#'
#' @export
process_gene_expression <- function(expression_dirs, sample_ids, whitelist_barcodes = NULL, use_filtered = TRUE) {
  
  # Handle single sample input by converting to lists
  if (!is.list(expression_dirs)) expression_dirs <- list(expression_dirs)
  if (!is.list(sample_ids)) sample_ids <- list(sample_ids)
  
  # If whitelist_barcodes is NULL, create a list of NULL values
  if (is.null(whitelist_barcodes)) {
    whitelist_barcodes <- vector("list", length(expression_dirs))
  }
  
  # Ensure inputs have the same length
  if (!all(length(expression_dirs) == length(whitelist_barcodes), length(whitelist_barcodes) == length(sample_ids))) {
    stop("All input lists (expression_dirs, whitelist_barcodes, sample_ids) must have the same length.", call. = FALSE)
  }
  
  # Helper function to process one sample
  process_sample <- function(expression_dir, sample_id, whitelist_barcode) {
    # Determine directory (filtered or raw)
    data_dir <- if (use_filtered) paste0(expression_dir, "/filtered") else paste0(expression_dir, "/raw")
    
    # Define paths
    expression_matrix_dir <- paste0(data_dir, "/matrix.mtx")
    expression_barcodes_dir <- paste0(data_dir, "/barcodes.tsv")
    expression_features_dir <- paste0(data_dir, "/features.tsv")
    filtered_barcodes_dir <- paste0(expression_dir, "/filtered/barcodes.tsv")  # For barcode filtration
    
    # Check for required files
    if (!file.exists(expression_matrix_dir)) stop("No expression matrix file found for sample: ", sample_id, call. = FALSE)
    if (!file.exists(expression_barcodes_dir)) stop("No barcodes file found for sample: ", sample_id, call. = FALSE)
    if (!file.exists(expression_features_dir)) stop("No features file found for sample: ", sample_id, call. = FALSE)
    
    # Read gene expression data
    cat("├── Processing gene expression data for sample: ", sample_id, "\n")
    g_mtx <- Matrix::readMM(expression_matrix_dir)
    g_brc <- data.table::fread(expression_barcodes_dir, header = FALSE)
    g_feature <- data.table::fread(expression_features_dir, header = FALSE)
    
    # Set row and column names
    rownames(g_mtx) <- g_feature$V1
    colnames(g_mtx) <- g_brc$V1
    
    # Apply barcode filtration
    if (!is.null(whitelist_barcode)) {
      cat("│  ├──  Applying provided whitelist for sample: ", sample_id, "\n")
    } else if (use_filtered && file.exists(filtered_barcodes_dir)) {
      whitelist_barcode <- data.table::fread(filtered_barcodes_dir, header = FALSE)$V1
      cat("│  ├──  Using filtered barcodes for sample: ", sample_id, "\n")
    }
    
    if (!is.null(whitelist_barcode)) {
      initial_barcodes <- ncol(g_mtx)
      
      # Filter matrix
      g_mtx <- g_mtx[, colnames(g_mtx) %in% whitelist_barcode, drop = FALSE]
      final_barcodes <- ncol(g_mtx)
      
      # Validate barcode filtration
      if (final_barcodes == 0) {
        stop("All barcodes were removed after filtration for sample: ", sample_id, call. = FALSE)
      }
      
      cat("│  ├──  Filtered barcodes applied for sample: ", sample_id, 
          " (Remaining barcodes: ", final_barcodes, ")\n")
    } else {
      cat("│  ├──  No barcode filtration applied for sample: ", sample_id, "\n")
    }
    
    # Collect summary information
    summary_row <- data.frame(
      Sample = sample_id,
      Barcodes = ncol(g_mtx),
      Genes = nrow(g_mtx),
      stringsAsFactors = FALSE
    )
    
    # Return processed matrix and summary
    g_mtx <- as(g_mtx, "CsparseMatrix")
    cat("│  └── Finished processing gene expression data for sample: ", sample_id, "\n")
    return(list(gene_expression = g_mtx, summary = summary_row))
  }
  
  # Use mapply to process all samples
  results <- mapply(
    process_sample,
    expression_dirs,
    sample_ids,
    whitelist_barcodes,
    SIMPLIFY = FALSE
  )
  
  # Extract results and summaries
  final_results <- lapply(results, function(x) x$gene_expression)
  summary_table <- do.call(rbind, lapply(results, function(x) x$summary))
  
  # Print summary table
  cat("\nSummary of Processed Samples:\n")
  print(summary_table)
  
  # Return results
  if (length(final_results) == 1) {
    return(final_results[[1]])  # Return the single result as a sparse matrix
  } else {
    names(final_results) <- sample_ids
    return(final_results)
  }
}


gene <- process_gene_expression(expression_dirs = paste0(list.files("/home/arsham79/scratch/Sandrine-s-project/results/star-output-new-whitelist/", full.names = TRUE), "/Solo.out/Gene/") |> as.list(), 
                          sample_ids = sub(pattern = "(^.*)_MPS.*$", 
                                           replacement = "\\1", 
                                           x = list.files("/home/arsham79/scratch/Sandrine-s-project/results/star-output-new-whitelist/")) |> as.list()
)




#######################################################################
################################ Velocyto #############################
#######################################################################

#' process_velocyto_data
#'
#' This function processes spliced and unspliced counts data from Velocyto outputs. It reads the respective matrices, 
#' sets appropriate row and column names, applies barcode filtration using a provided whitelist or the filtered barcodes file, 
#' and optionally merges counts across all samples into single large spliced and unspliced matrices.
#'
#' @param velocyto_dirs A character vector or list of strings representing paths to Velocyto directories.
#'                      Each directory should contain subdirectories (`filtered` or `raw`) with required files.
#' @param sample_ids A character vector or list of unique sample IDs corresponding to each directory in `velocyto_dirs`.
#' @param whitelist_barcodes A list of character vectors, each containing barcode whitelist(s) for the corresponding sample.
#'                            If `NULL` (default), the function uses the filtered barcodes file if `use_filtered` is `TRUE`.
#' @param use_filtered A logical flag (default `TRUE`) indicating whether to use the `filtered` data for barcode filtration.
#'                     If `FALSE`, the `raw` data is used, and no barcode filtration is applied unless a whitelist is provided.
#' @param merge_counts A logical flag (default `FALSE`) indicating whether to merge all spliced and unspliced matrices
#'                     across samples into two large matrices. If `TRUE`, the function returns a single spliced and
#'                     unspliced matrix instead of individual matrices for each sample.
#'
#' @return A list containing processed Velocyto data. If `merge_counts = FALSE`, the list contains spliced and unspliced 
#'         matrices for each sample. If `merge_counts = TRUE`, the list contains two large combined matrices:
#'         - `spliced`: A single large matrix combining all spliced counts.
#'         - `unspliced`: A single large matrix combining all unspliced counts.
#'
#' @examples
#' # Process data for multiple samples without merging
#' velocyto_dirs <- list("path/to/sample1", "path/to/sample2")
#' sample_ids <- c("sample1", "sample2")
#' result <- process_velocyto_data(velocyto_dirs, sample_ids, merge_counts = FALSE)
#'
#' # Process data for multiple samples and merge counts
#' result <- process_velocyto_data(velocyto_dirs, sample_ids, merge_counts = TRUE)
#'
#' # Process a single sample with a provided whitelist
#' velocyto_dir <- "path/to/sample1"
#' sample_id <- "sample1"
#' whitelist_barcode <- c("barcode1", "barcode2")
#' result <- process_velocyto_data(list(velocyto_dir), list(sample_id), list(whitelist_barcode), merge_counts = FALSE)
#'
#' @export
process_velocyto_data <- function(velocyto_dirs, sample_ids, whitelist_barcodes = NULL, use_filtered = TRUE, merge_counts = FALSE) {
  
  # Handle single sample input by converting to lists
  if (!is.list(velocyto_dirs)) velocyto_dirs <- list(velocyto_dirs)
  if (!is.list(sample_ids)) sample_ids <- list(sample_ids)
  
  # If whitelist_barcodes is NULL, create a list of NULL values
  if (is.null(whitelist_barcodes)) {
    whitelist_barcodes <- vector("list", length(velocyto_dirs))
  }
  
  # Ensure inputs have the same length
  if (!all(length(velocyto_dirs) == length(whitelist_barcodes), length(whitelist_barcodes) == length(sample_ids))) {
    stop("All input lists (velocyto_dirs, whitelist_barcodes, sample_ids) must have the same length.", call. = FALSE)
  }
  
  # Helper function to process one sample
  process_sample <- function(velocyto_dir, sample_id, whitelist_barcode) {
    # Determine directory (filtered or raw)
    data_dir <- if (use_filtered) paste0(velocyto_dir, "/filtered") else paste0(velocyto_dir, "/raw")
    
    # Define paths
    spliced_dir <- paste0(data_dir, "/spliced.mtx")
    unspliced_dir <- paste0(data_dir, "/unspliced.mtx")
    barcodes_dir <- paste0(data_dir, "/barcodes.tsv")
    features_dir <- paste0(data_dir, "/features.tsv")
    filtered_barcodes_dir <- paste0(velocyto_dir, "/filtered/barcodes.tsv")  # For barcode filtration
    
    # Check for required files
    if (!file.exists(spliced_dir)) stop("No spliced matrix file found for sample: ", sample_id, call. = FALSE)
    if (!file.exists(unspliced_dir)) stop("No unspliced matrix file found for sample: ", sample_id, call. = FALSE)
    if (!file.exists(barcodes_dir)) stop("No barcodes file found for sample: ", sample_id, call. = FALSE)
    if (!file.exists(features_dir)) stop("No features file found for sample: ", sample_id, call. = FALSE)
    
    # Read Velocyto data
    cat("├── Processing Velocyto data for sample: ", sample_id, "\n")
    spliced_mtx <- Matrix::readMM(spliced_dir)
    unspliced_mtx <- Matrix::readMM(unspliced_dir)
    stab_barcode <- data.table::fread(barcodes_dir, header = FALSE)
    stab_features <- data.table::fread(features_dir, header = FALSE)
    
    # Set row and column names
    rownames(spliced_mtx) <- stab_features$V1
    rownames(unspliced_mtx) <- stab_features$V1
    colnames(spliced_mtx) <- stab_barcode$V1
    colnames(unspliced_mtx) <- stab_barcode$V1
    
    # Apply barcode filtration
    if (!is.null(whitelist_barcode)) {
      cat("│  ├──  Applying provided whitelist for sample: ", sample_id, "\n")
    } else if (use_filtered && file.exists(filtered_barcodes_dir)) {
      whitelist_barcode <- data.table::fread(filtered_barcodes_dir, header = FALSE)$V1
      cat("│  ├──  Using filtered barcodes for sample: ", sample_id, "\n")
    }
    
    if (!is.null(whitelist_barcode)) {
      initial_barcodes <- ncol(spliced_mtx)
      
      # Filter both spliced and unspliced matrices
      spliced_mtx <- spliced_mtx[, colnames(spliced_mtx) %in% whitelist_barcode, drop = FALSE]
      unspliced_mtx <- unspliced_mtx[, colnames(unspliced_mtx) %in% whitelist_barcode, drop = FALSE]
      
      final_barcodes <- ncol(spliced_mtx)
      
      # Validate barcode filtration
      if (final_barcodes == 0) {
        stop("All barcodes were removed after filtration for sample: ", sample_id, call. = FALSE)
      }
      
      cat("│  ├──  Filtered barcodes applied for sample: ", sample_id, 
          " (Remaining barcodes: ", final_barcodes, ")\n")
    } else {
      cat("│  ├──  No barcode filtration applied for sample: ", sample_id, "\n")
    }
    
    # Append sample ID to column names
    colnames(spliced_mtx) <- paste0(colnames(spliced_mtx), "-", sample_id)
    colnames(unspliced_mtx) <- paste0(colnames(unspliced_mtx), "-", sample_id)
    
    # Collect summary information
    summary_row <- data.frame(
      Sample = sample_id,
      Barcodes = ncol(spliced_mtx),
      Genes = nrow(spliced_mtx),
      stringsAsFactors = FALSE
    )
    
    # Return processed matrices and summary
    cat("│  └── Finished processing Velocyto data for sample: ", sample_id, "\n")
    return(list(
      spliced = as(spliced_mtx, "CsparseMatrix"),
      unspliced = as(unspliced_mtx, "CsparseMatrix"),
      summary = summary_row
    ))
  }
  
  # Use mapply to process all samples
  results <- mapply(
    process_sample,
    velocyto_dirs,
    sample_ids,
    whitelist_barcodes,
    SIMPLIFY = FALSE
  )
  
  # Extract summaries
  summary_table <- do.call(rbind, lapply(results, function(x) x$summary))
  
  # Merge spliced and unspliced matrices if merge_counts is TRUE
  if (merge_counts) {
    spliced_combined <- do.call(cbind, lapply(results, function(x) x$spliced))
    unspliced_combined <- do.call(cbind, lapply(results, function(x) x$unspliced))
    
    cat("\nSummary of Merged Velocyto Counts:\n")
    print(summary_table)
    
    return(list(
      spliced = spliced_combined,
      unspliced = unspliced_combined
    ))
  }
  
  # Print summary table if not merging
  cat("\nSummary of Processed Samples:\n")
  print(summary_table)
  
  # Return individual results if not merging
  final_results <- lapply(results, function(x) list(spliced = x$spliced, unspliced = x$unspliced))
  if (length(final_results) == 1) {
    return(final_results[[1]])  # Return the single result as a list of spliced/unspliced
  } else {
    names(final_results) <- sample_ids
    return(final_results)  # Return results for multiple samples
  }
}
