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
  if (!all(length(STARsolo_SJ_dirs) == length(white_barcode_lists), 
           length(white_barcode_lists) == length(sample_ids))) {
    stop("All input lists (STARsolo_SJ_dirs, white_barcode_lists, sample_ids) must have the same length.", call. = FALSE)
  }
  
  # Helper function to process one sample
  process_sj_sample <- function(STARsolo_SJ_dir, white_barcode_list, sample_id) {
    # Define paths
    mtx_dir <- file.path(STARsolo_SJ_dir, "raw", "matrix.mtx")
    feature_dir <- file.path(STARsolo_SJ_dir, "../../SJ.out.tab")
    barcodes_dir <- file.path(STARsolo_SJ_dir, "raw", "barcodes.tsv")
    internal_whitelist_dir <- file.path(STARsolo_SJ_dir, "..", "Gene", "filtered", "barcodes.tsv")
    
    # Check for required files
    if (!file.exists(mtx_dir)) {
      stop("No abundance matrix in STARsolo SJ directory for sample: ", sample_id, call. = FALSE)
    }
    if (!file.exists(feature_dir)) {
      stop("No feature matrix in STARsolo SJ directory for sample: ", sample_id, call. = FALSE)
    }
    if (!file.exists(barcodes_dir)) {
      stop("No barcode file in STARsolo SJ directory for sample: ", sample_id, call. = FALSE)
    }
    
    # Read splicing data
    cat("├── Processing sample: ", sample_id, "\n")
    mtx <- Matrix::readMM(mtx_dir)
    raw_brc <- data.table::fread(barcodes_dir, header = FALSE, showProgress = FALSE)
    feature <- data.table::fread(
      feature_dir,
      select = c(1, 2, 3, 4, 5, 6),
      col.names = c('chr', 'start', 'end', 'strand', "intron_motif", 'is_annot'),
      showProgress = FALSE
    )
    
    # Use internal whitelist if enabled and no external whitelist is provided
    if (use_internal_whitelist && is.null(white_barcode_list)) {
      if (file.exists(internal_whitelist_dir)) {
        white_barcode_list <- data.table::fread(internal_whitelist_dir, header = FALSE, showProgress = FALSE)$V1
        cat("│  ├── Using STARsolo internal whitelist for sample: ", sample_id, "\n")
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
      
      cat("│  ├── Trimmed junction abundance matrix for sample: ", sample_id, 
          " (", final_barcodes, " barcodes remaining)\n")
    } else {
      cat("│  ├── No barcode filtration applied for sample: ", sample_id, "\n")
    }
    
    # Append sample_id to column names
    colnames(mtx) <- paste0(colnames(mtx), "-", sample_id)
    
    # Save splicing modality
    m1 <- list(
      eventdata = feature,
      junction_ab = as(mtx, "CsparseMatrix")
    )
    
    # Collect summary information
    summary_row <- data.frame(
      Sample   = sample_id,
      Barcodes = ncol(mtx),
      Events   = nrow(mtx),
      stringsAsFactors = FALSE
    )
    
    cat("│  └── Finished processing sample: ", sample_id, "\n")
    return(list(result = m1, summary = summary_row))
  }
  
  # Process all samples
  results <- mapply(
    process_sj_sample,
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
  
  # Always return a named list, even for a single sample
  names(final_results) <- unlist(sample_ids)
  return(final_results)
}



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
  
  # loading required libraries
  load_libraries_for_m1 <- function() {
    # List of required libraries
    libraries <- c("data.table", "dplyr", "Matrix")
    
    # Attempt to load each library and handle errors
    results <- sapply(libraries, function(lib) {
      if (!require(lib, character.only = TRUE)) {
        stop(paste("The package", lib, "is not installed. Please install it using install.packages(\"", lib, "\").", sep = " "))
      }
    })
  }
  
  suppressPackageStartupMessages(load_libraries_for_m1())
  
  # Remove `sample_id` and deduplicate events
  all_junctions <- unique(all_in_one_eventdata[, row_names_mtx])
  all_in_one_eventdata[, sample_id := NULL]
  temp_eventdata <- all_in_one_eventdata[match(all_junctions, row_names_mtx), ]
  
  # grouping based on first and last coordinates
  temp_eventdata_grouped <- temp_eventdata %>% 
    group_by(start_cor_id) %>% 
    mutate(start_cor_group_id = cur_group_id()) %>% 
    mutate(start_cor_group_count = n()) %>% 
    ungroup() %>% 
    group_by(end_cor_id) %>% 
    mutate(end_cor_group_id = cur_group_id()) %>% 
    mutate(end_cor_group_count = n()) %>%
    as.data.table()
  
  # creating the event data for start coordinates groups
  temp_eventdata_grouped_start <- data.table::copy(temp_eventdata_grouped)
  temp_eventdata_grouped_start <- temp_eventdata_grouped_start[,c("index" ,"end_cor_group_id", "end_cor_group_count") := NULL]
  temp_eventdata_grouped_start <- temp_eventdata_grouped_start[start_cor_group_count != 1,]
  temp_eventdata_grouped_start[, start_cor_group_id := paste0(start_cor_group_id, "_S")]
  temp_eventdata_grouped_start[, row_names_mtx_new := paste0(row_names_mtx, "_S")]
  setnames(x = temp_eventdata_grouped_start, 
           old = c("row_names_mtx", "row_names_mtx_new"), 
           new = c("raw_row_names_mtx", "row_names_mtx"))
  
  # creating the event data for end coordinates groups
  temp_eventdata_grouped_end <- data.table::copy(temp_eventdata_grouped)
  temp_eventdata_grouped_end <- temp_eventdata_grouped_end[, c("index" ,"start_cor_group_id", "start_cor_group_count") := NULL]
  temp_eventdata_grouped_end <- temp_eventdata_grouped_end[end_cor_group_count != 1,]
  temp_eventdata_grouped_end[, end_cor_group_id := paste0(end_cor_group_id, "_E")]
  temp_eventdata_grouped_end[, row_names_mtx_new := paste0(row_names_mtx, "_E")]
  setnames(x = temp_eventdata_grouped_end, 
           old = c("row_names_mtx", "row_names_mtx_new"), 
           new = c("raw_row_names_mtx", "row_names_mtx"))
  

  # Combine start and end groups
  if(ncol(temp_eventdata_grouped_end) == ncol(temp_eventdata_grouped_start)){
    
    colnames(temp_eventdata_grouped_end) <- colnames(temp_eventdata_grouped_start)
    eventdata <- rbind(temp_eventdata_grouped_start, temp_eventdata_grouped_end)
    
  }else{
    
    stop("The eventdata for start and end grouping have not the equal numbner of columns")
  }

  # Rename columns for clarity
  setnames(x = eventdata, 
           old = c("start_cor_group_id", "start_cor_group_count"), 
           new = c("group_id", "group_count"))
  
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

#' multigedi_make_gene
#'
#' This function processes gene expression data from a given directory and creates a sparse matrix
#' for gene expression. It supports barcode filtration using a provided whitelist or the filtered barcodes file.
#'
#' @param expression_dirs A character vector or list of strings representing paths to directories containing
#'                        the gene expression matrix (`matrix.mtx`), barcodes, and features files.
#' @param sample_ids A character vector or list of unique sample IDs corresponding to each directory in `expression_dirs`.
#' @param whitelist_barcodes A list of character vectors, each containing barcode whitelist(s) for the corresponding sample.
#'                            If `NULL` (default), the function uses the filtered barcodes file if available.
#' @param use_internal_whitelist A logical flag (default `TRUE`) indicating whether to use the `filtered` data for barcode filtration.
#'                     If `FALSE`, no barcode filtration is applied unless a whitelist is provided.
#' @return A list containing processed gene expression data for each sample. If a single sample is provided,
#'         the function returns the processed data as a sparse matrix. For multiple samples, a named list is returned.
#'
#' @examples
#' # Single sample processing with external whitelist
#' expression_dir <- "path/to/sample1"
#' sample_id <- "sample1"
#' whitelist_barcode <- c("barcode1", "barcode2")
#' result <- multigedi_make_gene(expression_dir, sample_id, list(whitelist_barcode), use_internal_whitelist = FALSE)
#'
#' # Multiple samples processing with default filtered data
#' expression_dirs <- list("path/to/sample1", "path/to/sample2")
#' sample_ids <- c("sample1", "sample2")
#' result <- multigedi_make_gene(expression_dirs, sample_ids, NULL, use_internal_whitelist = TRUE)
#'
#' @export
multigedi_make_gene <- function(expression_dirs, sample_ids, whitelist_barcodes = NULL, use_internal_whitelist = TRUE) {
  
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
  process_ex_sample <- function(expression_dir, sample_id, whitelist_barcode) {
    # Determine directory (filtered or raw)
    data_dir <- if (use_internal_whitelist) paste0(expression_dir, "/filtered") else paste0(expression_dir, "/raw")
    
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
    g_brc <- data.table::fread(expression_barcodes_dir, header = FALSE, showProgress = FALSE)$V1
    g_feature <- data.table::fread(expression_features_dir, header = FALSE, showProgress = FALSE)
    
    # Set row and column names
    rownames(g_mtx) <- g_feature$V1
    colnames(g_mtx) <- g_brc
    
    # Apply barcode filtration
    if (!is.null(whitelist_barcode)) {
      cat("│  ├──  Applying provided whitelist for sample: ", sample_id, "\n")
    } else if (use_internal_whitelist && file.exists(filtered_barcodes_dir)) {
      whitelist_barcode <- data.table::fread(filtered_barcodes_dir, header = FALSE, showProgress = FALSE)$V1
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
    
    # Append sample_id to column names
    colnames(g_mtx) <- paste0(colnames(g_mtx), "-", sample_id)
    
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
    process_ex_sample,
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





#######################################################################
################################ Velocyto #############################
#######################################################################

#' multigedi_make_velo
#'
#' This function processes spliced and unspliced counts data from Velocyto outputs. It reads the respective matrices, 
#' sets appropriate row and column names, applies barcode filtration using a provided whitelist or the filtered barcodes file, 
#' and optionally merges counts across all samples into single large spliced and unspliced matrices.
#'
#' @param velocyto_dirs A character vector or list of strings representing paths to Velocyto directories.
#'                      Each directory should contain subdirectories (`filtered` or `raw`) with required files.
#' @param sample_ids A character vector or list of unique sample IDs corresponding to each directory in `velocyto_dirs`.
#' @param whitelist_barcodes A list of character vectors, each containing barcode whitelist(s) for the corresponding sample.
#'                            If `NULL` (default), the function uses the filtered barcodes file if `use_internal_whitelist` is `TRUE`.
#' @param use_internal_whitelist A logical flag (default `TRUE`) indicating whether to use the `filtered` data for barcode filtration.
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
#' result <- multigedi_make_velo(velocyto_dirs, sample_ids, merge_counts = FALSE)
#'
#' # Process data for multiple samples and merge counts
#' result <- multigedi_make_velo(velocyto_dirs, sample_ids, merge_counts = TRUE)
#'
#' # Process a single sample with a provided whitelist
#' velocyto_dir <- "path/to/sample1"
#' sample_id <- "sample1"
#' whitelist_barcode <- c("barcode1", "barcode2")
#' result <- multigedi_make_velo(list(velocyto_dir), list(sample_id), list(whitelist_barcode), merge_counts = FALSE)
#'
#' @export
multigedi_make_velo <- function(velocyto_dirs, sample_ids, whitelist_barcodes = NULL, use_internal_whitelist = TRUE, merge_counts = FALSE) {
  
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
  process_velo_sample <- function(velocyto_dir, sample_id, whitelist_barcode) {
    # Determine directory (filtered or raw)
    data_dir <- if (use_internal_whitelist) paste0(velocyto_dir, "/filtered") else paste0(velocyto_dir, "/raw")
    
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
    stab_barcode <- data.table::fread(barcodes_dir, header = FALSE, showProgress = FALSE)$V1
    stab_features <- data.table::fread(features_dir, header = FALSE, showProgress = FALSE)
    
    # Set row and column names
    rownames(spliced_mtx) <- stab_features$V1
    rownames(unspliced_mtx) <- stab_features$V1
    colnames(spliced_mtx) <- stab_barcode
    colnames(unspliced_mtx) <- stab_barcode
    
    # Apply barcode filtration
    if (!is.null(whitelist_barcode)) {
      cat("│  ├──  Applying provided whitelist for sample: ", sample_id, "\n")
    } else if (use_internal_whitelist && file.exists(filtered_barcodes_dir)) {
      whitelist_barcode <- data.table::fread(filtered_barcodes_dir, header = FALSE, showProgress = FALSE)$V1
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
    process_velo_sample,
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


#######################################################################
############################## countsplit #############################
#######################################################################

#' multigedi_countsplit
#'
#' Splits a matrix into training and testing datasets using the countsplit method.
#'
#' This function requires the \code{countsplit} package. If it is not installed, the function will attempt to install it.
#'
#' @param m1_inclusion_matrix A numeric matrix to be split.
#' @param folds A positive numeric value specifying the number of folds. Default is 2.
#' @param epsilon A numeric vector of length 2 specifying the epsilon values. Default is \code{c(0.5, 0.5)}.
#' @param object_names A character string specifying the base name for output train/test objects.
#' @examples
#' \dontrun{
#' m1_inclusion_matrix <- matrix(runif(100), nrow = 10)
#' multigedi_countsplit(m1_inclusion_matrix, folds = 2, epsilon = c(0.5, 0.5), object_names = "example")
#' }
#' @import countsplit
#' @export
multigedi_countsplit <- function(m1_inclusion_matrix, folds = 2, epsilon = c(0.5, 0.5), object_names = "m1") {
  # Ensure countsplit package is installed and loaded
  if (!requireNamespace("countsplit", quietly = TRUE)) {
    message("The 'countsplit' package is not installed. Attempting to install it now...")
    install.packages("countsplit")
    if (!requireNamespace("countsplit", quietly = TRUE)) {
      stop("Failed to install the 'countsplit' package. Please install it manually and try again.")
    }
  }
  
  # Load the countsplit library
  library(countsplit)
  
  # Input validation
  if (missing(m1_inclusion_matrix) || is.null(m1_inclusion_matrix)) {
    stop("Error: 'm1_inclusion_matrix' is required and cannot be NULL.")
  }
  
  if (!is.character(object_names) || length(object_names) != 1) {
    stop("Error: 'object_names' must be a single character string.")
  }
  
  if (!is.numeric(folds) || folds <= 0) {
    stop("Error: 'folds' must be a positive numeric value.")
  }
  
  if (!is.numeric(epsilon) || length(epsilon) != 2) {
    stop("Error: 'epsilon' must be a numeric vector of length 2.")
  }
  
  # Try-catch block for better error management
  tryCatch({
    # Perform countsplit operation
    countsplit_obj <- countsplit(X = m1_inclusion_matrix, folds = folds, epsilon = epsilon)
    train_data <- countsplit_obj[[1]]
    test_data <- countsplit_obj[[2]]
    
    # Assign train and test data to global environment
    assign(paste0(object_names, "_train"), train_data, envir = .GlobalEnv)
    assign(paste0(object_names, "_test"), test_data, envir = .GlobalEnv)
    
    message(sprintf(
      "Train and test datasets have been assigned to R objects: '%s_train' and '%s_test'.",
      object_names, object_names
    ))
    
    message("Countsplit operation completed successfully!")
  }, error = function(e) {
    # Handle errors
    stop(paste("An error occurred during countsplit:", e$message))
  })
}




#######################################################################
################################ make M2 ##############################
#######################################################################

#' multigedi_make_m2
#'
#' Creates the M2 matrix from a given m1_inclusion_matrix and eventdata, ensuring the proper processing of group indices and matrix operations.
#'
#' @param m1_inclusion_matrix A sparse matrix to be modified and used for creating the M2 matrix.
#' @param eventdata A data.table containing event information with at least `group_id` and an index column.
#' @return A sparse matrix M2 with the dummy row removed and proper adjustments made.
#' @examples
#' library(data.table)
#' library(Matrix)
#'
#' m1_inclusion_matrix <- Matrix(data = runif(20), nrow = 5, ncol = 4, sparse = TRUE)
#' eventdata <- data.table(group_id = c("A", "B", "C", "D", "E"))
#' M2 <- multigedi_make_m2(m1_inclusion_matrix, eventdata)
#' @export
multigedi_make_m2 <- function(m1_inclusion_matrix, eventdata) {
  # Input validation
  if (missing(m1_inclusion_matrix) || !inherits(m1_inclusion_matrix, "sparseMatrix")) {
    stop("Error: 'm1_inclusion_matrix' must be a sparse matrix and cannot be NULL.")
  }
  if (missing(eventdata) || !is.data.table(eventdata)) {
    stop("Error: 'eventdata' must be a data.table and cannot be NULL.")
  }
  if (!"group_id" %in% colnames(eventdata)) {
    stop("Error: 'eventdata' must contain a 'group_id' column.")
  }
  

  # Add an index column to eventdata
  eventdata[, i := .I]
  
  # Create a dummy row and append to m1_inclusion_matrix
  dummy <- Matrix(data = 1, ncol = ncol(m1_inclusion_matrix), nrow = 1, sparse = TRUE, dimnames = list("dummy", colnames(m1_inclusion_matrix)))
  m1_inclusion_matrix <- rbind(m1_inclusion_matrix, dummy)
  
  message("Step 1 | Modifying the m1_inclusion_matrix")
  
  # Add dummy group to group_ids
  dummy_group <- data.table(i = nrow(m1_inclusion_matrix), group_id = 'dummy')
  group_ids <- eventdata[, .(i, group_id)]
  group_ids <- rbind(group_ids, dummy_group)
  
  rm(dummy_group)  # Remove intermediate variable
  
  # Add group index and initialize variables
  num_cells <- ncol(m1_inclusion_matrix)
  groups_start_vector <- eventdata[, unique(group_id)]
  
  message("Step 2 | Creating M2")
  
  # Convert m1_inclusion_matrix to data.table
  m1 <- summary(m1_inclusion_matrix) %>% as.data.table()
  setnames(m1, 'x', 'x_1')
  
  # Merge group information
  m1 <- merge(m1, group_ids, by = 'i')
  m1[, x_tot := sum(x_1), .(group_id, j)]
  m_tot <- m1[, .(group_id, j, x_tot)] %>% unique()
  
  # Filter and merge relevant data
  m_tot <- m_tot[x_tot > 0]
  m_tot <- merge(m_tot, group_ids, by = 'group_id', allow.cartesian = TRUE)
  m_tot <- merge(m_tot, m1, by = c('group_id', 'i', 'j', 'x_tot'), all.x = TRUE)
  m_tot[is.na(x_1), x_1 := 0]
  m_tot[, x_2 := x_tot - x_1]
  
  # Create sparse matrix for M2_train
  M2_train <- m_tot[, sparseMatrix(i = i, j = j, x = x_2)]
  
  message("Step 3 | Finalizing M2 creation")
  
  # Set row and column names
  rownames(M2_train) <- rownames(m1_inclusion_matrix)
  colnames(M2_train) <- colnames(m1_inclusion_matrix)
  
  # Remove dummy row from M2_train
  M2 <- M2_train[-nrow(M2_train), ]
  
  message("All done!")
  return(M2)
}



##############################################################################
################################ get_embeddings ##############################
##############################################################################

#' @title Compute Embeddings with SVD, PC, and UMAP
#'
#' @description
#' `multigedi_get_embeddings` performs SVD, computes Principal Components,
#' and optionally runs UMAP on the given GEDI model. The function is designed
#' to handle both the joint modality and individual modalities defined in the
#' `GEDI_model`.
#'
#' @param GEDI_model A list or object containing GEDI model data. Must contain
#' `GEDI_model$aux$Modalities` to identify modalities.
#' @param SVD Logical. If `TRUE`, SVD (via `svd.multigedi`) will be performed. Default is `TRUE`.
#' @param PC Logical. If `TRUE`, principal components will be computed from the SVD result. Default is `TRUE`.
#' @param UMAP Logical. If `TRUE`, UMAP (via `uwot::umap`) will be performed. Default is `TRUE`.
#' @param umap_min_res Numeric. Minimum distance parameter for UMAP. Default is `0.01`.
#' @param umap_method Character. The distance metric for UMAP. Default is `"euclidean"`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is `FALSE`.
#'
#' @return
#' A named list of embeddings for each modality (including `Joint`).
#' Each modality contains the results of SVD, PC, and optionally UMAP.
#'
#' @examples
#' \donttest{
#' # Example usage
#' # Suppose `my_gedi_model` is a valid GEDI model object
#' # results <- multigedi_get_embeddings(my_gedi_model, verbose = TRUE)
#' }
#'
#' @export
multigedi_get_embeddings <- function(
    GEDI_model, 
    SVD          = TRUE, 
    PC           = TRUE, 
    UMAP         = TRUE, 
    umap_min_res = 0.01, 
    umap_method  = "euclidean",
    verbose      = FALSE
) {
  # Check input validity
  if (!"aux" %in% names(GEDI_model) || !"Modalities" %in% names(GEDI_model$aux)) {
    stop("Invalid GEDI_model: Missing 'aux$Modalities'.")
  }
  
  # Get all modalities including the joint modality
  modality_names <- GEDI_model$aux$Modalities
  joint_modality <- "Joint"
  all_modalities <- c(joint_modality, modality_names)
  
  # Initialize variables
  embedding_list <- list()
  
  # Iterate through each modality (joint first, then specific modalities)
  for (modality_name in all_modalities) {
    if (verbose) message("Processing modality: ", modality_name)
    
    # Step 1: Perform SVD
    if (SVD) {
      if (verbose) message("  Performing SVD...")
      svd_res <- tryCatch(
        svd.multigedi(
          object     = GEDI_model, 
          Modalities = if (modality_name == joint_modality) NULL else modality_name
        ),
        error = function(e) {
          stop("SVD failed for modality ", modality_name, ": ", e$message)
        }
      )
      embedding_list[[modality_name]][["SVD"]] <- svd_res
      if (verbose) message("  SVD completed.")
    } else {
      stop("SVD is required but was skipped. Cannot proceed.")
    }
    
    # Step 2: Compute Principal Components (PC)
    if (PC) {
      if (!exists("svd_res")) {
        stop("PC computation requires SVD results, which are missing.")
      }
      if (verbose) message("  Computing Principal Components...")
      pc_res <- svd_res$v %*% diag(svd_res$d)
      embedding_list[[modality_name]][["PC"]] <- pc_res
      if (verbose) message("  PC computation completed.")
    } else {
      stop("PC is required but was skipped. Cannot proceed.")
    }
    
    # Step 3: Perform UMAP
    if (UMAP) {
      if (!exists("pc_res")) {
        stop("UMAP computation requires PC results, which are missing.")
      }
      if (verbose) message("  Performing UMAP...")
      umap_res <- tryCatch(
        uwot::umap(
          X       = pc_res, 
          metric  = umap_method, 
          min_dist = umap_min_res
        ),
        error = function(e) {
          stop("UMAP failed for modality ", modality_name, ": ", e$message)
        }
      )
      embedding_list[[modality_name]][["UMAP"]] <- umap_res
      if (verbose) message("  UMAP completed.")
    } else {
      if (verbose) message("  UMAP was skipped.")
    }
  }
  
  if (verbose) message("All modalities processed successfully.")
  return(embedding_list)
}



#' Create Gene-Event Data by Overlapping Event Data with GTF Annotation
#'
#' @description
#' This function reads in a GTF file, extracts gene annotations, and merges them with
#' event-level genomic intervals provided by the user. The final data table contains the
#' original event intervals and the corresponding gene information (for example, `gene_id`
#' and `gene_name`).
#'
#' @details
#' 1. Read the GTF: Uses `rtracklayer::readGFF()` to load GTF data and convert it to
#'    a data table.
#' 2. Subset for Genes: Keeps only rows where `type == "gene"`, retaining columns
#'    for chromosome, start, end, strand, gene_id, and gene_name.
#' 3. Strand Conversion: Merges the GTF data with a small lookup table to replace
#'    `+` and `-` with numeric strand indicators `1` and `2` (matching STAR).
#' 4. Overlaps: With both data sets keyed, uses `foverlaps()` from `data.table` to
#'    find intervals in `eventdata` that fall fully within gene boundaries.
#'
#' @param eventdata A data table of genomic intervals for events. Must contain columns:
#'   * `chr`: Chromosome name (e.g., "chr1").
#'   * `start`: Start coordinate of the event.
#'   * `end`: End coordinate of the event.
#'   * `strand`: Numeric strand indicator (`1` or `2`).
#' @param GTF_file_direction A character string specifying the path to a GTF file. The file
#'   must contain at least these columns: `seqid`, `start`, `end`, `strand`, `gene_id`,
#'   and `gene_name`.
#'
#' @return A data table containing overlapping event intervals with added gene metadata
#'   (such as `gene_id` and `gene_name`). The columns returned will include both event-level
#'   and gene-level information.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' library(rtracklayer)
#'
#' ed <- data.table(
#'   chr = c("chr1", "chr1"),
#'   start = c(1000, 5000),
#'   end = c(2000, 6000),
#'   strand = c(1, 2)
#' )
#'
#' gtf_path <- "path/to/genes.gtf"
#' result_dt <- multigedi_make_eventdata_plus(ed, gtf_path)
#' head(result_dt)
#' }
#'
#' @export
multigedi_make_eventdata_plus <- function(eventdata, GTF_file_direction) {
  
  # Ensure necessary libraries are installed
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    BiocManager::install("rtracklayer")
  }

  
  # Read GTF file and convert to data.table
  GTF <- rtracklayer::readGFF(GTF_file_direction)
  GTF <- data.table::as.data.table(GTF)
  
  # Filter for 'gene' entries
  ref_gtf <- GTF[GTF$type == "gene", ]
  
  # Select needed columns
  ref_gtf <- ref_gtf[, .(seqid, start, end, strand, gene_id, gene_name)]
  
  # Rename seqid to chr and ensure character type
  data.table::setnames(ref_gtf, "seqid", "chr")
  ref_gtf[, chr := as.character(chr)]
  
  # Add 'chr' prefix if missing
  ref_gtf[!grepl("^chr", chr), chr := paste0("chr", chr)]
  
  # Convert strand from + / - to 1 / 2
  temp_change_strand_dt <- data.table::data.table(new_strand = c(1, 2), strand = c("+", "-"))
  ref_gtf <- base::merge(ref_gtf, temp_change_strand_dt, by = "strand", sort = FALSE)
  
  # Replace old strand with numeric one
  ref_gtf[, strand := NULL]
  data.table::setnames(ref_gtf, "new_strand", "strand")
  
  # Standardize eventdata chr column
  eventdata[, chr := as.character(chr)]
  eventdata[!grepl("^chr", chr), chr := paste0("chr", chr)]
  
  # Set keys for foverlaps
  data.table::setkey(eventdata, chr, strand, start, end)
  data.table::setkey(ref_gtf, chr, strand, start, end)
  
  # Perform overlap
  new_eventdata <- data.table::foverlaps(eventdata, ref_gtf, type = "within")
  new_eventdata <- na.omit(new_eventdata)
  
  return(new_eventdata)
}




##############################################################################
################################ get_deviance ##############################
##############################################################################


###' Calculate Deviance for Inclusion and Exclusion Matrices
###'
###' This function computes the deviance for inclusion and exclusion matrices by calling a precompiled C++ function. 
###' It ensures matrices are properly formatted and compatible before proceeding.
###'
###' @param m1_matrix A matrix representing the inclusion matrix. Rows are events, columns are barcodes.
###' @param m2_matrix A matrix representing the exclusion matrix. Rows are events, columns are barcodes.
###' @param min_row_sum A numeric value specifying the minimum row sum threshold for filtering events. Defaults to 50.
###' @param ... Additional arguments to be passed.
###'
###' @return A data.table containing the events and their corresponding deviance values.
###' @export
multigedi_get_deviance <- function(m1_matrix, m2_matrix, min_row_sum = 50, ...) {
  
  # Load necessary libraries
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The 'data.table' package is required but not installed.")
  }
  
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("The 'Rcpp' package is required but not installed.")
  }
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required but not installed.")
  }
  
  # Compile the source code
  Rcpp::sourceCpp("./src/calcDeviances.cpp")
  
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
  cat("There are", length(libraries), "libraries detected...\n")
  
  # Initialize deviance sum vector
  sum_deviances <- numeric(nrow(m1_matrix))
  names(sum_deviances) <- rownames(m1_matrix)
  
  for (lib in libraries) {
    filter <- which(meta[, ID] == lib)
    M1_sub <- m1_matrix[, filter, drop = FALSE]
    M2_sub <- m2_matrix[, filter, drop = FALSE]
    
    # Calculate deviances using the C++ function
    deviance_values <- tryCatch({
      calcDeviances(M1_sub, M2_sub)
    }, error = function(e) {
      stop("Error in calcDeviances function: ", e$message)
    })
    
    deviance_values <- c(deviance_values)
    names(deviance_values) <- rownames(M1_sub)
    sum_deviances <- sum_deviances + deviance_values
    cat("Calculating the deviances for sample", lib, "has been completed!\n")
  }
  
  rez <- data.table::data.table(events = names(sum_deviances), sum_deviance = as.numeric(sum_deviances))
  return(rez)
}





##############################################################################
################################ fast_row_var ##############################
##############################################################################

#' Calculate Row Variance for a Sparse Matrix
#'
#' @description
#' This function calculates the row variance for a given sparse matrix by calling a C++ function via Rcpp.
#'
#' @param sparse_matrix A sparse matrix (of class \code{Matrix}) for which row variances are computed.
#' @param return_vector Logical. If \code{TRUE}, returns a numeric vector of variances; if \code{FALSE}, returns a \code{data.table} with row names and variances. Default is \code{TRUE}.
#' @param ... Additional parameters (currently not used).
#'
#' @return A numeric vector of row variances if \code{return_vector = TRUE}, otherwise a \code{data.table} with columns \code{features} (row names) and \code{rowVar} (variance values).
#'
#' @details
#' This function sources the C++ code from \code{./src/row_variance.cpp} using \code{Rcpp::sourceCpp} and calls the C++ function \code{get_row_variance}. Ensure that the C++ source file exists and contains the exported function using the \code{// [[Rcpp::export]]} tag.
#'
#' @examples
#' \dontrun{
#'   library(Matrix)
#'   # Create a sparse matrix
#'   m <- Matrix::rsparsematrix(100, 10, density = 0.1)
#'   # Calculate row variances as a vector
#'   row_var <- multigedi_get_row_variance(m)
#'
#'   # Calculate row variances as a data.table
#'   row_var_dt <- multigedi_get_row_variance(m, return_vector = FALSE)
#' }
#'
#' @export
multigedi_get_row_variance <- function(sparse_matrix, return_vector = TRUE, ...) {
  
  # Error control: Check for necessary package and class
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("The 'Rcpp' package is required but not installed.")
  }
  
  if (!inherits(sparse_matrix, "Matrix")) {
    stop("Input 'sparse_matrix' must be of class 'Matrix'.")
  }
  
  # Define the C++ source file path
  cpp_source_path <- "./src/row_variance.cpp"
  if (!file.exists(cpp_source_path)) {
    stop("The C++ source file does not exist at the specified path: ", cpp_source_path)
  }
  
  # Source the C++ code with error handling
  tryCatch({
    Rcpp::sourceCpp(cpp_source_path, verbose = FALSE)
  }, error = function(e) {
    stop("Error sourcing the C++ code from '", cpp_source_path, "': ", e$message)
  })
  
  # Call the C++ function to compute row variances
  rez <- tryCatch({
    get_row_variance(M = sparse_matrix)
  }, error = function(e) {
    stop("Error calling the C++ function 'get_row_variance': ", e$message)
  })
  
  # Return the result in the desired format
  if (return_vector) {
    return(rez[, 1])
  } else {
    if (is.null(rownames(sparse_matrix))) {
      stop("Row names are missing in 'sparse_matrix'. They are required when 'return_vector' is FALSE.")
    }
    rez_dt <- data.table::data.table(features = rownames(sparse_matrix), rowVar = rez[, 1])
    return(rez_dt)
  }
}





##################################################################################
################################ find varible genes ##############################
##################################################################################

#' Find Variable Genes Using Deviance and Variance Metrics
#'
#' @description
#' This function identifies variable genes from a sparse gene expression matrix.
#' It provides two methods: a Seurat-like standardization approach or a deviance-based method
#' calculated per library. The deviance method includes filtering low-count genes, calculating
#' negative binomial deviances per library, and combining these with a row variance metric computed
#' via a C++ function.
#'
#' @param gene_expression_matrix A sparse gene expression matrix (of class \code{Matrix}) with gene names as row names.
#' @param like_seurat Logical. If \code{TRUE}, uses a Seurat-like method (via \code{standardizeSparse_variance_loess}).
#'   If \code{FALSE} (default), uses a deviance-based approach.
#' @param ... Additional parameters (currently not used).
#'
#' @return A \code{data.table} containing the gene names (column \code{events}) and the computed metrics.
#'   For the deviance method, the output contains columns \code{sum_deviance} and \code{variance}.
#'
#' @details
#' When \code{like_seurat = TRUE}, the function calls the C++ function \code{standardizeSparse_variance_loess}
#' to compute a standardized variance vector for the genes. When \code{like_seurat = FALSE}, the function:
#' \enumerate{
#'   \item Filters out genes with zero row sum.
#'   \item Parses library IDs from the column names of the gene expression matrix.
#'   \item Loops over each library and computes negative binomial deviances using the C++ function
#'         \code{calcNBDeviancesWithThetaEstimation}.
#'   \item Computes row variance using \code{multigedi_get_row_variance}.
#'   \item Merges the deviance and variance information into a single \code{data.table}.
#' }
#'
#' The function sources required C++ files \code{deviance_gene.cpp} and \code{hvf_gene_expression.cpp}.
#' Ensure that these files exist in the \code{./src/} directory and contain the appropriate
#' \code{// [[Rcpp::export]]} tags.
#'
#' @examples
#' \dontrun{
#'   library(Matrix)
#'   library(data.table)
#'   # Create a sample sparse gene expression matrix
#'   gene_expression_matrix <- Matrix::rsparsematrix(100, 10, density = 0.1)
#'   # Use the deviance-based method
#'   result <- multigedi_find_variable_genes(gene_expression_matrix)
#'   # Alternatively, use the Seurat-like method
#'   result_seurat <- multigedi_find_variable_genes(gene_expression_matrix, like_seurat = TRUE)
#' }
#'
#' @export
multigedi_find_variable_genes <- function(gene_expression_matrix, like_seurat = FALSE, ...) {
  
  # Check required libraries
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The 'data.table' package is required but not installed.")
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("The 'Rcpp' package is required but not installed.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required but not installed.")
  }
  
  # Verify that gene_expression_matrix is a sparse Matrix
  if (!inherits(gene_expression_matrix, "Matrix")) {
    stop("The 'gene_expression_matrix' must be a sparse matrix of class 'Matrix'.")
  }
  
  # Source necessary C++ files for deviance and gene expression functions
  cpp_files <- c("./src/deviance_gene.cpp", "./src/hvf_gene_expression.cpp")
  for (cpp_file in cpp_files) {
    if (!file.exists(cpp_file)) {
      stop("The C++ source file does not exist: ", cpp_file)
    }
    tryCatch({
      Rcpp::sourceCpp(cpp_file, verbose = FALSE)
    }, error = function(e) {
      stop("Error sourcing the C++ file '", cpp_file, "': ", e$message)
    })
  }
  
  if (like_seurat) {
    cat("The method we are using is like Seurat...\n")
    if (!exists("standardizeSparse_variance_loess")) {
      stop("The function 'standardizeSparse_variance_loess' is not available. Check your C++ source files.")
    }
    rez_vector <- tryCatch({
      standardizeSparse_variance_loess(X = gene_expression_matrix)
    }, error = function(e) {
      stop("Error in standardizeSparse_variance_loess: ", e$message)
    })
    rez <- data.table::data.table(events = rownames(gene_expression_matrix),
                                  standardize_variance = rez_vector)
  } else {
    cat("The method we are using is like deviance summarion per library...\n")
    
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
    cat("There are", length(libraries), "libraries detected...\n")
    
    # Initialize deviance sum vector with gene names
    sum_deviances <- numeric(nrow(gene_expression_matrix))
    names(sum_deviances) <- rownames(gene_expression_matrix)
    
    # Loop over each library to compute deviances
    for (lib in libraries) {
      filter <- which(meta[, ID] == lib)
      gene_expression_matrix_sub <- gene_expression_matrix[, filter, drop = FALSE]
      
      # Calculate deviances using the C++ function
      deviance_values <- tryCatch({
        calcNBDeviancesWithThetaEstimation(gene_expression_matrix_sub)
      }, error = function(e) {
        stop("Error in calcNBDeviancesWithThetaEstimation function: ", e$message)
      })
      
      deviance_values <- c(deviance_values)
      names(deviance_values) <- rownames(gene_expression_matrix_sub)
      sum_deviances <- sum_deviances + deviance_values
      cat("Calculating the deviances for sample", lib, "has been completed!\n")
    }
    
    # Compute row variance using the previously defined function
    row_var <- tryCatch({
      multigedi_get_row_variance(sparse_matrix = gene_expression_matrix)
    }, error = function(e) {
      stop("Error in multigedi_get_row_variance: ", e$message)
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






