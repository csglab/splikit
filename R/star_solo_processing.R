#' make_junction_ab
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
#' @param verbose Logical (default `FALSE`). If `TRUE`, prints detailed progress messages during processing.
#' @param keep_multi_mapped_junctions Logical (default `FALSE`). If `TRUE`, retains multi-mapped junctions (flag 0).
#'                                     If `FALSE`, only keeps uniquely mapped junctions (flag 1).
#' @param ... Additional parameters for future extensions.
#' @return
#' A list containing processed splicing modality data for each sample. If a single sample is provided,
#' the function returns the processed data as a list. For multiple samples, a named list is returned.
#' Each sample's result contains:
#'
#' \describe{
#'   \item{eventdata}{A `data.table` containing feature metadata.}
#'   \item{junction_ab}{A sparse matrix of junction abundance.}
#' }
#'
#' @importFrom Matrix readMM
#' @importFrom data.table fread
#' @export
make_junction_ab <- function(STARsolo_SJ_dirs, white_barcode_lists = NULL, sample_ids, use_internal_whitelist = TRUE, verbose = FALSE, keep_multi_mapped_junctions = FALSE, ...) {

  # Handle single sample input by converting to lists
  if (!is.list(STARsolo_SJ_dirs)) STARsolo_SJ_dirs <- list(STARsolo_SJ_dirs)
  if (!is.list(sample_ids)) sample_ids <- list(sample_ids)

  # Convert white_barcode_lists to list if it's not NULL and not already a list
  if (!is.null(white_barcode_lists) && !is.list(white_barcode_lists)) {
    white_barcode_lists <- list(white_barcode_lists)
  }

  # Automatically set use_internal_whitelist to FALSE if white_barcode_lists is provided
  if (!is.null(white_barcode_lists)) {
    use_internal_whitelist <- FALSE
    if (verbose) message("External whitelist provided: automatically disabling internal whitelist")
  }

  # If white_barcode_lists is NULL, create a list of NULL values
  if (is.null(white_barcode_lists)) {
    white_barcode_lists <- vector("list", length(STARsolo_SJ_dirs))
  }

  # Ensure inputs have the same length
  if (!all(length(STARsolo_SJ_dirs) == length(white_barcode_lists),
           length(white_barcode_lists) == length(sample_ids))) {
    stop("All input lists (STARsolo_SJ_dirs, white_barcode_lists, sample_ids) must have the same length.\n",
         "Lengths: STARsolo_SJ_dirs=", length(STARsolo_SJ_dirs),
         ", white_barcode_lists=", length(white_barcode_lists),
         ", sample_ids=", length(sample_ids), call. = FALSE)
  }

  # Helper function to process one sample
  process_sj_sample <- function(STARsolo_SJ_dir, white_barcode_list, sample_id) {
    # Define paths
    mtx_dir <- file.path(STARsolo_SJ_dir, "raw", "matrix.mtx")
    feature_dir <- file.path(STARsolo_SJ_dir, "../../SJ.out.tab")
    barcodes_dir <- file.path(STARsolo_SJ_dir, "raw", "barcodes.tsv")
    internal_whitelist_dir <- file.path(STARsolo_SJ_dir, "..", "Gene", "filtered", "barcodes.tsv")

    # Debug: Print paths if verbose
    if (verbose) {
      message("|-- Debug paths for sample: ", sample_id)
      message("    Matrix file: ", mtx_dir)
      message("    Feature file: ", feature_dir)
      message("    Barcodes file: ", barcodes_dir)
      message("    Internal whitelist: ", internal_whitelist_dir)
    }

    # Check for required files
    if (!file.exists(mtx_dir)) {
      stop("No abundance matrix in STARsolo SJ directory for sample: ", sample_id,
           "\nLooked at: ", mtx_dir, call. = FALSE)
    }
    if (!file.exists(feature_dir)) {
      stop("No feature matrix in STARsolo SJ directory for sample: ", sample_id,
           "\nLooked at: ", feature_dir, call. = FALSE)
    }
    if (!file.exists(barcodes_dir)) {
      stop("No barcode file in STARsolo SJ directory for sample: ", sample_id,
           "\nLooked at: ", barcodes_dir, call. = FALSE)
    }

    # Read splicing data with error handling
    if (verbose) message("|-- Processing sample: ", sample_id)

    tryCatch({
      mtx <- Matrix::readMM(mtx_dir)
      if (verbose) message("    Matrix dimensions: ", nrow(mtx), "x", ncol(mtx))
    }, error = function(e) {
      stop("Error reading matrix file for sample ", sample_id, ": ", e$message, call. = FALSE)
    })

    tryCatch({
      raw_brc <- data.table::fread(barcodes_dir, header = FALSE, showProgress = FALSE)
      if (verbose) message("    Barcodes read: ", nrow(raw_brc))
    }, error = function(e) {
      stop("Error reading barcodes file for sample ", sample_id, ": ", e$message, call. = FALSE)
    })

    tryCatch({
      feature <- data.table::fread(
        feature_dir,
        select = c(1, 2, 3, 4, 5, 6, 7),
        col.names = c('chr', 'start', 'end', 'strand', "intron_motif", 'is_annot', 'unique_mapped'),
        showProgress = FALSE
      )
      if (verbose) message("    Features read: ", nrow(feature))
    }, error = function(e) {
      stop("Error reading feature file for sample ", sample_id, ": ", e$message, call. = FALSE)
    })

    # Validate dimensions match
    if(nrow(mtx) != nrow(feature)) {
      stop("Dimension mismatch for sample ", sample_id, ": matrix has ", nrow(mtx),
           " rows but feature file has ", nrow(feature), " rows", call. = FALSE)
    }

    if(ncol(mtx) != nrow(raw_brc)) {
      stop("Dimension mismatch for sample ", sample_id, ": matrix has ", ncol(mtx),
           " columns but barcodes file has ", nrow(raw_brc), " rows", call. = FALSE)
    }

    # Use internal whitelist if enabled and no external whitelist is provided
    if (use_internal_whitelist && is.null(white_barcode_list)) {
      if (file.exists(internal_whitelist_dir)) {
        tryCatch({
          white_barcode_list <- data.table::fread(internal_whitelist_dir, header = FALSE, showProgress = FALSE)$V1
          if (verbose) message("    |-- Using STARsolo internal whitelist for sample: ", sample_id,
                               " (", length(white_barcode_list), " barcodes)")
        }, error = function(e) {
          stop("Error reading internal whitelist for sample ", sample_id, ": ", e$message, call. = FALSE)
        })
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

    # Filter multi-mapped junctions if requested
    if(!keep_multi_mapped_junctions) {
      old_junction_numbers <- nrow(feature)

      # Filter for unique mapped junctions only
      feature <- feature[unique_mapped > 0, ]

      # Subset matrix to match filtered features
      if(nrow(feature) > 0) {
        mtx <- mtx[feature$row_names_mtx, , drop = FALSE]
      } else {
        warning("No unique mapped junctions found for sample ", sample_id)
        mtx <- mtx[integer(0), , drop = FALSE]  # Empty matrix with correct structure
      }

      eliminated_junctions <- old_junction_numbers - nrow(feature)
      if (verbose && eliminated_junctions > 0) {
        message("    |-- Eliminated ", eliminated_junctions, " junctions due to only multi-mapped records")
      }
    }

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

      if (verbose) message("|  |-- Trimmed junction abundance matrix for sample: ", sample_id,
                           " (", final_barcodes, " barcodes remaining)")
    } else {
      if (verbose) message("|  |-- No barcode filtration applied for sample: ", sample_id)
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

    if (verbose) message("|  +-- Finished processing sample: ", sample_id)
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
  if (verbose) {
    message("\nSummary of Processed Samples in M1 matrix:")
    print(summary_table)
  }

  # Always return a named list, even for a single sample
  names(final_results) <- unlist(sample_ids)
  return(final_results)
}

#' Load the toy SJ object
#'
#' Loads a toy splice junction object used for examples or testing.
#'
#' @return An R object from the toy_SJ_object.RDS file.
#' @export
load_toy_SJ_object <- function() {
  file <- system.file("extdata", "toy_SJ_object.RDS", package = "splikit")
  readRDS(file)
}

#' make_m1
#'
#' This function processes junction abundance data from multiple samples to create a splicing modality inclusion matrix (M1).
#' It merges event data, handles start and end coordinate groups, ensures matrix compatibility, and includes robust error handling.
#'
#' The function requires the following libraries: `data.table`, and `Matrix`.
#'
#' @param junction_ab_object A named list where each element represents a sample's junction abundance data.
#'                           Each element must contain `eventdata` and a sparse matrix.
#' @param min_counts Numeric (default 1). Minimum count threshold for filtering events. 
#'                   Events with total counts below this threshold will be removed.
#' @param verbose Logical (default `FALSE`). If `TRUE`, prints detailed progress messages during processing.
#' @return
#' A list containing the processed data from all samples:
#'
#' \describe{
#'   \item{m1_inclusion_matrix}{A matrix representing the processed inclusion values for all events across all samples.}
#'   \item{event_data}{A `data.table` containing the merged and grouped metadata for each event.}
#' }
#'
#' @examples
#' # Example usage
#' junction_abundance_object <- load_toy_SJ_object()
#' result <- make_m1(junction_abundance_object)
#' m1_matrix <- result$m1_inclusion_matrix
#' event_metadata <- result$event_data
#'
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom data.table setDT copy setnames := data.table is.data.table .N .GRP
#' @export

make_m1 <- function(junction_ab_object, min_counts = 1, verbose = FALSE) {
  # Validate inputs
  if (!is.list(junction_ab_object) || length(junction_ab_object) == 0) {
    stop("`junction_ab_object` must be a non-empty list.", call. = FALSE)
  }

  if (!is.numeric(min_counts) || length(min_counts) != 1 || min_counts < 0) {
    stop("`min_counts` must be a single non-negative number.", call. = FALSE)
  }

  if (verbose) message("Starting M1 matrix creation...")

  # Extract and combine all eventdata
  all_in_one_eventdata <- tryCatch({
    eventdata_list <- lapply(junction_ab_object, function(x) {
      if (!"eventdata" %in% names(x)) {
        stop("Each element must contain 'eventdata' component.")
      }
      return(x$eventdata)
    })
    do.call(rbind, eventdata_list)
  }, error = function(e) {
    stop("Error combining `eventdata`: ", e$message, call. = FALSE)
  })

  # Validate required columns
  required_cols <- c("row_names_mtx", "start_cor_id", "end_cor_id")
  missing_cols <- setdiff(required_cols, colnames(all_in_one_eventdata))
  if (length(missing_cols) > 0) {
    stop("`eventdata` must contain columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (verbose) message("Combined eventdata from ", length(junction_ab_object), " samples")

  # Remove sample_id and deduplicate events
  all_junctions <- unique(all_in_one_eventdata$row_names_mtx)
  if (verbose) message("Found ", length(all_junctions), " unique junctions")

  # Create a copy and remove sample_id if it exists
  temp_eventdata <- data.table::copy(all_in_one_eventdata)
  if ("sample_id" %in% colnames(temp_eventdata)) {
    temp_eventdata[, sample_id := NULL]
  }

  # Get unique events
  temp_eventdata <- temp_eventdata[match(all_junctions, row_names_mtx), ]

  # Convert to data.table if not already
  if (!data.table::is.data.table(temp_eventdata)) {
    data.table::setDT(temp_eventdata)
  }

  if (verbose) message("Creating coordinate groups...")

  # Group by start and end coordinates
  tryCatch({
    temp_eventdata[, `:=`(
      start_cor_group_id = .GRP,
      start_cor_group_count = .N
    ), by = start_cor_id]

    temp_eventdata[, `:=`(
      end_cor_group_id = .GRP,
      end_cor_group_count = .N
    ), by = end_cor_id]
  }, error = function(e) {
    stop("Error creating coordinate groups: ", e$message, call. = FALSE)
  })

  temp_eventdata_grouped <- temp_eventdata

  if (verbose) {
    message("Start coordinate groups: ", max(temp_eventdata_grouped$start_cor_group_id, na.rm = TRUE))
    message("End coordinate groups: ", max(temp_eventdata_grouped$end_cor_group_id, na.rm = TRUE))
  }

  # Create event data for start coordinate groups
  temp_eventdata_grouped_start <- data.table::copy(temp_eventdata_grouped)
  temp_eventdata_grouped_start <- temp_eventdata_grouped_start[, c("index", "end_cor_group_id", "end_cor_group_count") := NULL]

  # Filter for alternative splicing events only (group count > 1)
  temp_eventdata_grouped_start <- temp_eventdata_grouped_start[start_cor_group_count > 1, ]

  if (nrow(temp_eventdata_grouped_start) > 0) {
    temp_eventdata_grouped_start[, start_cor_group_id := paste0(start_cor_group_id, "_S")]
    temp_eventdata_grouped_start[, row_names_mtx_new := paste0(row_names_mtx, "_S")]
    data.table::setnames(x = temp_eventdata_grouped_start,
                         old = c("row_names_mtx", "row_names_mtx_new"),
                         new = c("raw_row_names_mtx", "row_names_mtx"))
  }

  # Create event data for end coordinate groups
  temp_eventdata_grouped_end <- data.table::copy(temp_eventdata_grouped)
  temp_eventdata_grouped_end <- temp_eventdata_grouped_end[, c("index", "start_cor_group_id", "start_cor_group_count") := NULL]

  # Filter for alternative splicing events only (group count > 1)
  temp_eventdata_grouped_end <- temp_eventdata_grouped_end[end_cor_group_count > 1, ]

  if (nrow(temp_eventdata_grouped_end) > 0) {
    temp_eventdata_grouped_end[, end_cor_group_id := paste0(end_cor_group_id, "_E")]
    temp_eventdata_grouped_end[, row_names_mtx_new := paste0(row_names_mtx, "_E")]
    data.table::setnames(x = temp_eventdata_grouped_end,
                         old = c("row_names_mtx", "row_names_mtx_new"),
                         new = c("raw_row_names_mtx", "row_names_mtx"))
  }

  if (verbose) message("Start coordinate alternative events: ", nrow(temp_eventdata_grouped_start))
  if (verbose) message("End coordinate alternative events: ", nrow(temp_eventdata_grouped_end))

  # Combine start and end groups
  if (nrow(temp_eventdata_grouped_start) == 0 && nrow(temp_eventdata_grouped_end) == 0) {
    stop("No alternative splicing events found (all groups have count = 1)", call. = FALSE)
  }

  if (nrow(temp_eventdata_grouped_start) > 0 && nrow(temp_eventdata_grouped_end) > 0) {
    if (ncol(temp_eventdata_grouped_end) == ncol(temp_eventdata_grouped_start)) {
      colnames(temp_eventdata_grouped_end) <- colnames(temp_eventdata_grouped_start)
      eventdata <- rbind(temp_eventdata_grouped_start, temp_eventdata_grouped_end)
    } else {
      stop("The eventdata for start and end grouping have different numbers of columns", call. = FALSE)
    }
  } else if (nrow(temp_eventdata_grouped_start) > 0) {
    eventdata <- temp_eventdata_grouped_start
  } else {
    eventdata <- temp_eventdata_grouped_end
    # Rename columns to match start naming convention
    data.table::setnames(x = eventdata,
                         old = c("end_cor_group_id", "end_cor_group_count"),
                         new = c("group_id", "group_count"))
  }

  # Rename columns for clarity (only if they exist)
  if ("start_cor_group_id" %in% colnames(eventdata) && "start_cor_group_count" %in% colnames(eventdata)) {
    data.table::setnames(x = eventdata,
                         old = c("start_cor_group_id", "start_cor_group_count"),
                         new = c("group_id", "group_count"))
  }

  if (verbose) message("Combined eventdata has ", nrow(eventdata), " alternative splicing events")

  # Create a table for all events
  all_events <- eventdata[, .(raw_row_names_mtx, row_names_mtx)]
  all_events[, group_type := sub("^.*(.)$", "\\1", row_names_mtx)]

  # Process and merge all junction abundance matrices
  process_junction_ab_matrices <- function(junction_ab_object, all_events, eventdata, min_counts, verbose) {
    if (verbose) message("Processing junction abundance matrices...")

    # Initialize the final merged matrix
    m1_raw <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(nrow(eventdata), 0),
      dimnames = list(eventdata$row_names_mtx, NULL)
    )

    for (j in seq_along(junction_ab_object)) {
      if (verbose) message("Processing sample ", j, " of ", length(junction_ab_object))

      tryCatch({
        # Extract current matrix and events
        little_m1 <- junction_ab_object[[j]]

        if (!"junction_ab" %in% names(little_m1)) {
          stop("Sample ", j, " missing 'junction_ab' component")
        }

        current_matrix <- little_m1[["junction_ab"]]
        current_events <- data.table::data.table(raw_row_names_mtx = rownames(current_matrix))

        # Process start and end events
        start_events <- all_events[group_type == "S"]
        end_events <- all_events[group_type == "E"]

        combined_matrix_parts <- list()

        if (nrow(start_events) > 0) {
          matched_start <- merge(current_events, start_events, by = "raw_row_names_mtx", sort = FALSE)
          if (nrow(matched_start) > 0) {
            start_matrix <- current_matrix[matched_start$raw_row_names_mtx, , drop = FALSE]
            rownames(start_matrix) <- matched_start$row_names_mtx
            combined_matrix_parts[["start"]] <- start_matrix
          }
        }

        if (nrow(end_events) > 0) {
          matched_end <- merge(current_events, end_events, by = "raw_row_names_mtx", sort = FALSE)
          if (nrow(matched_end) > 0) {
            end_matrix <- current_matrix[matched_end$raw_row_names_mtx, , drop = FALSE]
            rownames(end_matrix) <- matched_end$row_names_mtx
            combined_matrix_parts[["end"]] <- end_matrix
          }
        }

        # Combine matrices
        if (length(combined_matrix_parts) > 0) {
          merged_matrix <- do.call(rbind, combined_matrix_parts)
        } else {
          # Create empty matrix with correct dimensions
          merged_matrix <- Matrix::sparseMatrix(
            i = integer(0), j = integer(0),
            dims = c(0, ncol(current_matrix)),
            dimnames = list(character(0), colnames(current_matrix))
          )
        }

        # Add missing events as zero rows
        missing_events <- setdiff(eventdata$row_names_mtx, rownames(merged_matrix))
        if (length(missing_events) > 0) {
          missing_matrix <- Matrix::sparseMatrix(
            i = integer(0), j = integer(0),
            dims = c(length(missing_events), ncol(current_matrix)),
            dimnames = list(missing_events, colnames(current_matrix))
          )
          merged_matrix <- rbind(merged_matrix, missing_matrix)
        }

        # Reorder rows to match eventdata
        final_matrix <- merged_matrix[eventdata$row_names_mtx, , drop = FALSE]

        # Combine into the final matrix
        m1_raw <- cbind(m1_raw, final_matrix)

        # Log progress
        if(verbose) cat("Processed matrix:", names(junction_ab_object)[j], "\n")

      }, error = function(e) {
        stop("Error processing sample ", j, " (", names(junction_ab_object)[j], "): ", e$message, call. = FALSE)
      })
    }

    return(m1_raw)
  }

  # Process all matrices
  m1 <- process_junction_ab_matrices(junction_ab_object, all_events, eventdata, min_counts, verbose)

  if (verbose) message("Applying count threshold filtering...")

  # Filter based on rowSums threshold
  row_sums <- Matrix::rowSums(m1)
  events_to_keep <- row_sums >= min_counts

  if (sum(events_to_keep) == 0) {
    stop("No events pass the minimum count threshold of ", min_counts, call. = FALSE)
  }

  # Filter matrix and eventdata
  m1_filtered <- m1[events_to_keep, , drop = FALSE]
  eventdata_filtered <- eventdata[events_to_keep, ]

  if (verbose) message("Filtered from ", nrow(m1), " to ", nrow(m1_filtered), " events")
  if (verbose) message("Events removed: ", sum(!events_to_keep))

  # Verify matrix row order
  if (!identical(rownames(m1_filtered), eventdata_filtered$row_names_mtx)) {
    m1_filtered <- m1_filtered[eventdata_filtered$row_names_mtx, , drop = FALSE]
  }

  if (verbose) message("Finished processing M1.")

  # Return results with summary
  result <- list(
    m1_inclusion_matrix = m1_filtered,
    event_data = eventdata_filtered,
    summary = list(
      total_events_input = length(all_junctions),
      alternative_events_found = nrow(eventdata),
      events_passing_threshold = nrow(eventdata_filtered),
      min_counts_threshold = min_counts,
      samples_processed = length(junction_ab_object),
      total_cells = ncol(m1_filtered)
    )
  )

  if (verbose) {
    message("\nSummary:")
    message("  Input junctions: ", result$summary$total_events_input)
    message("  Alternative splicing events: ", result$summary$alternative_events_found)
    message("  Events passing threshold: ", result$summary$events_passing_threshold)
    message("  Total cells: ", result$summary$total_cells)
  }

  return(result)
}

#' make_m2 (Integrated with Automatic Batching)
#'
#' Creates the M2 matrix from a given m1_inclusion_matrix and eventdata with intelligent
#' memory management. Automatically detects when the operation would exceed memory limits
#' and switches to a batched sparse matrix approach.
#'
#' @param m1_inclusion_matrix A sparse matrix to be modified and used for creating the M2 matrix.
#' @param eventdata A data.table containing event information with at least `group_id` and an index column.
#' @param batch_size An integer specifying the number of groups to process per batch (default: 5000).
#'   Only used when batched processing is triggered.
#' @param memory_threshold A numeric value representing the maximum number of rows allowed in the
#'   summary before switching to batched processing (default: 2e9, which is ~93% of 2^31).
#' @param force_fast A logical flag to force fast processing regardless of size estimates (default: FALSE).
#'   WARNING: This may cause memory errors on large datasets.
#' @param n_threads Number of threads for parallel processing in batched operations (default: 1).
#'   Only used when batched processing is triggered. Values > 1 require parallel package.
#' @param verbose A logical flag for detailed progress reporting (default: FALSE).
#'
#' @return A sparse matrix M2 with the dummy row removed and proper adjustments made.
#'
#' @examples
#' junction_abundance_object <- load_toy_SJ_object()
#' m1_obj <- make_m1(junction_ab_object = junction_abundance_object)
#'
#' # obtaining the m1 and eventdata
#' m1 <- m1_obj$m1_inclusion_matrix
#' eventdata <- m1_obj$event_data
#' m2 <- make_m2(m1_inclusion_matrix = m1, eventdata = eventdata)
#'
#' @import Matrix
#' @import data.table
#' @importFrom data.table is.data.table := data.table as.data.table
#'
#' @export
make_m2 <- function(m1_inclusion_matrix, eventdata, batch_size = 5000,
                    memory_threshold = 2e9, force_fast = FALSE,
                    n_threads = 1, verbose = FALSE) {

  # Input validation
  if (missing(m1_inclusion_matrix) || !inherits(m1_inclusion_matrix, "sparseMatrix")) {
    stop("Error: 'm1_inclusion_matrix' must be a sparse matrix and cannot be NULL.", call. = FALSE)
  }
  if (missing(eventdata) || !is.data.table(eventdata)) {
    stop("Error: 'eventdata' must be a data.table and cannot be NULL.", call. = FALSE)
  }
  if (!"group_id" %in% colnames(eventdata)) {
    stop("Error: 'eventdata' must contain a 'group_id' column.", call. = FALSE)
  }
  if (!is.numeric(batch_size) || batch_size <= 0) {
    stop("Error: 'batch_size' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(memory_threshold) || memory_threshold <= 0) {
    stop("Error: 'memory_threshold' must be a positive number.", call. = FALSE)
  }
  if (!is.logical(force_fast) || length(force_fast) != 1) {
    stop("Error: 'force_fast' must be a logical value (TRUE or FALSE).", call. = FALSE)
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1 || n_threads < 1) {
    stop("Error: 'n_threads' must be a positive integer.", call. = FALSE)
  }
  n_threads <- as.integer(n_threads)

  # Check for parallel package if multi-threading is requested
  if (n_threads > 1 && !requireNamespace("parallel", quietly = TRUE)) {
    stop("Error: 'parallel' package is required for n_threads > 1 but is not installed.", call. = FALSE)
  }

  if (verbose) message("Starting M2 matrix creation...")

  # Add an index column to eventdata
  eventdata[, i := .I]

  # Create a dummy row and append to m1_inclusion_matrix
  tryCatch({
    dummy <- Matrix::Matrix(data = 1, ncol = ncol(m1_inclusion_matrix),
                            nrow = 1, sparse = TRUE,
                            dimnames = list("dummy", colnames(m1_inclusion_matrix)))
    m1_inclusion_matrix <- rbind(m1_inclusion_matrix, dummy)
  }, error = function(e) {
    stop("Error creating dummy row: ", e$message, call. = FALSE)
  })

  if (verbose) message("+-- Step 1 | Modifying the m1_inclusion_matrix")

  # Add dummy group to group_ids
  dummy_group <- data.table::data.table(i = nrow(m1_inclusion_matrix), group_id = "dummy")
  group_ids <- eventdata[, .(i, group_id)]
  group_ids <- rbind(group_ids, dummy_group)
  rm(dummy_group)

  # Extract basic information for size estimation
  num_barcodes <- ncol(m1_inclusion_matrix)
  num_events <- nrow(eventdata)
  unique_groups <- unique(group_ids$group_id)
  num_groups <- length(unique_groups)

  if (verbose) {
    message("|-- Dataset Information:")
    message("|   |-- Number of barcodes: ", formatC(num_barcodes, format = "d", big.mark = ","))
    message("|   |-- Number of events: ", formatC(num_events, format = "d", big.mark = ","))
    message("|   +-- Number of groups: ", formatC(num_groups, format = "d", big.mark = ","))
  }

  # Estimate memory requirements for the fast approach
  if (verbose) message("+-- Step 2 | Analyzing memory requirements")

  # Get the number of non-zero elements in the sparse matrix
  nnz_elements <- Matrix::nnzero(m1_inclusion_matrix)

  if (verbose) {
    message("|   |-- Non-zero elements in matrix: ", formatC(nnz_elements, format = "d", big.mark = ","))
    message("|   |-- Memory threshold: ", formatC(memory_threshold, format = "d", big.mark = ","), " rows")
  }

  # Estimate the size of operations based on matrix sparsity and group structure
  # This is a conservative estimate based on the summary size and potential Cartesian products

  # Calculate group statistics for better estimation
  tryCatch({
    group_sizes <- group_ids[, .(group_size = .N), by = group_id]
    max_group_size <- max(group_sizes$group_size)
    avg_group_size <- mean(group_sizes$group_size)
  }, error = function(e) {
    stop("Error calculating group statistics: ", e$message, call. = FALSE)
  })

  # Safe multiplication to avoid integer overflow
  # Convert to numeric to handle large numbers safely
  estimated_operations <- tryCatch({
    as.numeric(nnz_elements) * as.numeric(max_group_size)
  }, error = function(e) {
    # If there's still an error, assume it's very large and use batched approach
    Inf
  })

  # Handle NA or overflow cases
  if (is.na(estimated_operations) || is.infinite(estimated_operations)) {
    estimated_operations <- Inf
    if (verbose) message("|   |-- WARNING: Calculation overflow detected - dataset is very large")
  }

  if (verbose) {
    message("|   |-- Maximum group size: ", formatC(max_group_size, format = "d", big.mark = ","))
    message("|   |-- Average group size: ", round(avg_group_size, 2))
    if (is.finite(estimated_operations)) {
      message("|   |-- Estimated operation size: ", formatC(estimated_operations, format = "d", big.mark = ","))
    } else {
      message("|   |-- Estimated operation size: Very large (overflow detected)")
    }
  }

  # Decide on processing approach - use batched if overflow or exceeds threshold
  # unless forced to use fast approach
  should_use_batched <- is.na(estimated_operations) ||
    is.infinite(estimated_operations) ||
    estimated_operations > memory_threshold

  use_batched_approach <- should_use_batched && !force_fast

  if (force_fast && should_use_batched) {
    if (verbose) {
      if (is.finite(estimated_operations)) {
        message("|   |-- WARNING: force_fast=TRUE but estimated size (",
                formatC(estimated_operations, format = "d", big.mark = ","),
                ") exceeds threshold!")
      } else {
        message("|   |-- WARNING: force_fast=TRUE but dataset size suggests overflow risk!")
      }
      message("|   |-- Proceeding with fast approach as requested - monitor memory usage")
    }
  }

  if (use_batched_approach) {
    if (verbose) {
      if (is.finite(estimated_operations)) {
        message("|   |-- Estimated size exceeds memory threshold!")
      } else {
        message("|   |-- Dataset size suggests high overflow risk!")
      }
      message("|   +-- Automatically switching to batched processing approach")
    }

    # Call the batched processing function
    result <- .make_m2_batched(m1_inclusion_matrix, group_ids, batch_size,
                               unique_groups, n_threads, verbose)
  } else {
    if (verbose) {
      if (force_fast) {
        message("|   |-- Using fast processing approach (forced by user)")
      } else {
        message("|   |-- Memory requirements within limits")
        message("|   +-- Using fast processing approach")
      }
    }

    # Call the fast processing function
    result <- .make_m2_fast(m1_inclusion_matrix, group_ids, verbose)
  }

  if (verbose) message("Finished M2 matrix creation.")
  return(result)
}

#' Fast M2 Processing (Internal Function)
#'
#' Implements the original fast approach using data.table operations.
#' This approach creates the full operation in memory at once.
#'
#' @param m1_inclusion_matrix The M1 inclusion matrix.
#' @param group_ids Data table with group IDs.
#' @param verbose Logical for verbose output.
#'
#' @return A sparse M2 matrix.
#' @noRd
.make_m2_fast <- function(m1_inclusion_matrix, group_ids, verbose) {

  if (verbose) message("+-- Step 3 | Creating M2 (fast approach)")

  tryCatch({
    # Convert m1_inclusion_matrix to data.table
    m1 <- summary(m1_inclusion_matrix) |> data.table::as.data.table()
    data.table::setnames(m1, "x", "x_1")

    if (verbose) message("|   |-- Converted matrix to data.table format")

    # Merge group information
    m1 <- merge(m1, group_ids, by = "i")
    m1[, x_tot := sum(x_1), .(group_id, j)]
    m_tot <- m1[, .(group_id, j, x_tot)] |> unique()

    if (verbose) message("|   |-- Calculated group totals")

    # Filter and merge relevant data
    m_tot <- m_tot[x_tot > 0]
    m_tot <- merge(m_tot, group_ids, by = "group_id", allow.cartesian = TRUE)
    m_tot <- merge(m_tot, m1, by = c("group_id", "i", "j", "x_tot"), all.x = TRUE)
    m_tot[is.na(x_1), x_1 := 0]
    m_tot[, x_2 := x_tot - x_1]

    if (verbose) message("|   |-- Created final calculation table")

    # Create sparse matrix for M2_train
    M2_train <- m_tot[, Matrix::sparseMatrix(i = i, j = j, x = x_2)]

    if (verbose) message("|   +-- Generated M2 sparse matrix")

  }, error = function(e) {
    if (grepl("2\\^31", e$message) || grepl("vecseq", e$message)) {
      stop("Memory limit exceeded in fast approach. The dataset is too large for the fast method. ",
           "Please try reducing the data size or contact support.", call. = FALSE)
    } else {
      stop("Error in fast M2 processing: ", e$message, call. = FALSE)
    }
  })

  if (verbose) message("+-- Step 4 | Finalizing M2 creation")

  # Set row and column names
  rownames(M2_train) <- rownames(m1_inclusion_matrix)
  colnames(M2_train) <- colnames(m1_inclusion_matrix)

  # Remove dummy row from M2_train
  M2 <- M2_train[-nrow(M2_train), ]

  if (verbose) message("+-- All done!")
  return(M2)
}

#' Batched M2 Processing (Internal Function)
#'
#' Implements the batched triplet combination approach for memory-efficient processing.
#' Processes groups in batches using lapply/mclapply and combines results efficiently.
#'
#' @param m1_inclusion_matrix The M1 inclusion matrix.
#' @param group_ids Data table with group IDs.
#' @param batch_size Number of groups per batch.
#' @param unique_groups Vector of unique group identifiers.
#' @param n_threads Number of threads for parallel processing.
#' @param verbose Logical for verbose output.
#'
#' @return A sparse M2 matrix.
#' @noRd
.make_m2_batched <- function(m1_inclusion_matrix, group_ids, batch_size,
                             unique_groups, n_threads, verbose) {

  if (verbose) message("+-- Step 3 | Creating M2 (batched processing approach)")

  # Prepare batch processing
  n_groups <- length(unique_groups)
  n_batches <- ceiling(n_groups / batch_size)

  if (verbose) message("|   |-- Processing ", n_groups, " groups in ", n_batches, " batches")

  if (n_threads > 1) {
    if (verbose) message("|   |-- Using parallel processing with mclapply (", n_threads, " threads)")
  } else {
    if (verbose) message("|   |-- Using sequential processing with lapply")
  }

  # Store target dimensions for final sparse matrix creation
  target_nrow <- nrow(m1_inclusion_matrix)
  target_ncol <- ncol(m1_inclusion_matrix)

  # Convert m1_inclusion_matrix to data.table for merging (once, outside the apply)
  tryCatch({
    m1 <- summary(m1_inclusion_matrix) |> data.table::as.data.table()
    data.table::setnames(m1, "x", "x_1")
    m1 <- merge(m1, group_ids, by = "i")

    # Pre-calculate group totals for all groups
    m1[, x_tot := sum(x_1), .(group_id, j)]
    group_cols <- unique(m1[x_tot > 0, .(group_id, j, x_tot)])

    if (verbose) message("|   |-- Pre-calculated group totals for all groups")

  }, error = function(e) {
    stop("Error in initial data preparation for batched processing: ", e$message, call. = FALSE)
  })

  # Create batch processing function
  process_batch <- function(batch_idx) {
    tryCatch({
      start_idx <- (batch_idx - 1) * batch_size + 1
      end_idx <- min(batch_idx * batch_size, n_groups)
      current_groups <- unique_groups[start_idx:end_idx]

      if (verbose) {
        message("|   |---- Processing batch ", batch_idx, " of ", n_batches,
                " (groups ", start_idx, " to ", end_idx, ")")
      }

      # Get subset for current batch
      batch_group_cols <- group_cols[group_id %in% current_groups]
      batch_group_ids <- group_ids[group_id %in% current_groups]

      # Create Cartesian product for this batch only
      batch_grid <- merge(batch_group_cols, batch_group_ids,
                          by = "group_id", allow.cartesian = TRUE)

      # Merge with original data
      batch_result <- merge(batch_grid,
                            m1[group_id %in% current_groups, .(group_id, i, j, x_1)],
                            by = c("group_id", "i", "j"), all.x = TRUE)

      # Calculate x_2
      batch_result[is.na(x_1), x_1 := 0]
      batch_result[, x_2 := x_tot - x_1]

      # Return only the triplet data (i, j, x_2)
      return(batch_result[, .(i, j, x_2)])

    }, error = function(e) {
      stop("Error processing batch ", batch_idx, ": ", e$message, call. = FALSE)
    })
  }

  # Process all batches using lapply or mclapply
  batch_triplets <- tryCatch({
    if (n_threads > 1) {
      parallel::mclapply(1:n_batches, process_batch, mc.cores = n_threads)
    } else {
      lapply(1:n_batches, process_batch)
    }
  }, error = function(e) {
    stop("Error in batch processing: ", e$message, call. = FALSE)
  })

  # Combine all triplets and create final sparse matrix
  if (verbose) message("|   +---- Combining all triplets and creating final sparse matrix")

  tryCatch({
    # Combine all triplet data.tables
    all_triplets <- data.table::rbindlist(batch_triplets)

    # Clean up batch triplets to free memory
    rm(batch_triplets)
    gc()

    # Create the final sparse matrix from combined triplets
    M2_train <- Matrix::sparseMatrix(
      i = all_triplets$i,
      j = all_triplets$j,
      x = all_triplets$x_2,
      dims = c(target_nrow, target_ncol)
    )

    # Clean up triplets
    rm(all_triplets)
    gc()

  }, error = function(e) {
    stop("Error combining triplets and creating final matrix: ", e$message, call. = FALSE)
  })

  if (verbose) message("+-- Step 4 | Finalizing M2 creation")

  # Set row and column names
  rownames(M2_train) <- rownames(m1_inclusion_matrix)
  colnames(M2_train) <- colnames(m1_inclusion_matrix)

  # Remove the dummy row
  M2 <- M2_train[-nrow(M2_train), ]

  if (verbose) message("+-- All done!")
  return(M2)
}

#' make_eventdata_plus
#' @description
#' This function reads in a GTF file, extracts gene annotations, and merges them with
#' event-level genomic intervals provided by the user. The final data table contains the
#' original event intervals and the corresponding gene information (for example, `gene_id`
#' and `gene_name`).
#'
#' @details
#' 1. Read the GTF: Uses `data.table::fread` to load GTF data and convert it to
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
#' @export
make_eventdata_plus <- function(eventdata, GTF_file_direction) {

  # Read GTF file as plain text using fread
  GTF <- data.table::fread(
    GTF_file_direction,
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attribute"),
    sep = "\t",
    header = FALSE,
    quote = "",
    showProgress = FALSE
  )

  # Filter for 'gene' entries
  ref_gtf <- GTF[type == "gene"]

  # Extract gene_id and gene_name from the attribute column
  ref_gtf[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attribute)]
  ref_gtf[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attribute)]

  # Select and rename relevant columns
  ref_gtf <- ref_gtf[, .(chr = seqid, start, end, strand, gene_id, gene_name)]

  # Add 'chr' prefix if missing
  ref_gtf[!grepl("^chr", chr), chr := paste0("chr", chr)]
  ref_gtf[, chr := as.character(chr)]

  # Convert strand from + / - to 1 / 2
  strand_map <- data.table(strand = c("+", "-"), new_strand = c(1, 2))
  ref_gtf <- merge(ref_gtf, strand_map, by = "strand", sort = FALSE)
  ref_gtf[, strand := NULL]
  setnames(ref_gtf, "new_strand", "strand")

  # Prepare eventdata
  eventdata[, chr := as.character(chr)]
  eventdata[!grepl("^chr", chr), chr := paste0("chr", chr)]

  # Set keys for foverlaps
  data.table::setkey(eventdata, chr, strand, start, end)
  data.table::setkey(ref_gtf, chr, strand, start, end)

  # Overlap join
  new_eventdata <- data.table::foverlaps(eventdata, ref_gtf, type = "within")
  new_eventdata <- na.omit(new_eventdata)

  return(new_eventdata)
}

#' make_gene_count
#'
#' @description
#' Constructs sparse gene expression matrices from one or more directories containing 10X Genomics-style output.
#' The function supports barcode filtering using either an external whitelist or the internally provided filtered barcode file.
#'
#' @param expression_dirs A character vector or list of strings. Each element must be a path to a directory containing the
#'   gene expression matrix files: \code{matrix.mtx}, \code{barcodes.tsv}, and \code{features.tsv} (or \code{genes.tsv}).
#'
#' @param sample_ids A character vector or list of unique sample identifiers, one for each element in \code{expression_dirs}.
#'   These are used to name outputs in the returned list when multiple samples are provided.
#'
#' @param whitelist_barcodes A list of character vectors. Each list element corresponds to a sample and contains the
#'   barcodes to retain for that sample. If \code{NULL} (default), the function will attempt to use the internal filtered
#'   barcode file (e.g., \code{barcodes.tsv} or \code{barcodes_filtered.tsv}) if available.
#'
#' @param use_internal_whitelist Logical (default \code{TRUE}). If \code{TRUE} and \code{whitelist_barcodes} is \code{NULL},
#'   the function will attempt to use the default filtered barcode list from the input directory.
#'   If \code{FALSE}, no internal filtration will be applied unless a whitelist is explicitly provided.
#' @param verbose Logical. If \code{TRUE}, prints progress and informational messages. Default is \code{FALSE}.
#'
#' @return
#' If a single sample is provided, returns a sparse matrix of class \code{"dgCMatrix"} with genes as rows and barcodes as columns.
#' If multiple samples are provided, returns a named list of sparse matrices, one per sample ID.
#'
#' @details
#' The function is designed for bulk or single-cell gene expression processing from 10X-style output folders.
#' Each input directory should contain the standard \code{matrix.mtx}, \code{features.tsv}/\code{genes.tsv}, and \code{barcodes.tsv}
#' files. Barcodes can be filtered using either a provided whitelist or by relying on the \code{filtered} barcode files
#' output by tools like CellRanger.
#'
#' If neither an external whitelist nor an internal filtered barcode file is available, all barcodes from the
#' raw matrix will be retained.
#'
#' @section Dependencies:
#' Requires the \pkg{Matrix} package for sparse matrix handling and potentially \pkg{data.table} for efficient I/O.
#'
#' @export
make_gene_count <- function(expression_dirs, sample_ids, whitelist_barcodes = NULL, use_internal_whitelist = TRUE, verbose = FALSE) {

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
    if (verbose) message("|-- Processing gene expression data for sample: ", sample_id)
    g_mtx <- Matrix::readMM(expression_matrix_dir)
    g_brc <- data.table::fread(expression_barcodes_dir, header = FALSE, showProgress = FALSE)$V1
    g_feature <- data.table::fread(expression_features_dir, header = FALSE, showProgress = FALSE)

    # Set row and column names
    rownames(g_mtx) <- g_feature$V1
    colnames(g_mtx) <- g_brc

    # Apply barcode filtration
    if (!is.null(whitelist_barcode)) {
      if (verbose) message("|  |--  Applying provided whitelist for sample: ", sample_id)
    } else if (use_internal_whitelist && file.exists(filtered_barcodes_dir)) {
      whitelist_barcode <- data.table::fread(filtered_barcodes_dir, header = FALSE, showProgress = FALSE)$V1
      if (verbose) message("|  |--  Using filtered barcodes for sample: ", sample_id)
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

      if (verbose) message("|  |--  Filtered barcodes applied for sample: ", sample_id,
                           " (Remaining barcodes: ", final_barcodes, ")")
    } else {
      if (verbose) message("|  |--  No barcode filtration applied for sample: ", sample_id)
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
    if (verbose) message("|  +-- Finished processing gene expression data for sample: ", sample_id)
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
  if (verbose) {
    message("\nSummary of Processed Samples:")
    print(summary_table)
  }

  # Return results
  if (length(final_results) == 1) {
    return(final_results[[1]])  # Return the single result as a sparse matrix
  } else {
    names(final_results) <- sample_ids
    return(final_results)
  }
}

#' Process Spliced and Unspliced Counts from Velocyto Outputs
#'
#' @description
#' Parses and processes spliced and unspliced gene expression matrices from one or more Velocyto output directories.
#' The function applies barcode filtering using an external whitelist or filtered barcodes file, and optionally merges
#' the results across samples into unified matrices.
#'
#' @param velocyto_dirs A character vector or list of strings. Each element should be a path to a Velocyto output directory.
#'   Each directory must contain subdirectories (typically \code{filtered} or \code{raw}) with the required matrix files:
#'   \code{spliced.mtx}, \code{unspliced.mtx}, \code{barcodes.tsv}, and \code{genes.tsv} or \code{features.tsv}.
#'
#' @param sample_ids A character vector or list of unique sample identifiers corresponding to each entry in \code{velocyto_dirs}.
#'
#' @param whitelist_barcodes A list of character vectors. Each element should provide a whitelist of barcodes to retain for the
#'   corresponding sample. If \code{NULL} (default), the function will attempt to use the internally provided filtered barcodes
#'   when \code{use_internal_whitelist = TRUE}.
#'
#' @param use_internal_whitelist Logical (default \code{TRUE}). If \code{TRUE}, and \code{whitelist_barcodes} is \code{NULL},
#'   the function uses the filtered barcode file (if present in the directory). If \code{FALSE}, all barcodes from the
#'   raw matrix will be used unless a whitelist is explicitly provided.
#'
#' @param merge_counts Logical (default \code{FALSE}). If \code{TRUE}, spliced and unspliced matrices across all samples are merged
#'   into two combined matrices (one for spliced, one for unspliced). If \code{FALSE}, the results are returned per sample.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress and informational messages. Default is \code{FALSE}.
#'
#' @return
#' A list containing processed gene expression matrices:
#' \itemize{
#'   \item If \code{merge_counts = FALSE}, returns a named list of sample-specific matrices. Each entry contains:
#'     \describe{
#'       \item{\code{spliced}}{Sparse matrix of spliced transcript counts.}
#'       \item{\code{unspliced}}{Sparse matrix of unspliced transcript counts.}
#'     }
#'   \item If \code{merge_counts = TRUE}, returns a list with two elements:
#'     \describe{
#'       \item{\code{spliced}}{Merged sparse matrix of spliced counts across all samples.}
#'       \item{\code{unspliced}}{Merged sparse matrix of unspliced counts across all samples.}
#'     }
#' }
#'
#' @details
#' The function assumes that each Velocyto directory follows the 10X-like structure typically produced by tools like
#' Loom or Velocyto CLI. Barcode filtering ensures that only high-quality or selected barcodes are retained
#' for downstream RNA velocity analysis.
#'
#' When merging matrices, barcodes are prefixed with their corresponding sample ID to avoid collisions and preserve traceability.
#'
#' @section Dependencies:
#' Requires the \pkg{Matrix} package for sparse matrix operations and \pkg{data.table} for efficient file parsing.
#'
#' @export

make_velo_count <- function(velocyto_dirs, sample_ids, whitelist_barcodes = NULL, use_internal_whitelist = TRUE, merge_counts = FALSE, verbose = FALSE) {

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
    if (verbose) message("|-- Processing Velocyto data for sample: ", sample_id)
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
      if (verbose) message("|  |--  Applying provided whitelist for sample: ", sample_id)
    } else if (use_internal_whitelist && file.exists(filtered_barcodes_dir)) {
      whitelist_barcode <- data.table::fread(filtered_barcodes_dir, header = FALSE, showProgress = FALSE)$V1
      if (verbose) message("|  |--  Using filtered barcodes for sample: ", sample_id)
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

      if (verbose) message("|  |--  Filtered barcodes applied for sample: ", sample_id,
          " (Remaining barcodes: ", final_barcodes, ")")
    } else {
      if (verbose) message("|  |--  No barcode filtration applied for sample: ", sample_id)
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
    if (verbose) message("|  +-- Finished processing Velocyto data for sample: ", sample_id)
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

    if (verbose) {
      message("\nSummary of Merged Velocyto Counts:")
      print(summary_table)
    }

    return(list(
      spliced = spliced_combined,
      unspliced = unspliced_combined
    ))
  }

  # Print summary table if not merging
  if (verbose) {
    message("\nSummary of Processed Samples:")
    print(summary_table)
  }

  # Return individual results if not merging
  final_results <- lapply(results, function(x) list(spliced = x$spliced, unspliced = x$unspliced))
  if (length(final_results) == 1) {
    return(final_results[[1]])  # Return the single result as a list of spliced/unspliced
  } else {
    names(final_results) <- sample_ids
    return(final_results)  # Return results for multiple samples
  }
}

#' Load the toy M1/M2 object
#'
#' Loads a toy object of M1 and M2 used for examples or testing.
#'
#' @return An R object from the toy_m1_m2_obj.rds file.
#' @export
load_toy_M1_M2_object <- function() {
  file <- system.file("extdata", "toy_m1_m2_obj.rds", package = "splikit")
  readRDS(file)
}
