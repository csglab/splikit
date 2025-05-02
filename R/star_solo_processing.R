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

make_junction_ab <- function(STARsolo_SJ_dirs, white_barcode_lists = NULL, sample_ids, use_internal_whitelist = TRUE) {

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

make_m1 <- function(junction_ab_object) {
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
  all_junctions <- unique(all_in_one_eventdata$row_names_mtx)
  all_in_one_eventdata[, sample_id := NULL]
  temp_eventdata <- all_in_one_eventdata[match(all_junctions, row_names_mtx), ]

  # grouping based on first and last coordinates using data.table
  if (!is.data.table(temp_eventdata)) {
    setDT(temp_eventdata)
  }

  temp_eventdata[, `:=`(
    start_cor_group_id = .GRP,
    start_cor_group_count = .N
  ), by = start_cor_id]

  temp_eventdata[, `:=`(
    end_cor_group_id = .GRP,
    end_cor_group_count = .N
  ), by = end_cor_id]

  temp_eventdata_grouped <- temp_eventdata

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
  } else {
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

#' make_m2
#'
#' Creates the M2 matrix from a given m1_inclusion_matrix and eventdata, ensuring the proper processing of group indices and matrix operations.
#'
#' @param m1_inclusion_matrix A sparse matrix to be modified and used for creating the M2 matrix.
#' @param eventdata A data.table containing event information with at least `group_id` and an index column.
#' @return A sparse matrix M2 with the dummy row removed and proper adjustments made.
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
make_m2 <- function(m1_inclusion_matrix, eventdata) {
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

  cat("┌── Step 1 | Modifying the m1_inclusion_matrix\n")

  # Add dummy group to group_ids
  dummy_group <- data.table(i = nrow(m1_inclusion_matrix), group_id = 'dummy')
  group_ids <- eventdata[, .(i, group_id)]
  group_ids <- rbind(group_ids, dummy_group)

  rm(dummy_group)  # Remove intermediate variable

  # Add group index and initialize variables
  num_cells <- ncol(m1_inclusion_matrix)
  groups_start_vector <- eventdata[, unique(group_id)]

  cat("├── Step 2 | Creating M2\n")

  # Convert m1_inclusion_matrix to data.table
  m1 <- summary(m1_inclusion_matrix) |> as.data.table()
  setnames(m1, 'x', 'x_1')

  # Merge group information
  m1 <- merge(m1, group_ids, by = 'i')
  m1[, x_tot := sum(x_1), .(group_id, j)]
  m_tot <- m1[, .(group_id, j, x_tot)] |> unique()

  # Filter and merge relevant data
  m_tot <- m_tot[x_tot > 0]
  m_tot <- merge(m_tot, group_ids, by = 'group_id', allow.cartesian = TRUE)
  m_tot <- merge(m_tot, m1, by = c('group_id', 'i', 'j', 'x_tot'), all.x = TRUE)
  m_tot[is.na(x_1), x_1 := 0]
  m_tot[, x_2 := x_tot - x_1]

  # Create sparse matrix for M2_train
  M2_train <- m_tot[, sparseMatrix(i = i, j = j, x = x_2)]

  cat("├── Step 3 | Finalizing M2 creation\n")

  # Set row and column names
  rownames(M2_train) <- rownames(m1_inclusion_matrix)
  colnames(M2_train) <- colnames(m1_inclusion_matrix)

  # Remove dummy row from M2_train
  M2 <- M2_train[-nrow(M2_train), ]

  cat("└── All done!")
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
make_gene_count <- function(expression_dirs, sample_ids, whitelist_barcodes = NULL, use_internal_whitelist = TRUE) {

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

make_velo_count <- function(velocyto_dirs, sample_ids, whitelist_barcodes = NULL, use_internal_whitelist = TRUE, merge_counts = FALSE) {

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

