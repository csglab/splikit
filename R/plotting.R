# plotting.R --------------------------------------------------------------
# Transcript-exclusive splice-junction plots. Works with both Ensembl-style
# (gene_name/transcript_name) and RefSeq-style (gene_id/transcript_id) GTFs.

# --- internal helpers ---------------------------------------------------------

# Extract a value from a GTF attribute column by trying each key in turn.
# Returns NA for elements where none of the keys match. Designed so a caller
# can write  .gtf_attr(attr, "gene_name", "gene_id", "gene")  and get the
# first non-NA hit per row.
.gtf_attr <- function(x, ...) {
  keys <- c(...)
  out  <- rep(NA_character_, length(x))
  for (k in keys) {
    if (!any(is.na(out))) break
    re   <- sprintf('.*%s "([^"]+)".*', k)
    todo <- is.na(out)
    hit  <- todo & grepl(re, x)
    if (any(hit))
      out[hit] <- sub(re, "\\1", x[hit])
  }
  out
}

# Load only the GTF rows relevant to `target_gene`. Accepts a path or an
# already-loaded data.table; normalises column names to (seqname, ..., attr).
# On paths, tries fast grep prefilters against `gene_name` and `gene_id`
# (Ensembl vs RefSeq) before falling back to a full fread.
.load_gene_gtf <- function(gtf, target_gene) {
  gtf_cols <- c("seqname", "source", "type", "start", "end",
                "score", "strand", "frame", "attr")

  if (data.table::is.data.table(gtf)) {
    g <- data.table::copy(gtf)
    if (!"seqname" %in% names(g) && "seqid"     %in% names(g))
      data.table::setnames(g, "seqid", "seqname")
    if (!"attr"    %in% names(g) && "attribute" %in% names(g))
      data.table::setnames(g, "attribute", "attr")
    return(g)
  }

  path <- path.expand(gtf)

  for (key in c("gene_name", "gene_id")) {
    g <- tryCatch(
      data.table::fread(
        cmd = sprintf("grep -F '%s \"%s\"' %s",
                      key, target_gene, shQuote(path)),
        sep = "\t", header = FALSE, col.names = gtf_cols, quote = "",
        showProgress = FALSE
      ),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(g) && nrow(g)) return(g)
  }

  # last resort: full read
  data.table::fread(path, sep = "\t", header = FALSE,
                    col.names = gtf_cols, quote = "",
                    skip = "#", showProgress = FALSE)
}

# From a raw gene-scope GTF slice, derive the exon table used downstream.
# Returns a data.table with seqname/strand/transcript_id/transcript_name/
# exon_number/start/end. Handles RefSeq GTFs (no gene_name/transcript_name)
# by falling back to gene_id/transcript_id.
.gene_exons <- function(g, target_gene) {
  g[, gene_name := .gtf_attr(attr, "gene_name", "gene_id", "gene")]
  g <- g[gene_name == target_gene]

  if (!nrow(g))
    stop(sprintf("No GTF rows matched gene '%s'.", target_gene), call. = FALSE)

  g[, gene_id := .gtf_attr(attr, "gene_id")]

  exons <- g[type == "exon"]
  exons[, transcript_id   := .gtf_attr(attr, "transcript_id")]
  exons[, transcript_name := .gtf_attr(attr, "transcript_name",
                                             "transcript_id")]
  exons[, exon_number     := as.integer(.gtf_attr(attr, "exon_number"))]

  keep <- c("seqname", "strand", "transcript_id", "transcript_name",
            "exon_number", "start", "end")
  exons <- exons[, ..keep]
  if (!nrow(exons))
    stop(sprintf("No exon rows for gene '%s'.", target_gene), call. = FALSE)

  attr(exons, "gene_id") <- unique(stats::na.omit(g$gene_id))[1]
  exons
}

# Per-transcript GTF introns + gene-wide exclusivity.
.tx_junctions <- function(exons) {
  data.table::setorder(exons, transcript_id, start)
  j <- exons[, {
    if (.N >= 2) .(j_start = end[-.N] + 1L, j_end = start[-1] - 1L)
    else         .(j_start = integer(0),    j_end = integer(0))
  }, by = .(transcript_id, transcript_name, strand, seqname)]

  counts <- j[, .(n_tx = data.table::uniqueN(transcript_id)),
              by = .(j_start, j_end)]
  j <- merge(j, counts, by = c("j_start", "j_end"))
  j[, exclusive := n_tx == 1L]
  j
}

# Assemble the y-axis map (which transcripts are drawn, in what order).
# `override_tx` (when non-NULL) pins the plot to those transcripts;
# `show_exclusive` restricts to those owning >=1 exclusive junction;
# otherwise all transcripts of the gene are drawn.
.order_transcripts <- function(junctions, all_tx, show_exclusive, override_tx) {
  tx_order <- junctions[, .(has_excl = any(exclusive)), by = transcript_name
                        ][order(-has_excl, transcript_name)]

  if (!is.null(override_tx)) {
    missing <- setdiff(override_tx, all_tx)
    if (length(missing))
      stop("Transcript(s) not in gene: ",
           paste(missing, collapse = ", "), call. = FALSE)
    extra <- setdiff(override_tx, tx_order$transcript_name)
    if (length(extra))
      tx_order <- rbind(tx_order,
                        data.table::data.table(transcript_name = extra,
                                               has_excl = FALSE))
    tx_order <- tx_order[transcript_name %in% override_tx]
    tx_order <- tx_order[match(override_tx, transcript_name)]
  } else if (isTRUE(show_exclusive)) {
    excl_subset <- tx_order[has_excl == TRUE]
    if (nrow(excl_subset)) {
      tx_order <- excl_subset
    } else {
      message("No transcript owns an exclusive junction under the current ",
              "filter; showing all transcripts instead.")
    }
  } else {
    tx_order <- rbind(
      tx_order,
      data.table::data.table(
        transcript_name = setdiff(all_tx, tx_order$transcript_name),
        has_excl = FALSE))
  }

  tx_order
}

# Build the ggplot object.
.build_plot <- function(exons, junctions, tx_map, target_gene,
                        curvature, title_suffix) {
  intron_spans <- exons[, .(x = min(start), xend = max(end)),
                        by = .(transcript_name, y)]
  chrom <- unique(exons$seqname)[1]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = intron_spans,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          color = "grey55", linewidth = 0.35) +
    ggplot2::geom_rect(data = exons,
                       ggplot2::aes(xmin = start, xmax = end,
                                    ymin = y - 0.28, ymax = y + 0.28),
                       fill = "grey85", color = "black", linewidth = 0.25)

  shared <- junctions[exclusive == FALSE]
  if (nrow(shared))
    p <- p + ggplot2::geom_curve(
      data = shared,
      ggplot2::aes(x = j_start, xend = j_end, y = y + 0.28, yend = y + 0.28),
      curvature = curvature, linewidth = 0.25, color = "grey70")

  solo <- junctions[exclusive == TRUE]
  if (nrow(solo))
    p <- p + ggplot2::geom_curve(
      data = solo,
      ggplot2::aes(x = j_start, xend = j_end, y = y + 0.28, yend = y + 0.28),
      curvature = curvature, linewidth = 0.7, color = "black")

  p +
    ggplot2::scale_y_continuous(breaks = tx_map$y, labels = tx_map$y_label,
                                expand = ggplot2::expansion(add = 0.9)) +
    ggplot2::labs(
      x = sprintf("%s position (bp)", chrom), y = NULL,
      title = sprintf("%s - transcript-exclusive junctions%s",
                      target_gene, title_suffix),
      caption = paste0(
        "Black arcs: junctions exclusive to a single transcript. ",
        "Grey arcs: shared junctions.")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 8, lineheight = 0.9),
      plot.caption       = ggplot2::element_text(hjust = 0, size = 8,
                                                 color = "grey30")
    )
}

# Build the `$info` data.table returned to the caller. One row per
# (transcript, junction) drawn. `eventdata` (optional) supplies per-junction
# annotations (row_names_mtx, is_annot) for the event-based function.
.build_info <- function(junctions, exons, target_gene, gene_id_val,
                        eventdata = NULL) {
  info <- junctions[, .(
    gene_name             = target_gene,
    gene_id               = gene_id_val,
    transcript_name, transcript_id,
    chr                   = seqname,
    strand,
    j_start, j_end,
    j_width               = j_end - j_start + 1L,
    exclusive,
    n_tx_with_junction    = n_tx
  )]

  if (!is.null(eventdata) && all(c("i.start","i.end") %in% names(eventdata))) {
    ed <- eventdata[, .(i.start, i.end,
                        row_names_mtx = if ("row_names_mtx" %in% names(eventdata))
                          row_names_mtx else NA_character_,
                        is_annot      = if ("is_annot" %in% names(eventdata))
                          is_annot else NA_integer_)]
    ed_u <- ed[, .(row_names_mtx = row_names_mtx[1],
                   is_annot      = any(is_annot == 1L, na.rm = TRUE)),
               by = .(i.start, i.end)]
    info <- merge(info, ed_u,
                  by.x = c("j_start","j_end"),
                  by.y = c("i.start","i.end"),
                  all.x = TRUE, sort = FALSE)
    info[, observed_in_eventdata := TRUE]
  } else {
    info[, observed_in_eventdata := NA]
    info[, row_names_mtx         := NA_character_]
    info[, is_annot              := NA]
  }

  data.table::setcolorder(info, c(
    "gene_name", "gene_id", "transcript_name", "transcript_id",
    "chr", "strand", "j_start", "j_end", "j_width",
    "exclusive", "n_tx_with_junction",
    "observed_in_eventdata", "row_names_mtx", "is_annot"
  ))
  data.table::setorder(info, transcript_name, j_start)
  info
}

.ensure_ggplot <- function(fn_name) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(sprintf("Package 'ggplot2' is required for %s(). ", fn_name),
         "Install it with install.packages('ggplot2').", call. = FALSE)
}


# --- exported functions ------------------------------------------------------

#' Plot transcript-exclusive splice junctions
#'
#' @description
#' Draws a gene-model view of a gene's transcripts (one row per transcript)
#' with splice junctions rendered as arcs above each row. Junctions used by a
#' single transcript of the gene ("transcript-exclusive") are drawn as solid
#' black arcs; shared junctions as thin grey arcs. Works with both
#' Ensembl-style and RefSeq-style GTFs: if a row has no `gene_name`
#' attribute, `gene_id` is used; if no `transcript_name`, `transcript_id`
#' is used.
#'
#' @return An S3 object of class `"splikit_junction_plot"` (a list) with
#'   components:
#'   \describe{
#'     \item{`plot`}{A `ggplot` object.}
#'     \item{`info`}{A `data.table` with one row per drawn (transcript,
#'           junction) pair and columns \code{gene_name}, \code{gene_id},
#'           \code{transcript_name}, \code{transcript_id}, \code{chr},
#'           \code{strand}, \code{j_start}, \code{j_end}, \code{j_width},
#'           \code{exclusive}, \code{n_tx_with_junction},
#'           \code{observed_in_eventdata}, \code{row_names_mtx},
#'           \code{is_annot}. For GTF-only calls the last three are `NA`.}
#'     \item{`exons`, `junctions`, `tx_order`}{The underlying tables used to
#'           build the plot, retained for advanced users.}
#'   }
#'   Printing the object renders the plot and then prints `info`.
#'
#' @param gtf Either a path to an Ensembl- or RefSeq-style GTF, or a
#'   pre-loaded `data.table` of GTF rows.
#' @param target_gene Gene symbol (matches `gene_name`, or `gene_id` if no
#'   `gene_name` is present).
#' @param show_exclusive Logical (default `TRUE`). Restrict to transcripts
#'   that own >= 1 exclusive junction.
#' @param transcript Optional character vector of transcript names to pin
#'   the plot to (overrides `show_exclusive`).
#' @param curvature Numeric, arc-height knob for `ggplot2::geom_curve()`.
#' @param out_file Optional character. If given, the plot is also written
#'   to this file with `ggplot2::ggsave()`.
#'
#' @section Required packages:
#' Requires the \pkg{ggplot2} package (declared in `Suggests`).
#'
#' @examples
#' \dontrun{
#'   plot_exclusive_junctions("Mus_musculus.GRCm39.110.gtf", "Tpm1")
#' }
#'
#' @export
plot_exclusive_junctions <- function(gtf,
                                     target_gene,
                                     show_exclusive = TRUE,
                                     transcript     = NULL,
                                     curvature      = -0.2,
                                     out_file       = NULL) {

  .ensure_ggplot("plot_exclusive_junctions")

  g        <- .load_gene_gtf(gtf, target_gene)
  exons    <- .gene_exons(g, target_gene)
  gene_id_val <- attr(exons, "gene_id")
  junctions <- .tx_junctions(exons)

  all_tx   <- unique(exons$transcript_name)
  tx_order <- .order_transcripts(junctions, all_tx,
                                 show_exclusive = show_exclusive,
                                 override_tx    = transcript)

  tx_id_map <- unique(exons[, .(transcript_name, transcript_id)])
  tx_map    <- merge(data.table::data.table(transcript_name = tx_order$transcript_name,
                                            y = seq_len(nrow(tx_order))),
                     tx_id_map, by = "transcript_name", sort = FALSE)
  data.table::setorder(tx_map, y)
  tx_map[, y_label := ifelse(transcript_name == transcript_id,
                             transcript_name,
                             paste0(transcript_name, "\n", transcript_id))]

  if (!nrow(tx_map))
    stop(sprintf("No transcripts to plot for '%s' with the given filters.",
                 target_gene), call. = FALSE)

  exons     <- merge(exons,     tx_map[, .(transcript_name, y)], by = "transcript_name")
  junctions <- merge(junctions, tx_map[, .(transcript_name, y)], by = "transcript_name")

  title_suffix <- if (!is.null(transcript))
                    sprintf(" (%s)", paste(transcript, collapse = ", "))
                  else if (isTRUE(show_exclusive))
                    " (exclusive transcripts only)"
                  else ""

  p <- .build_plot(exons, junctions, tx_map, target_gene,
                   curvature, title_suffix)

  if (!is.null(out_file))
    ggplot2::ggsave(out_file, p,
                    width  = 10,
                    height = max(3, 0.45 * nrow(tx_map) + 2),
                    dpi    = 150)

  info <- .build_info(junctions, exons, target_gene, gene_id_val)

  structure(
    list(plot = p, info = info,
         exons = exons, junctions = junctions, tx_order = tx_map),
    class = c("splikit_junction_plot", "list")
  )
}

#' Multi-page PDF of transcript-exclusive junctions
#'
#' @description
#' Writes a multi-page PDF: page 1 = full gene view (all transcripts,
#' exclusive arcs in black); subsequent pages = one per exclusive-owning
#' transcript, with its exclusive junction in black.
#'
#' @param gtf,target_gene,curvature See [plot_exclusive_junctions()].
#' @param out_pdf Character path for the PDF.
#' @param width,height PDF page dimensions. `height = NULL` (default) sizes
#'   to the full-gene page.
#'
#' @return Invisibly, a list with `exclusive_transcripts` and `n_pages`.
#'
#' @section Required packages:
#' Requires the \pkg{ggplot2} package (declared in `Suggests`).
#'
#' @export
plot_exclusive_junctions_pdf <- function(gtf, target_gene, out_pdf,
                                         curvature = -0.2,
                                         width = 10, height = NULL) {
  .ensure_ggplot("plot_exclusive_junctions_pdf")

  full <- plot_exclusive_junctions(gtf, target_gene,
                                   show_exclusive = FALSE,
                                   curvature      = curvature)
  excl_tx <- unique(full$info[exclusive == TRUE, transcript_name])

  if (is.null(height))
    height <- max(5, 0.28 * nrow(full$tx_order) + 2)

  grDevices::pdf(out_pdf, width = width, height = height)
  on.exit(grDevices::dev.off())

  print(full$plot)
  for (tx in excl_tx) {
    r <- plot_exclusive_junctions(gtf, target_gene,
                                  transcript = tx,
                                  curvature  = curvature)
    print(r$plot)
  }
  invisible(list(exclusive_transcripts = excl_tx,
                 n_pages = length(excl_tx) + 1L))
}

#' Plot transcript-exclusive splice junctions observed in eventdata
#'
#' @description
#' Sibling of [plot_exclusive_junctions()] that restricts drawn arcs to
#' junctions present in a `splikit` eventdata table. Exon structure and
#' exclusivity are still derived from the GTF, so a black arc means the
#' junction is both transcript-exclusive in the annotation and observed in
#' the data.
#'
#' @details
#' Works with both Ensembl-style and RefSeq-style GTFs. If `eventdata`
#' lacks a `gene_name` column the function runs [make_eventdata_plus()]
#' internally (requires `GTF` to be a path in that case). Observed
#' junctions not matching any annotated intron are dropped with a message
#' and their count is reported as `novel_junction_count`.
#'
#' @param target_gene Gene symbol to plot.
#' @param GTF Either a GTF path or a pre-loaded `data.table`.
#' @param eventdata A `splikit` eventdata `data.table`.
#' @param show_exclusive,transcript,curvature See [plot_exclusive_junctions()].
#'
#' @return An S3 object of class `"splikit_junction_plot"` with the same
#'   structure as in [plot_exclusive_junctions()], plus
#'   `novel_junction_count` for the unannotated observed junctions. The
#'   `info$observed_in_eventdata` column is `TRUE` for all rows, and
#'   `row_names_mtx` / `is_annot` are populated from `eventdata` when
#'   available.
#'
#' @section Required packages:
#' Requires the \pkg{ggplot2} package (declared in `Suggests`).
#'
#' @export
plot_exclusive_junctions_event <- function(target_gene,
                                           GTF,
                                           eventdata,
                                           show_exclusive = TRUE,
                                           transcript     = NULL,
                                           curvature      = -0.2) {

  .ensure_ggplot("plot_exclusive_junctions_event")

  if (!data.table::is.data.table(eventdata))
    eventdata <- data.table::as.data.table(eventdata)

  if (!"gene_name" %in% names(eventdata)) {
    if (!is.character(GTF) || length(GTF) != 1L)
      stop("`eventdata` has no `gene_name` column; to run ",
           "make_eventdata_plus() `GTF` must be a single path.",
           call. = FALSE)
    message("eventdata has no `gene_name` column; running make_eventdata_plus().")
    eventdata <- make_eventdata_plus(eventdata, GTF_file_direction = GTF)
  }

  ed <- eventdata[gene_name == target_gene]
  if (!nrow(ed))
    stop(sprintf("No eventdata rows for gene '%s'.", target_gene),
         call. = FALSE)
  if (!all(c("i.start", "i.end") %in% names(ed)))
    stop("eventdata must contain `i.start` and `i.end` junction coordinates.",
         call. = FALSE)

  obs_junctions <- unique(ed[, .(j_start = i.start, j_end = i.end)])

  g        <- .load_gene_gtf(GTF, target_gene)
  exons    <- .gene_exons(g, target_gene)
  gene_id_val <- attr(exons, "gene_id")
  gtf_j    <- .tx_junctions(exons)

  junctions <- merge(gtf_j, obs_junctions, by = c("j_start", "j_end"))

  matched_obs     <- unique(junctions[, .(j_start, j_end)])
  novel_junctions <- data.table::fsetdiff(obs_junctions, matched_obs)
  novel_count     <- nrow(novel_junctions)
  if (novel_count > 0L)
    message(sprintf(
      "%d observed junction(s) did not match any annotated intron - dropped.",
      novel_count))

  if (!nrow(junctions))
    stop(sprintf("No observed junction matched any annotated transcript of '%s'.",
                 target_gene), call. = FALSE)

  all_tx   <- unique(junctions$transcript_name)
  tx_order <- .order_transcripts(junctions, all_tx,
                                 show_exclusive = show_exclusive,
                                 override_tx    = transcript)

  if (!nrow(tx_order))
    stop(sprintf("No transcripts to plot for '%s' with the given filters.",
                 target_gene), call. = FALSE)

  tx_id_map <- unique(exons[, .(transcript_name, transcript_id)])
  tx_map    <- merge(data.table::data.table(transcript_name = tx_order$transcript_name,
                                            y = seq_len(nrow(tx_order))),
                     tx_id_map, by = "transcript_name", sort = FALSE)
  data.table::setorder(tx_map, y)
  tx_map[, y_label := ifelse(transcript_name == transcript_id,
                             transcript_name,
                             paste0(transcript_name, "\n", transcript_id))]

  exons     <- merge(exons,     tx_map[, .(transcript_name, y)], by = "transcript_name")
  junctions <- merge(junctions, tx_map[, .(transcript_name, y)], by = "transcript_name")

  title_suffix <- if (!is.null(transcript))
                    sprintf(" (%s)", paste(transcript, collapse = ", "))
                  else if (isTRUE(show_exclusive))
                    " (exclusive transcripts only)"
                  else ""

  p <- .build_plot(exons, junctions, tx_map, target_gene,
                   curvature, title_suffix)
  # override caption for the event-based variant
  p <- p + ggplot2::labs(caption = paste0(
    "Arcs = junctions observed in eventdata. ",
    "Black: exclusive to one transcript; grey: shared."))

  info <- .build_info(junctions, exons, target_gene, gene_id_val, eventdata = ed)

  structure(
    list(plot = p, info = info,
         exons = exons, junctions = junctions, tx_order = tx_map,
         novel_junction_count = novel_count),
    class = c("splikit_junction_plot", "list")
  )
}


# --- S3 print method ---------------------------------------------------------

#' Print method for splikit junction-plot results
#'
#' Renders the stored plot and then prints the `info` data.table. Called
#' automatically at the REPL when a function like
#' [plot_exclusive_junctions()] is invoked without assignment.
#'
#' @param x A `splikit_junction_plot` object.
#' @param n Rows of `info` to preview. `NULL` (default) shows all.
#' @param ... Unused.
#'
#' @method print splikit_junction_plot
#' @export
print.splikit_junction_plot <- function(x, n = NULL, ...) {
  print(x$plot)
  cat("\n--- info (", nrow(x$info), " rows) ---\n", sep = "")
  if (is.null(n)) print(x$info) else print(utils::head(x$info, n))
  if (!is.null(x$novel_junction_count) && x$novel_junction_count > 0)
    cat("novel_junction_count:", x$novel_junction_count,
        "(observed junctions not matching any annotated intron)\n")
  invisible(x)
}
