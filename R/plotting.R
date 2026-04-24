#' Plot transcript-exclusive splice junctions
#'
#' @description
#' Draws a gene-model view of a gene's transcripts (one row per transcript)
#' with splice junctions rendered as arcs above each row. Junctions used by a
#' single transcript of the gene ("transcript-exclusive") are highlighted in
#' solid black; junctions shared with other transcripts are drawn as thin grey
#' arcs. The y-axis labels show `transcript_name` over `transcript_id`
#' (two lines), so the Ensembl id is visible next to the usual transcript name.
#'
#' @details
#' Exclusivity is computed against every transcript of the gene in the GTF,
#' regardless of which transcripts end up being drawn. That means the black
#' arcs always reflect gene-wide exclusivity, even when the view is restricted
#' via `show_exclusive = TRUE` or `transcript = "..."`.
#'
#' Plotting behaviour by flag:
#' * `show_exclusive = TRUE` (default): restrict the plot to transcripts that
#'   own at least one exclusive junction.
#' * `show_exclusive = FALSE`: draw every transcript of the gene.
#' * `transcript = "Tpm1-205"` (or a vector of names): pin the plot to those
#'   transcripts; overrides `show_exclusive`.
#'
#' @param gtf Either a path (character) to an Ensembl-style GTF file, or a
#'   `data.table` already loaded from one. When a path is given, a fast
#'   `grep`-based prefilter is tried first (Linux/macOS); it falls back to a
#'   full `fread()` if the prefilter is unavailable or returns nothing. When a
#'   `data.table` is passed it must contain the nine standard GTF columns
#'   (either `seqname`/`attr` or `seqid`/`attribute` naming is accepted).
#' @param target_gene Character. Gene name (matches the `gene_name` attribute
#'   in the GTF), e.g. `"Tpm1"`.
#' @param show_exclusive Logical (default `TRUE`). If `TRUE`, only transcripts
#'   that own at least one exclusive junction are plotted; if `FALSE`, every
#'   transcript of the gene is drawn.
#' @param transcript Character vector of transcript names to pin the plot to
#'   (e.g. `"Tpm1-205"` or `c("Tpm1-220", "Tpm1-221")`). When supplied this
#'   overrides `show_exclusive`. Exclusivity is still computed against the
#'   whole gene.
#' @param curvature Numeric. Arc-height knob passed to
#'   `ggplot2::geom_curve()`. Values near `0` yield flat arcs; more negative
#'   values yield taller arcs. Default `-0.2`.
#' @param out_file Optional character. If given, the plot is written to this
#'   file via `ggplot2::ggsave()` (format inferred from the extension).
#'
#' @return A list with elements `plot` (the `ggplot` object), `junctions`
#'   (per-transcript junction table with an `exclusive` column), `exons` (the
#'   exon table used for drawing), and `tx_order` (transcript-to-y mapping
#'   with Ensembl ids).
#'
#' @section Required packages:
#' Requires the \pkg{ggplot2} package, listed in `Suggests`. Install it before
#' calling this function.
#'
#' @examples
#' \dontrun{
#'   res <- plot_exclusive_junctions(
#'     gtf         = "Mus_musculus.GRCm39.110.gtf",
#'     target_gene = "Tpm1",
#'     out_file    = "Tpm1_exclusive_junctions.png"
#'   )
#'   # pin to a single transcript
#'   plot_exclusive_junctions(gtf = res_gtf_dt, target_gene = "Tpm1",
#'                            transcript = "Tpm1-220")
#' }
#'
#' @export
plot_exclusive_junctions <- function(gtf,
                                     target_gene,
                                     show_exclusive = TRUE,
                                     transcript     = NULL,
                                     curvature      = -0.2,
                                     out_file       = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_exclusive_junctions(). ",
         "Install it with install.packages('ggplot2').", call. = FALSE)
  }

  # --- 1. Obtain the gene's GTF rows ---
  g <- .load_gene_gtf(gtf, target_gene)

  if (!nrow(g))
    stop(sprintf("No GTF rows matched gene_name = '%s'.", target_gene), call. = FALSE)

  get_attr <- function(x, key) {
    re  <- sprintf('.*%s "([^"]+)".*', key)
    hit <- grepl(re, x)
    out <- rep(NA_character_, length(x))
    out[hit] <- sub(re, "\\1", x[hit])
    out
  }

  g[, gene_name := get_attr(attr, "gene_name")]
  g <- g[gene_name == target_gene]

  exons <- g[type == "exon"]
  exons[, transcript_id   := get_attr(attr, "transcript_id")]
  exons[, transcript_name := get_attr(attr, "transcript_name")]
  exons[, exon_number     := as.integer(get_attr(attr, "exon_number"))]
  exons <- exons[, .(seqname, strand, transcript_id, transcript_name,
                     exon_number, start, end)]

  if (!nrow(exons))
    stop(sprintf("No exon rows found for gene '%s'.", target_gene), call. = FALSE)

  # --- 2. Per-transcript junctions (sort by genomic position, strand-agnostic) ---
  data.table::setorder(exons, transcript_id, start)
  junctions <- exons[, {
    if (.N >= 2) .(j_start = end[-.N] + 1L, j_end = start[-1] - 1L)
    else         .(j_start = integer(0),    j_end = integer(0))
  }, by = .(transcript_id, transcript_name, strand, seqname)]

  # --- 3. Exclusive = used by exactly one transcript (gene-wide) ---
  j_counts  <- junctions[, .(n_tx = data.table::uniqueN(transcript_id)),
                         by = .(j_start, j_end)]
  junctions <- merge(junctions, j_counts, by = c("j_start", "j_end"))
  junctions[, exclusive := n_tx == 1L]
  junctions[, targeted := exclusive]

  # --- 4. Pick / order transcripts ---
  all_tx   <- unique(exons$transcript_name)
  tx_order <- junctions[, .(has_excl = any(exclusive)), by = transcript_name
                        ][order(-has_excl, transcript_name)]

  if (!is.null(transcript)) {
    missing <- setdiff(transcript, all_tx)
    if (length(missing))
      stop(sprintf("Transcript(s) not in %s: %s",
                   target_gene, paste(missing, collapse = ", ")), call. = FALSE)
    extra <- setdiff(transcript, tx_order$transcript_name)
    if (length(extra))
      tx_order <- rbind(tx_order,
                        data.table::data.table(transcript_name = extra,
                                               has_excl = FALSE))
    tx_order <- tx_order[transcript_name %in% transcript]
    tx_order <- tx_order[match(transcript, transcript_name)]
  } else if (isTRUE(show_exclusive)) {
    tx_order <- tx_order[has_excl == TRUE]
  } else {
    tx_order <- rbind(tx_order,
                      data.table::data.table(
                        transcript_name = setdiff(all_tx, tx_order$transcript_name),
                        has_excl = FALSE))
  }

  tx_id_map <- unique(exons[, .(transcript_name, transcript_id)])
  tx_map    <- merge(data.table::data.table(transcript_name = tx_order$transcript_name,
                                            y = seq_len(nrow(tx_order))),
                     tx_id_map, by = "transcript_name", sort = FALSE)
  data.table::setorder(tx_map, y)
  tx_map[, y_label := paste0(transcript_name, "\n", transcript_id)]

  exons     <- merge(exons,     tx_map[, .(transcript_name, y)], by = "transcript_name")
  junctions <- merge(junctions, tx_map[, .(transcript_name, y)], by = "transcript_name")

  if (!nrow(tx_map))
    stop(sprintf("No transcripts to plot for '%s' with the given filters.",
                 target_gene), call. = FALSE)

  intron_spans <- exons[, .(x = min(start), xend = max(end)),
                        by = .(transcript_name, y)]
  chrom <- unique(exons$seqname)[1]

  # --- 5. Build plot ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = intron_spans,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          color = "grey55", linewidth = 0.35) +
    ggplot2::geom_rect(data = exons,
                       ggplot2::aes(xmin = start, xmax = end,
                                    ymin = y - 0.28, ymax = y + 0.28),
                       fill = "grey85", color = "black", linewidth = 0.25)

  if (nrow(junctions[targeted == FALSE]) > 0)
    p <- p + ggplot2::geom_curve(
      data = junctions[targeted == FALSE],
      ggplot2::aes(x = j_start, xend = j_end,
                   y = y + 0.28, yend = y + 0.28),
      curvature = curvature, linewidth = 0.25, color = "grey70")

  if (nrow(junctions[targeted == TRUE]) > 0)
    p <- p + ggplot2::geom_curve(
      data = junctions[targeted == TRUE],
      ggplot2::aes(x = j_start, xend = j_end,
                   y = y + 0.28, yend = y + 0.28),
      curvature = curvature, linewidth = 0.7, color = "black")

  subtitle_bit <- if (!is.null(transcript))
                    sprintf(" (%s)", paste(transcript, collapse = ", "))
                  else if (isTRUE(show_exclusive))
                    " (exclusive transcripts only)"
                  else ""

  p <- p +
    ggplot2::scale_y_continuous(breaks = tx_map$y, labels = tx_map$y_label,
                                expand = ggplot2::expansion(add = 0.9)) +
    ggplot2::labs(
      x = sprintf("%s position (bp)", chrom), y = NULL,
      title = sprintf("%s - transcript-exclusive junctions%s",
                      target_gene, subtitle_bit),
      caption = "Black arcs: junctions used by a single transcript. Grey arcs: shared junctions."
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 8, lineheight = 0.9),
      plot.caption       = ggplot2::element_text(hjust = 0, size = 8,
                                                 color = "grey30")
    )

  if (!is.null(out_file))
    ggplot2::ggsave(out_file, p,
                    width  = 10,
                    height = max(3, 0.45 * nrow(tx_map) + 2),
                    dpi    = 150)

  invisible(list(plot = p, junctions = junctions, exons = exons, tx_order = tx_map))
}

#' Multi-page PDF of transcript-exclusive junctions
#'
#' @description
#' Writes a multi-page PDF summarising a gene's splicing variation. Page one
#' shows the full gene (every transcript) with exclusive junctions highlighted
#' in black. Each subsequent page shows one transcript that owns an exclusive
#' junction, with that junction rendered in black and its other junctions in
#' grey.
#'
#' @param gtf Same as in [plot_exclusive_junctions()] - a path or a pre-loaded
#'   GTF `data.table`.
#' @param target_gene Gene name, e.g. `"Tpm1"`.
#' @param out_pdf Character. Path to the PDF file to write.
#' @param curvature Numeric. Arc-height knob. Default `-0.2`.
#' @param width Numeric. PDF page width in inches. Default `10`.
#' @param height Numeric. PDF page height in inches. If `NULL` (default), it
#'   is sized automatically to fit the full-gene page (all transcripts).
#'
#' @return Invisibly returns a list with `exclusive_transcripts` (the names of
#'   the exclusive-owning transcripts, in the order they were drawn) and
#'   `n_pages` (number of pages written).
#'
#' @section Required packages:
#' Requires the \pkg{ggplot2} package, listed in `Suggests`.
#'
#' @examples
#' \dontrun{
#'   plot_exclusive_junctions_pdf(
#'     gtf         = "Mus_musculus.GRCm39.110.gtf",
#'     target_gene = "Tpm1",
#'     out_pdf     = "Tpm1_variations.pdf"
#'   )
#' }
#'
#' @export
plot_exclusive_junctions_pdf <- function(gtf, target_gene, out_pdf,
                                         curvature = -0.2,
                                         width = 10, height = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_exclusive_junctions_pdf(). ",
         "Install it with install.packages('ggplot2').", call. = FALSE)
  }

  full <- plot_exclusive_junctions(gtf, target_gene,
                                   show_exclusive = FALSE,
                                   curvature      = curvature)
  excl_tx <- unique(full$junctions[exclusive == TRUE, transcript_name])

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


# --- internal helper ---------------------------------------------------------

# Load only the GTF rows relevant to `target_gene`. Accepts a path or an
# already-loaded data.table; normalises column names to (seqname, ..., attr).
.load_gene_gtf <- function(gtf, target_gene) {
  gtf_cols <- c("seqname", "source", "type", "start", "end",
                "score", "strand", "frame", "attr")

  if (data.table::is.data.table(gtf)) {
    g <- data.table::copy(gtf)
  } else {
    # Try fast grep prefilter (Unix-like). Fall back to a full fread.
    g <- tryCatch(
      data.table::fread(
        cmd = sprintf("grep -F 'gene_name \"%s\"' %s",
                      target_gene, shQuote(path.expand(gtf))),
        sep = "\t", header = FALSE, col.names = gtf_cols, quote = "",
        showProgress = FALSE
      ),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    if (is.null(g) || !nrow(g)) {
      g <- data.table::fread(
        path.expand(gtf),
        sep = "\t", header = FALSE, col.names = gtf_cols, quote = "",
        skip = "#", showProgress = FALSE
      )
    }
  }

  # normalise alternate column names used elsewhere in splikit
  if (!"seqname" %in% names(g) && "seqid"     %in% names(g))
    data.table::setnames(g, "seqid", "seqname")
  if (!"attr"    %in% names(g) && "attribute" %in% names(g))
    data.table::setnames(g, "attribute", "attr")

  g
}
