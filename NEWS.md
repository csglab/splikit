# splikit 2.3.1

## CRAN Compliance
* Fixed `no visible binding for '..keep'` NOTE in `.gene_exons()`: replaced
  `exons[, ..keep]` with `exons[, keep, with = FALSE]`.
* Fixed `no visible binding for 'x'` NOTE in `.build_plot()`: added `x` to
  `globalVariables()` in `globals.R`.
* Wrapped `get_pseudo_correlation()` example in `\donttest{}` to avoid
  exceeding the CRAN 5-second example time limit.

# splikit 2.3.0

## New Features
* **RefSeq GTF support in plotting functions.**
  `plot_exclusive_junctions()` and `plot_exclusive_junctions_event()` now
  accept RefSeq-style GTFs that lack `gene_name` / `transcript_name`
  attributes. The internal attribute extractor tries `gene_name` ->
  `gene_id` -> `gene`, and `transcript_name` -> `transcript_id`, so a
  RefSeq annotation like NCBI's GCF_*_GRCh38 file works with no manual
  preprocessing.
* **Structured return value.**
  Both plotting functions now return an S3 object of class
  `splikit_junction_plot` with components `$plot` (ggplot),
  `$info` (tidy per-junction data.table), plus the underlying
  `exons` / `junctions` / `tx_order` tables. The `info` data.table
  includes `gene_name`, `gene_id`, `transcript_name`, `transcript_id`,
  `chr`, `strand`, `j_start`, `j_end`, `j_width`, `exclusive`,
  `n_tx_with_junction`, `observed_in_eventdata`, `row_names_mtx`, and
  `is_annot` for straightforward downstream analysis.
* **S3 `print()` method** for the new class: rendering the plot and
  printing `info` when the object is auto-printed at the REPL.

## Internals
* Extracted shared helpers (`.gtf_attr`, `.load_gene_gtf`, `.gene_exons`,
  `.tx_junctions`, `.order_transcripts`, `.build_plot`, `.build_info`)
  so the GTF-only and eventdata-based plots now share a single
  implementation path.

# splikit 2.2.1

## CRAN Compliance
* Addressed CRAN reviewer feedback on the 2.0.0 submission: replaced the
  un-suppressable `print(summary_table)` calls in `make_junction_ab()`,
  `make_gene_count()`, and `make_velo_count()` with
  `message(paste(capture.output(...), collapse = "\n"))`. The summary now
  respects `suppressMessages()` and stays behind the `verbose = TRUE`
  guard. A stray `if (verbose) cat(...)` in `make_m1()` was likewise
  converted to `message()` for consistency.
* Title and Description already had the redundant "A toolkit designed
  for" phrasing removed in 2.0.1; no further action needed there.

# splikit 2.2.0

## New Features
* **`plot_exclusive_junctions_event()`** - sibling of
  `plot_exclusive_junctions()` that sources the drawn arcs from a
  `splikit` eventdata table. Exon structure still comes from the GTF and
  exclusivity is computed gene-wide from the annotation, but only junctions
  whose `(i.start, i.end)` coordinates are observed in `eventdata` are
  drawn. If `eventdata` lacks a `gene_name` column the function runs
  `make_eventdata_plus()` internally using the supplied GTF path. The plot
  is rendered directly to the current graphics device (no `out_file`).

# splikit 2.1.0

## New Features
* **Transcript-exclusive junction plotting**
  - `plot_exclusive_junctions()` renders a gene model with one row per
    transcript, exon rectangles, intron/junction arcs, and Ensembl ids
    shown beneath transcript names on the y-axis. Junctions used by a
    single transcript of the gene are highlighted as solid black arcs;
    shared junctions render as thin grey arcs.
  - Flags: `show_exclusive` (restrict to exclusive-owning transcripts),
    `transcript` (pin the plot to one or more transcript names),
    `curvature` (arc-height knob).
  - `plot_exclusive_junctions_pdf()` writes a multi-page PDF: page 1 is
    the full gene view, subsequent pages show each exclusive transcript
    with its exclusive junction in black.
  - Accepts a GTF path or an already-loaded `data.table`. `ggplot2` is
    declared under `Suggests`.

## Bug Fixes
* **`make_eventdata_plus()` gene_name extraction**
  - The previous regex used `sub('.*gene_name "([^"]+)".*', ...)`, which
    returns the input unchanged when no match is found. Gene records
    without a `gene_name` attribute (e.g. ~310 rows in Ensembl mouse
    GRCm39.110) therefore ended up with a `gene_name` equal to the full
    attribute blob.
  - Fixed with a `grepl()`-guarded extractor that returns `NA` on
    no-match, with a fallback to `gene_id` so downstream matching by
    name keeps working. The function now emits a `message()` reporting
    how many records required the fallback.

# splikit 1.0.5

## 🚀 Major Enhancements

### Performance Improvements
* **Multi-threading Support**: Added parallel processing capabilities for gene expression analysis
  - `find_variable_genes()` now supports `n_threads` parameter for faster computation with `method = "sum_deviance"`
  - Consistent multi-threading implementation across all feature selection functions
  - Significant speedup for large-scale datasets

* **Intelligent Memory Management**: Enhanced memory efficiency for large datasets
  - `make_m2()` now includes automatic detection of memory limits
  - Smart switching between fast and batched processing modes
  - New parameters: `batch_size`, `memory_threshold`, `force_fast` for fine-tuned control
  - Prevents 32-bit integer overflow issues with large matrices

### New Features
* **Enhanced Pseudo-correlation Analysis**
  - `get_pseudo_correlation()` now supports both Cox-Snell and Nagelkerke R² metrics
  - Added `metric` parameter for selecting correlation type (default: "CoxSnell")
  - Full support for both sparse and dense matrix inputs for M1/M2 matrices
  - Improved numerical stability with singular matrix handling

* **Improved Data Processing**
  - `make_junction_ab()` enhanced with:
    - `verbose` parameter for detailed progress tracking
    - `keep_multi_mapped_junctions` option for including multi-mapped reads
    - Better memory management for large junction datasets
  - `make_m1()` now includes:
    - `min_counts` parameter for filtering low-count events
    - `verbose` option for progress monitoring
    - Improved handling of sparse events

### Developer Experience
* **Enhanced Logging**: All major functions now support verbose output
  - Detailed progress tracking for long-running operations
  - Memory usage reporting in `make_m2()`
  - Informative messages even with `verbose = FALSE` for critical steps

* **Documentation Website**: Added pkgdown support
  - Automatic documentation generation via GitHub Actions
  - Bootstrap 5 modern interface
  - Complete function reference with examples
  - Available at: https://arshammik.github.io/splikit/

## 🐛 Bug Fixes
* Fixed non-ASCII characters in verbose output messages (replaced with ASCII alternatives)
* Corrected parameter naming inconsistency (`multithread` → `multi_thread`)
* Removed unused Nagelkerke R² calculation in pseudo-correlation function
* Fixed gfortran dependency issues in Makevars configuration

## 📚 Documentation Updates
* Comprehensive updates to all function documentation with new parameters
* Added detailed vignettes:
  - Complete splikit manual for single-cell splicing analysis
  - STARsolo processing guide with step-by-step instructions
* Updated README with clearer examples and installation instructions
* Consistent Roxygen2 documentation across all functions

## 🔧 Technical Improvements
* Updated to version 1.0.5
* Improved C++ implementations for better performance
* Better error handling and input validation
* Memory-efficient sparse matrix operations
* Cleaned up package structure for CRAN submission readiness

## Breaking Changes
None - all changes are backward compatible. Existing code will continue to work with default parameters.

## Acknowledgments
This release includes contributions addressing issue #16, focusing on performance enhancements and multi-threading support. Special thanks to all users who provided feedback and testing.

## Installation
```r
# Install from GitHub
devtools::install_github("Arshammik/splikit")

# Load the package
library(splikit)
```

## What's Next
* S4 class implementation for object-oriented interface (planned for 2.0.0)
* Additional statistical methods for splicing analysis
* Enhanced visualization capabilities
* Further performance optimizations

---
*For questions or issues, please visit: https://github.com/Arshammik/splikit/issues*