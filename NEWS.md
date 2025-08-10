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