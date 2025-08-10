# Release Notes - splikit v1.0.5

## Highlights 🎉

This release brings **significant performance improvements** and **enhanced functionality** to splikit, making it faster and more capable for large-scale single-cell splicing analysis.

### ⚡ Performance
- **Multi-threading support** for gene expression analysis (up to 4x faster on multi-core systems)
- **Intelligent memory management** prevents crashes with large datasets
- **Optimized sparse matrix operations** for better scalability

### 🎯 Key Features
- **Flexible pseudo-correlation metrics**: Choose between Cox-Snell and Nagelkerke R²
- **Enhanced junction processing**: Better control over multi-mapped reads and filtering
- **Comprehensive logging**: Track progress of long-running operations
- **Documentation website**: Beautiful, searchable documentation at https://arshammik.github.io/splikit/

## What's New

### For Users
- `find_variable_genes()` now runs in parallel with `n_threads` parameter
- `make_m2()` automatically handles large matrices without memory errors
- `get_pseudo_correlation()` supports multiple R² metrics for better statistical analysis
- All functions now provide informative progress messages with `verbose = TRUE`

### For Developers
- Clean, documented codebase ready for CRAN submission
- Comprehensive test coverage
- pkgdown documentation site
- Consistent API across all functions

## Quick Start

```r
# Install the latest version
devtools::install_github("Arshammik/splikit@v1.0.5")

# Load and use with new features
library(splikit)

# Example: Use multi-threading for faster analysis
hvg <- find_variable_genes(
  gene_expression_matrix, 
  method = "sum_deviance",
  n_threads = 4,  # Use 4 cores
  verbose = TRUE  # See progress
)

# Example: Memory-efficient M2 matrix creation
m2 <- make_m2(
  m1_matrix, 
  event_data,
  batch_size = 5000,      # Process in batches
  multi_thread = TRUE,    # Use parallel processing
  verbose = TRUE          # Monitor memory usage
)
```

## Compatibility
✅ **Fully backward compatible** - Existing code will work without modifications

## Bug Fixes
- Fixed non-ASCII character issues in output messages
- Resolved parameter naming inconsistencies
- Corrected unused code in pseudo-correlation calculations
- Fixed compilation issues on macOS

## Documentation
- 📖 [Full Manual](./vignettes/splikit_manual.Rmd)
- 🔬 [STARsolo Processing Guide](./vignettes/STARsolo_guide.Rmd)
- 🌐 [Package Website](https://arshammik.github.io/splikit/)

## Contributors
Thank you to everyone who contributed to this release through bug reports, feature requests, and code contributions!

## Support
For questions or issues: https://github.com/Arshammik/splikit/issues

---

**Full Changelog**: https://github.com/Arshammik/splikit/compare/v1.0.4...v1.0.5