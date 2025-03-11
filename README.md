# Splicing Pipeline
[![Documentation](https://img.shields.io/badge/Docs-Learn%20More-blue.svg)](./docs/README.md)


This pipeline provides a streamlined approach for processing and analyzing splicing data from single-cell RNA sequencing. It begins with mapping raw reads using **STARsolo**, capturing splice junctions in the `SJ.out.tab` file. Next, junctions sharing intron coordinates are grouped into Local Junction Variants (LJVs), and counts are split into **M1** (junction-specific reads) and **M2** (reads supporting alternative junctions). A Poisson-based thinning method helps split these counts into independent train/test datasets, ensuring unbiased downstream analyses. Additional tools within the pipeline facilitate generating gene expression matrices, combining splicing and gene-level measurements, and preparing Velocyto (spliced/unspliced) data, offering a comprehensive workflow for single-cell splicing quantification.

### **Requirments**
---
- Unix-compatible OS.  
- [R version 3.0.0](http://www.r-project.org/) or later. 
- R libraries: [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [countsplit](https://cran.r-project.org/web/packages/countsplit/index.html), [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

### **Key Functions**  
---
- `multigedi_make_junction_ab`: Generate per-sample junction abundance.  
- `multigedi_make_m1`: Combine M1 matrices with event data.  
- `multigedi_countsplit`: Split M1 into train/test sets.  
- `multigedi_make_m2`: Create the M2 (exclusion) matrix.  
- `multigedi_make_gene` & `multigedi_make_velo`: Generate gene expression and Velocyto (spliced/unspliced) matrices.
