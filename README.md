# Splikit

*A toolkit for analyzing RNA splicing in single-cell RNA-seq (scRNA-seq) data.*

<p align="right">

<img src="./docs/logo.png" alt="Splikit logo" width="250"/>

</p>

[![Documentation](https://img.shields.io/badge/Docs-Learn%20More-blue.svg)](./docs/README.md)

### **Requirments**

---
- Unix-compatible OS.  
- [R version 3.0.0](http://www.r-project.org/) or later. 
- R libraries: [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html), [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html)

### **Key Functions**  
---

-   `multigedi_make_junction_ab`: Generate per-sample junction abundance.\
-   `multigedi_make_m1`: Combine M1 matrices with event data.\
-   `multigedi_countsplit`: Split M1 into train/test sets.\
-   `multigedi_make_m2`: Create the M2 (exclusion) matrix.\
-   `multigedi_make_gene` & `multigedi_make_velo`: Generate gene expression and Velocyto (spliced/unspliced) matrices.
