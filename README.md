# Splikit <img src="/docs/logo.png" align="right" height="220"/></a>

**Splikit** /ˈsplaɪ.kɪt/ is a toolkit designed for the analysis of high-dimensional single-cell
splicing data. It provides a framework to extract and work with ratio-based data structures
derived from single-cell RNA sequencing experiments. The package avoids the need for bulky
S4 objects by offering direct and efficient manipulation of matrices. Core functionalities
are implemented in C++ via Rcpp to ensure high performance and scalability on large datasets.

[![Documentation](https://img.shields.io/badge/Docs-Learn%20More-blue.svg)](./docs/README.md)

### **Requirments**

-   Unix-compatible OS.
-   [R version 3.0.0](http://www.r-project.org/) or later.
-   R libraries: [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html), [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
