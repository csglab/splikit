.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    "=======================================================================================================\n",
    "  Welcome to Splikit — Developed by the Computational and Statistical Genomics Laboratory\n",
    "  at McGill University, Montreal. Distributed under the MIT License.\n",
    "  For documentation and examples, see the package vignette.\n",
    "-------------------------------------------------------------------------------------------------------\n",
    "  Bienvenue dans Splikit — Developpe par le Laboratoire de genomique computationnelle et statistique\n",
    "  de l’Universite McGill, a Montreal. Distribue sous licence MIT.\n",
    "  Pour la documentation et des exemples, consultez la vignette du package.\n",
    "=======================================================================================================\n"
  )
}


.onLoad <- function(libname, pkgname) {
  required_pkgs <- c("Rcpp", "RcppEigen", "RcppArmadillo", "Matrix", "data.table")
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs)) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }
}
