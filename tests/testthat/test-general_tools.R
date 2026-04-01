test_that("get_pseudo_correlation works on valid data", {
  n_events <- 5
  n_cells <- 100
  set.seed(42)
  
  m1 <- matrix(rbinom(n_events * n_cells, size = 10, prob = 0.5), nrow = n_events)
  m2 <- matrix(rbinom(n_events * n_cells, size = 10, prob = 0.5), nrow = n_events)
  ZDB <- matrix(rnorm(n_events * n_cells), nrow = n_events, ncol = n_cells)
  
  rownames(m1) <- rownames(m2) <- rownames(ZDB) <- paste0("event_", 1:n_events)
  colnames(m1) <- colnames(m2) <- colnames(ZDB) <- paste0("cell_", 1:n_cells)
  
  # Test with default arguments
  res <- get_pseudo_correlation(ZDB, m1, m2)
  expect_true(data.table::is.data.table(res))
  expect_equal(nrow(res), n_events)
  expect_true(all(c("event", "pseudo_correlation", "null_distribution") %in% names(res)))
  expect_type(res$pseudo_correlation, "double")
  
  # Test with Nagelkerke metric
  res_nagel <- get_pseudo_correlation(ZDB, m1, m2, metric = "Nagelkerke")
  expect_true(data.table::is.data.table(res_nagel))
  expect_equal(nrow(res_nagel), n_events)
  
  # Test with sparse matrices
  m1_sp <- Matrix::Matrix(m1, sparse = TRUE)
  m2_sp <- Matrix::Matrix(m2, sparse = TRUE)
  res_sp <- get_pseudo_correlation(ZDB, m1_sp, m2_sp)
  expect_equal(res_sp$pseudo_correlation, res$pseudo_correlation)
  
  # Test with SplikitObject wrapper
  m1_obj <- list(m1 = m1_sp, m2 = m2_sp)
  class(m1_obj) <- "SplikitObject"
  res_obj <- get_pseudo_correlation(m1_obj, m1_inclusion = ZDB)
  expect_equal(res_obj$pseudo_correlation, res$pseudo_correlation)
})

test_that("get_pseudo_correlation error handling works", {
  m1 <- matrix(1, nrow = 2, ncol = 2)
  m2 <- matrix(1, nrow = 2, ncol = 2)
  ZDB <- matrix(1, nrow = 3, ncol = 2) # wrong rows
  
  expect_error(get_pseudo_correlation(ZDB, m1, m2), "must have the same number of rows")
  
  ZDB_bad <- data.frame(ZDB)
  expect_error(get_pseudo_correlation(ZDB_bad, m1, m2), "ZDB_matrix must be a dense matrix")
})

test_that("get_rowVar works correctly", {
  set.seed(42)
  dm <- matrix(rnorm(100), nrow = 10)
  
  # Calculate row variance using apply and var (sample variance)
  expected_sample_vars <- apply(dm, 1, var)
  # get_rowVar computes population variance
  expected_vars <- expected_sample_vars * (ncol(dm) - 1) / ncol(dm)
  
  # Using get_rowVar
  actual_vars <- get_rowVar(dm)
  
  expect_equal(actual_vars, expected_vars)
  
  # Sparse test
  sm <- Matrix::Matrix(dm, sparse = TRUE)
  actual_vars_sp <- get_rowVar(sm)
  expect_equal(actual_vars_sp, expected_vars)
  
  expect_error(get_rowVar(c(1, 2, 3)), "must be either a dense numeric matrix or a dgCMatrix")
})

test_that("get_silhouette_mean works correctly", {
  set.seed(42)
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  clusters <- as.integer(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2))
  
  # Test execution without error
  score <- get_silhouette_mean(X, clusters)
  expect_true(is.numeric(score))
  expect_length(score, 1)
  
  # Test parallel
  score_p <- get_silhouette_mean(X, clusters, n_threads = 2)
  expect_equal(score_p, score)
  
  # Dimension failures
  expect_error(get_silhouette_mean(X, clusters[1:5]))
})
