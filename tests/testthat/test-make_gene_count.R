make_gene_count_fixture <- function(include_em = TRUE) {
  root <- tempfile("splikit-gene-")
  raw_dir <- file.path(root, "raw")
  filtered_dir <- file.path(root, "filtered")
  dir.create(raw_dir, recursive = TRUE)
  dir.create(filtered_dir, recursive = TRUE)

  raw_counts <- Matrix::Matrix(
    matrix(c(1, 2, 3, 4, 5, 6), nrow = 2), sparse = TRUE
  )
  em_counts <- Matrix::Matrix(
    matrix(c(11, 12, 13, 14, 15, 16), nrow = 2), sparse = TRUE
  )
  filtered_counts <- raw_counts[, c(1, 3), drop = FALSE]

  Matrix::writeMM(raw_counts, file.path(raw_dir, "matrix.mtx"))
  if (include_em) {
    Matrix::writeMM(em_counts, file.path(raw_dir, "UniqueAndMult-EM.mtx"))
  }
  Matrix::writeMM(filtered_counts, file.path(filtered_dir, "matrix.mtx"))

  data.table::fwrite(
    data.table::data.table(barcode = c("AA", "BB", "CC")),
    file.path(raw_dir, "barcodes.tsv"), col.names = FALSE
  )
  data.table::fwrite(
    data.table::data.table(barcode = c("AA", "CC")),
    file.path(filtered_dir, "barcodes.tsv"), col.names = FALSE
  )

  features <- data.table::data.table(
    id = c("gene1", "gene2"),
    name = c("Gene 1", "Gene 2"),
    type = "Gene Expression"
  )
  data.table::fwrite(features, file.path(raw_dir, "features.tsv"),
                     sep = "\t", col.names = FALSE)
  data.table::fwrite(features, file.path(filtered_dir, "features.tsv"),
                     sep = "\t", col.names = FALSE)

  list(
    root = root,
    raw_counts = raw_counts,
    em_counts = em_counts,
    filtered_counts = filtered_counts
  )
}

test_that("make_gene_count selects raw, filtered, and custom count matrices", {
  fixture <- make_gene_count_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  raw_result <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = FALSE,
    matrix_source = "raw",
    matrix_file = "matrix.mtx"
  )
  expect_equal(unname(as.matrix(raw_result)),
               unname(as.matrix(fixture$raw_counts)))
  expect_identical(colnames(raw_result), c("AA-sample1", "BB-sample1", "CC-sample1"))

  filtered_result <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = FALSE,
    matrix_source = "filtered"
  )
  expect_equal(unname(as.matrix(filtered_result)),
               unname(as.matrix(fixture$filtered_counts)))
  expect_identical(colnames(filtered_result), c("AA-sample1", "CC-sample1"))

  em_result <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = FALSE,
    matrix_source = "raw",
    matrix_file = "UniqueAndMult-EM.mtx"
  )
  expect_equal(unname(as.matrix(em_result)),
               unname(as.matrix(fixture$em_counts)))
})

test_that("matrix source is independent from barcode filtering", {
  fixture <- make_gene_count_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  result <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = TRUE,
    matrix_source = "raw",
    matrix_file = "matrix.mtx"
  )

  expect_equal(unname(as.matrix(result)),
               unname(as.matrix(fixture$raw_counts[, c(1, 3)])))
  expect_identical(colnames(result), c("AA-sample1", "CC-sample1"))

  em_result <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = TRUE,
    matrix_source = "raw",
    matrix_file = "UniqueAndMult-EM.mtx"
  )
  expect_equal(unname(as.matrix(em_result)),
               unname(as.matrix(fixture$em_counts[, c(1, 3)])))
  expect_identical(colnames(em_result), c("AA-sample1", "CC-sample1"))
})

test_that("defaults prefer raw EM counts and fall back to raw matrix.mtx", {
  fixture <- make_gene_count_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  default_result <- make_gene_count(fixture$root, "sample1")
  expect_equal(unname(as.matrix(default_result)),
               unname(as.matrix(fixture$em_counts[, c(1, 3)])))
  expect_identical(colnames(default_result), c("AA-sample1", "CC-sample1"))

  fallback <- make_gene_count_fixture(include_em = FALSE)
  on.exit(unlink(fallback$root, recursive = TRUE), add = TRUE)
  fallback_result <- make_gene_count(fallback$root, "sample1")
  expect_equal(unname(as.matrix(fallback_result)),
               unname(as.matrix(fallback$raw_counts[, c(1, 3)])))
})

test_that("auto matrix source preserves historical directory selection", {
  fixture <- make_gene_count_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  auto_filtered <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = TRUE,
    matrix_source = "auto",
    matrix_file = "matrix.mtx"
  )
  auto_raw <- make_gene_count(
    fixture$root, "sample1",
    use_internal_whitelist = FALSE,
    matrix_source = "auto",
    matrix_file = "matrix.mtx"
  )

  expect_equal(unname(as.matrix(auto_filtered)),
               unname(as.matrix(fixture$filtered_counts)))
  expect_equal(unname(as.matrix(auto_raw)),
               unname(as.matrix(fixture$raw_counts)))
})

test_that("make_gene_count validates matrix selection", {
  fixture <- make_gene_count_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)

  expect_error(
    make_gene_count(fixture$root, "sample1", matrix_source = "unknown"),
    "arg.*should be one of"
  )
  expect_error(
    make_gene_count(fixture$root, "sample1", matrix_file = ""),
    "single non-empty filename"
  )
  expect_error(
    make_gene_count(
      fixture$root, "sample1",
      matrix_source = "raw", matrix_file = "missing.mtx"
    ),
    "missing\\.mtx"
  )
})
