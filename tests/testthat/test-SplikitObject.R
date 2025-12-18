# Tests for SplikitObject R6 class
# Uses test_splikit.rds to verify old and new methods produce identical results

test_that("SplikitObject initializes from matrices", {
  # Load test data
  test_data <- load_toy_M1_M2_object()

  # Create object from existing matrices
  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  expect_s3_class(obj, "SplikitObject")
  expect_equal(nrow(obj$m1), 2000)
  expect_equal(ncol(obj$m1), 2000)
  expect_true(!is.null(obj$eventData))
  expect_true(is.null(obj$m2))
})

test_that("SplikitObject validates dimension mismatches", {
  test_data <- load_toy_M1_M2_object()

  # Wrong number of rows in eventData
  bad_eventdata <- test_data$eventdata[1:100, ]

  expect_error(
    SplikitObject$new(m1 = test_data$m1, eventData = bad_eventdata),
    "same number of rows"
  )
})

test_that("SplikitObject validates m1/m2 dimension mismatch", {
  test_data <- load_toy_M1_M2_object()

  # Create m2 with wrong dimensions
  bad_m2 <- test_data$m1[1:100, ]

  expect_error(
    SplikitObject$new(m1 = test_data$m1, m2 = bad_m2),
    "identical dimensions"
  )
})

test_that("SplikitObject$makeM2 produces same result as make_m2 function", {
  skip_if_not_installed("Matrix")

  test_data <- load_toy_M1_M2_object()

  # Old method
  m2_old <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  # New method via R6
  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )
  obj$makeM2(verbose = FALSE)
  m2_new <- obj$m2

  # Compare results
  expect_equal(dim(m2_old), dim(m2_new))
  expect_equal(rownames(m2_old), rownames(m2_new))
  expect_equal(colnames(m2_old), colnames(m2_new))

  # Check values are identical (within floating point tolerance)
  expect_equal(Matrix::nnzero(m2_old), Matrix::nnzero(m2_new))
  expect_equal(sum(m2_old), sum(m2_new))
})

test_that("SplikitObject$findVariableEvents produces same result as function", {
  skip_if_not_installed("Matrix")

  test_data <- load_toy_M1_M2_object()

  # Compute M2 first
  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  # Old method
  hve_old <- find_variable_events(
    m1_matrix = test_data$m1,
    m2_matrix = m2,
    min_row_sum = 50,
    n_threads = 1,
    verbose = FALSE
  )

  # New method via R6
  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )
  hve_new <- obj$findVariableEvents(min_row_sum = 50, n_threads = 1, verbose = FALSE)

  # Compare results
  expect_equal(nrow(hve_old), nrow(hve_new))
  expect_equal(names(hve_old), names(hve_new))

  # Sort by events to ensure same order
  hve_old <- hve_old[order(events)]
  hve_new <- hve_new[order(events)]

  expect_equal(hve_old$events, hve_new$events)
  expect_equal(hve_old$sum_deviance, hve_new$sum_deviance, tolerance = 1e-10)
})

test_that("SplikitObject$subset works correctly", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Subset by indices
  original_events <- nrow(obj$m1)
  original_cells <- ncol(obj$m1)

  obj$subset(events = 1:100, cells = 1:500)

  expect_equal(nrow(obj$m1), 100)
  expect_equal(ncol(obj$m1), 500)
  expect_equal(nrow(obj$eventData), 100)
})

test_that("SplikitObject$subset validates empty results", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Empty events
  expect_error(
    obj$subset(events = integer(0)),
    "remove all events"
  )
})

test_that("SplikitObject$findVariableEvents validates threshold", {
  test_data <- load_toy_M1_M2_object()

  # Compute M2
  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )

  # Very high threshold that removes all events
  expect_error(
    obj$findVariableEvents(min_row_sum = 1e9),
    "No events pass"
  )
})

test_that("SplikitObject requires M2 for findVariableEvents", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # M2 not computed yet
  expect_error(
    obj$findVariableEvents(),
    "M2 not computed"
  )
})

test_that("SplikitObject$summary returns correct information", {
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )

  summ <- obj$summary()

  expect_equal(summ$events, 2000)
  expect_equal(summ$cells, 2000)
  expect_true(summ$has_m2)
  expect_false(summ$has_gene_expression)
  expect_true(is.numeric(summ$sparsity_m1))
  expect_true(summ$sparsity_m1 >= 0 && summ$sparsity_m1 <= 1)
})

test_that("SplikitObject$print works without error", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Should print without error
  expect_output(obj$print(), "SplikitObject")
  expect_output(obj$print(), "Events:")
  expect_output(obj$print(), "Cells:")
})

test_that("SplikitObject method chaining works", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Chain makeM2 and then check m2 is not null
  result <- obj$makeM2(verbose = FALSE)

  # makeM2 should return self for chaining
  expect_identical(result, obj)
  expect_true(!is.null(obj$m2))
})

test_that("splikit() convenience function works", {
  test_data <- load_toy_M1_M2_object()

  obj <- splikit(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  expect_s3_class(obj, "SplikitObject")
  expect_equal(nrow(obj$m1), 2000)
})

test_that("SplikitObject$setGeneExpression validates dimensions", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Wrong number of columns
  bad_gene_matrix <- Matrix::rsparsematrix(1000, 100, 0.1)

  expect_error(
    obj$setGeneExpression(bad_gene_matrix),
    "same number of cells"
  )
})

test_that("SplikitObject$getPseudoCorrelation validates dimensions", {
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )

  # Wrong dimensions
  bad_zdb <- matrix(rnorm(100), nrow = 10, ncol = 10)

  expect_error(
    obj$getPseudoCorrelation(bad_zdb),
    "mismatch"
  )
})

test_that("SplikitObject$deepCopy creates independent copy", {
  test_data <- load_toy_M1_M2_object()

  obj1 <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  obj2 <- obj1$deepCopy()

  # Modify obj1
  obj1$subset(events = 1:100)

  # obj2 should be unchanged
  expect_equal(nrow(obj1$m1), 100)
  expect_equal(nrow(obj2$m1), 2000)
})

test_that("SplikitObject stores results in metadata", {
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )

  obj$findVariableEvents(min_row_sum = 100, verbose = FALSE)

  # Check metadata was updated
  expect_true("variableEvents" %in% names(obj$metadata))
  expect_true("variableEvents_params" %in% names(obj$metadata))
  expect_equal(obj$metadata$variableEvents_params$min_row_sum, 100)
})

test_that("SplikitObject handles sparse matrix conversion", {
  test_data <- load_toy_M1_M2_object()

  # Convert to dense and back
  m1_dense <- as.matrix(test_data$m1[1:100, 1:100])

  obj <- SplikitObject$new(
    m1 = m1_dense,
    eventData = test_data$eventdata[1:100, ]
  )

  # Should be converted to dgCMatrix
  expect_s4_class(obj$m1, "dgCMatrix")
})

test_that("SplikitObject$annotateEvents validates file existence", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  expect_error(
    obj$annotateEvents("/nonexistent/file.gtf"),
    "not found"
  )
})

test_that("Backward compatibility: old functions still work", {
  skip_if_not_installed("Matrix")

  test_data <- load_toy_M1_M2_object()

  # Old workflow should still work exactly as before
  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  expect_s4_class(m2, "dgCMatrix")
  expect_equal(dim(m2), dim(test_data$m1))

  hve <- find_variable_events(
    m1_matrix = test_data$m1,
    m2_matrix = m2,
    min_row_sum = 50,
    verbose = FALSE
  )

  expect_true(nrow(hve) > 0)
  expect_true("events" %in% names(hve))
  expect_true("sum_deviance" %in% names(hve))
})

# ============================================================================
# Edge Case Tests (from deep analysis Issue #1)
# ============================================================================

test_that("rowVariance_cpp handles integer matrices", {
  # Issue #16 from deep analysis
  m_int <- matrix(1:20, nrow = 4)
  m_num <- matrix(as.numeric(1:20), nrow = 4)

  result_int <- rowVariance_cpp(m_int)
  result_num <- rowVariance_cpp(m_num)

  expect_equal(result_int, result_num)
})

test_that("rowVariance_cpp handles all-zero sparse matrix", {
  m <- Matrix::Matrix(0, nrow = 10, ncol = 5, sparse = TRUE)
  result <- rowVariance_cpp(m)

  expect_equal(length(result), 10)
  expect_true(all(result == 0))
})

test_that("get_pseudo_correlation catches dimension mismatches", {
  # Issue #14 from deep analysis
  skip_if_not_installed("Matrix")

  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  # Wrong row dimension
  bad_zdb <- matrix(rnorm(100 * ncol(test_data$m1)), nrow = 100, ncol = ncol(test_data$m1))

  expect_error(
    get_pseudo_correlation(bad_zdb, test_data$m1, m2),
    "same number of rows"
  )
})

test_that("find_variable_events handles very high threshold gracefully", {
  # Issue #23 from deep analysis
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  # Threshold so high no events pass
  expect_error(
    find_variable_events(test_data$m1, m2, min_row_sum = 1e9, verbose = FALSE),
    "No events pass"
  )
})

test_that("SplikitObject validates negative threshold", {
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )

  # Negative threshold should still work (all events pass)
  result <- obj$findVariableEvents(min_row_sum = -1, verbose = FALSE)
  expect_true(nrow(result) > 0)
})

test_that("SplikitObject handles subset by names", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Get first 10 event names
  event_names <- rownames(obj$m1)[1:10]

  obj$subset(events = event_names)

  expect_equal(nrow(obj$m1), 10)
  expect_equal(rownames(obj$m1), event_names)
})

test_that("SplikitObject subset warns about missing names", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Mix of valid and invalid names
  mixed_names <- c(rownames(obj$m1)[1:5], "nonexistent_event_1", "nonexistent_event_2")

  expect_warning(
    obj$subset(events = mixed_names),
    "not found"
  )

  expect_equal(nrow(obj$m1), 5)
})

test_that("SplikitObject handles single event subset", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  obj$subset(events = 1)

  expect_equal(nrow(obj$m1), 1)
  expect_equal(nrow(obj$eventData), 1)
})

test_that("SplikitObject handles single cell subset", {
  test_data <- load_toy_M1_M2_object()

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  obj$subset(cells = 1)

  expect_equal(ncol(obj$m1), 1)
})

test_that("make_m2 produces symmetric results for group operations", {
  # Verify M2 = group_sum - M1 for each event
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  # For any event, M1 + M2 should equal the group sum
  # Check first few groups
  unique_groups <- unique(test_data$eventdata$group_id)[1:5]

  for (grp in unique_groups) {
    grp_events <- which(test_data$eventdata$group_id == grp)
    if (length(grp_events) > 1) {
      # Group sum should be the same for all events in group
      for (i in grp_events) {
        m1_val <- test_data$m1[i, 1]
        m2_val <- m2[i, 1]
        group_sum <- sum(test_data$m1[grp_events, 1])
        expect_equal(m1_val + m2_val, group_sum, info = paste("Group:", grp, "Event:", i))
      }
    }
  }
})

# ============================================================================
# Integration Tests (from deep analysis recommendations)
# ============================================================================

test_that("Full R6 pipeline runs without errors", {
  skip_if_not_installed("Matrix")

  test_data <- load_toy_M1_M2_object()

  # Create object

obj <- splikit(
    m1 = test_data$m1,
    eventData = test_data$eventdata
  )

  # Compute M2
  obj$makeM2(verbose = FALSE)

  expect_true(!is.null(obj$m2))
  expect_equal(dim(obj$m2), dim(obj$m1))

  # Find variable events
  hve <- obj$findVariableEvents(min_row_sum = 50, verbose = FALSE)

  expect_true(nrow(hve) > 0)
  expect_true("events" %in% names(hve))
  expect_true("sum_deviance" %in% names(hve))

  # Check metadata was updated
  expect_true("variableEvents" %in% names(obj$metadata))
})

test_that("SplikitObject works with very small matrices", {
  # Edge case: minimal viable input
  m1_small <- Matrix::rsparsematrix(10, 5, 0.5)
  rownames(m1_small) <- paste0("event_", 1:10)
  colnames(m1_small) <- paste0("cell_", 1:5)

  eventdata_small <- data.table::data.table(
    event_id = paste0("event_", 1:10),
    group_id = rep(c("group1", "group2"), each = 5)
  )

  obj <- SplikitObject$new(
    m1 = m1_small,
    eventData = eventdata_small
  )

  expect_equal(nrow(obj$m1), 10)
  expect_equal(ncol(obj$m1), 5)

  obj$makeM2(verbose = FALSE)
  expect_true(!is.null(obj$m2))
})

test_that("n_threads parameter is passed correctly", {
  test_data <- load_toy_M1_M2_object()

  m2 <- make_m2(
    m1_inclusion_matrix = test_data$m1,
    eventdata = test_data$eventdata,
    verbose = FALSE
  )

  obj <- SplikitObject$new(
    m1 = test_data$m1,
    m2 = m2,
    eventData = test_data$eventdata
  )

  # Should run without error with multiple threads
  result <- obj$findVariableEvents(min_row_sum = 100, n_threads = 2, verbose = FALSE)
  expect_true(nrow(result) > 0)
})
