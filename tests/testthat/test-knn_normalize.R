test_that("normalize_with_neighbors changes the data", {
  object <- make_object()
  result <- normalize_with_neighbors(object, make_graph(), log_transform = FALSE, verbose = FALSE)
  counts <- SeuratObject::LayerData(object, assay = "ADT", layer = "counts")
  expect_false(isTRUE(all.equal(as.matrix(normalized_data(result)), as.matrix(counts), check.attributes = FALSE)))
  expect_equal(dim(normalized_data(result)), dim(counts))
  expect_identical(
    as.matrix(SeuratObject::LayerData(result, assay = "ADT", layer = "counts")),
    as.matrix(counts)
  ) # counts are untouched
})

test_that(".get_neighbors returns a square graph for each neighbor source", {
  object <- make_object(with_rna = TRUE)
  object[["knn"]] <- SeuratObject::as.Graph(make_graph())
  graphs <- list(
    KNN.normalization:::.get_neighbors(object, "prot", n_neighbors = 3, verbose = FALSE),
    KNN.normalization:::.get_neighbors(object, "rna", n_neighbors = 3, verbose = FALSE),
    KNN.normalization:::.get_neighbors(object, "use_existing_neighbors", n_neighbors = 3, graph_name = "knn")
  )
  for (graph in graphs) {
    expect_s4_class(graph, "dgCMatrix")
    expect_equal(dim(graph), c(100, 100))
    expect_true(all(Matrix::diag(graph) == 0))
  }
  # the graph computed from the data is symmetric
  expect_true(Matrix::isSymmetric(graphs[[1]]))
  expect_true(Matrix::isSymmetric(graphs[[2]]))
})

test_that("knn_normalize works with protein neighbors", {
  object <- make_object()
  result <- knn_normalize(object, calculate_neighbors_from = "prot", n_neighbors = 3, verbose = FALSE)
  expect_s4_class(result, "Seurat")
  expect_false(is.null(SeuratObject::Misc(result, "knn_normalization")))
  # the input is untouched (R copies on modify)
  expect_null(SeuratObject::Misc(object, "knn_normalization"))
  expect_true(all(is.finite(normalized_data(result))))
})

test_that("knn_normalize works with RNA neighbors and saves the size factors", {
  object <- make_object(with_rna = TRUE)
  result <- knn_normalize(
    object,
    calculate_neighbors_from = "rna", n_neighbors = 3, save_size_factors = TRUE, verbose = FALSE
  )
  expect_s4_class(result, "Seurat")
  expect_length(result$size_factor, 100)
  history <- SeuratObject::Misc(result, "knn_normalization")$size_factor_history
  expect_equal(nrow(history), 100)
  expect_equal(unname(result$size_factor), unname(apply(history, 1, prod)))
  expect_true("connectivities_KNN_normalization" %in% names(result@graphs))
})

test_that("knn_normalize works with an existing graph and prints progress", {
  object <- make_object()
  object[["knn"]] <- SeuratObject::as.Graph(make_graph())
  expect_message(
    result <- knn_normalize(object, "use_existing_neighbors", graph_name = "knn"),
    "Running KNN normalization"
  )
  expected <- normalize_with_neighbors(object, make_graph(), verbose = FALSE)
  expect_equal(normalized_data(result), normalized_data(expected))
})

test_that("the total counts are conserved by the size factors", {
  object <- make_object()
  result <- normalize_with_neighbors(object, make_graph(), log_transform = FALSE, verbose = FALSE)
  counts <- SeuratObject::LayerData(object, assay = "ADT", layer = "counts")
  expect_equal(sum(normalized_data(result)), sum(counts + 5))
})
