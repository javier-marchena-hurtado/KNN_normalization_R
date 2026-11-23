test_that("invalid arguments raise clear errors", {
  object <- make_object()
  expect_error(knn_normalize(1:3), "must be a Seurat object")
  expect_error(knn_normalize(object, "somewhere"), "should be one of")
  expect_error(knn_normalize(object, "prot", max_iterations = 0), "max_iterations must be a positive integer")
  expect_error(knn_normalize(object, "prot", n_neighbors = 1.5), "n_neighbors must be a positive integer")
  expect_error(knn_normalize(object, "prot", n_neighbors = 100), "must be smaller than the number of cells")
  expect_error(knn_normalize(object, "prot", pseudocount = -1), "pseudocount must be a finite number")
  expect_error(knn_normalize(object, "prot", log_transform = NA), "log_transform must be TRUE or FALSE")
})

test_that("missing assays and graphs raise clear errors", {
  object <- make_object()
  expect_error(knn_normalize(object, "prot", protein_assay = "nope"), "no protein assay called 'nope'")
  expect_error(knn_normalize(object, "rna", n_neighbors = 3), "requires an RNA assay")
  expect_error(knn_normalize(object, "use_existing_neighbors"), "requires graph_name")
  expect_error(knn_normalize(object, "use_existing_neighbors", graph_name = "nope"), "no graph called")
  object[["small"]] <- SeuratObject::as.Graph(make_graph(n_cells = 50))
  expect_error(
    knn_normalize(object, "use_existing_neighbors", graph_name = "small"),
    "has dimensions 50 x 50"
  )
  expect_error(normalize_with_neighbors(object, make_graph(n_cells = 50)), "has dimensions 50 x 50")
})

test_that("too few cells raise an error", {
  few <- SeuratObject::CreateSeuratObject(simulate_counts(20, 10, n_cells = 10), assay = "ADT")
  expect_error(knn_normalize(few, "prot"), "Too few cells")
})

test_that("weird data raises errors and warnings", {
  counts <- simulate_counts(20, 10)
  non_integer <- SeuratObject::CreateSeuratObject(methods::as(counts + 0.5, "CsparseMatrix"), assay = "ADT")
  expect_warning(
    knn_normalize(non_integer, "prot", n_neighbors = 3, verbose = FALSE),
    "not integer-valued"
  )
  negative <- counts
  negative[1, 1] <- -1
  expect_warning(
    expect_warning(
      knn_normalize(SeuratObject::CreateSeuratObject(methods::as(negative, "CsparseMatrix"), assay = "ADT"), "prot", n_neighbors = 3, verbose = FALSE),
      "negative values"
    ),
    NA
  )
  nan <- counts
  nan[1, 1] <- NA
  expect_error(
    knn_normalize(SeuratObject::CreateSeuratObject(methods::as(nan, "CsparseMatrix"), assay = "ADT"), "prot", n_neighbors = 3),
    "NaN or infinite"
  )
})

test_that("normalizing twice and overwriting size factors warns", {
  object <- make_object()
  once <- knn_normalize(object, "prot", n_neighbors = 3, save_size_factors = TRUE, verbose = FALSE)
  twice <- function() knn_normalize(once, "prot", n_neighbors = 3, save_size_factors = TRUE, verbose = FALSE)
  expect_warning(
    expect_warning(twice(), "already been KNN-normalized"),
    "size_factor already exists"
  )
})

test_that("the default number of neighbors scales with the number of cells", {
  default_n <- KNN.normalization:::.default_n_neighbors
  expect_equal(default_n(100), 15)
  expect_equal(default_n(1000), 50)
  expect_equal(default_n(4075), 204)
  expect_equal(default_n(1e6), 300)
})
