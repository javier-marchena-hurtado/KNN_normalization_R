# Test data: 100 cells, 20 proteins (and 100 genes). Same simulation as in the Python tests.
simulate_counts <- function(n_features, lambda, n_cells = 100, seed = 0) {
  set.seed(seed)
  methods::as(matrix(
    rpois(n_features * n_cells, lambda),
    nrow = n_features,
    dimnames = list(paste0("feature", seq_len(n_features)), paste0("cell", seq_len(n_cells)))
  ), "CsparseMatrix")
}

make_object <- function(with_rna = FALSE) {
  object <- SeuratObject::CreateSeuratObject(simulate_counts(20, 10), assay = "ADT")
  if (with_rna) {
    object[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = simulate_counts(100, 5, seed = 1))
  }
  object
}

# A KNN graph in which every cell is a neighbor of the next k cells.
make_graph <- function(n_cells = 100, k = 5) {
  cells <- paste0("cell", seq_len(n_cells))
  Matrix::sparseMatrix(
    i = rep(seq_len(n_cells), each = k),
    j = as.vector(vapply(seq_len(n_cells), function(i) (i + seq_len(k) - 1) %% n_cells + 1, numeric(k))),
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(cells, cells)
  )
}

normalized_data <- function(object, assay = "ADT") {
  SeuratObject::LayerData(object, assay = assay, layer = "data")
}
