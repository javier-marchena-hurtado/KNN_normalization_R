# Check that the fast normalize_with_neighbors gives the same results as the original (slow)
# algorithm, and the same results as the Python package for the same KNN graph.

calculate_size_factors <- function(x, target_cells, neighbor_cells, neighbors_per_cell) {
  # Original implementation: median ratios of all (cell, neighbor) pairs, recomputed from the data.
  n_cells <- nrow(x)
  median_ratios <- vapply(seq_along(target_cells), function(i) {
    stats::median(x[neighbor_cells[i], ] / x[target_cells[i], ])
  }, numeric(1))
  ratio_sums <- as.numeric(tapply(median_ratios, factor(target_cells, levels = seq_len(n_cells)), sum))
  ratio_sums[is.na(ratio_sums)] <- 0
  ifelse(neighbors_per_cell > 0, ratio_sums / pmax(neighbors_per_cell, 1), 1)
}

normalize_reference <- function(counts, connectivities, log_transform, pseudocount = 5,
                                max_iterations = 25, change_for_stop = 0.0005) {
  # Original implementation: recomputes the median ratios in each iteration.
  x <- t(counts) + pseudocount
  total_counts_before <- sum(x)
  history <- list()
  pairs <- KNN.normalization:::.neighbor_pairs(connectivities)
  for (iteration in seq_len(max_iterations)) {
    s <- calculate_size_factors(x, pairs$target_cells, pairs$neighbor_cells, pairs$neighbors_per_cell)
    s <- s / (sum(s * rowSums(x)) / total_counts_before)
    x <- x * s
    history[[iteration]] <- s
    if (iteration > 1 && max(abs(history[[iteration]] - history[[iteration - 1]])) < change_for_stop) break
  }
  if (log_transform) x <- log(x)
  list(x = x, history = do.call(cbind, history))
}

isolate_first_cell <- function(graph) {
  graph[1, ] <- 0
  Matrix::drop0(graph)
}

for (log_transform in c(TRUE, FALSE)) {
  for (change_for_stop in c(0.0005, 0)) { # 0 never stops early: runs all 25 iterations
    for (isolated_cell in c(FALSE, TRUE)) {
      test_that(sprintf(
        "fast normalization matches the reference (log=%s, change_for_stop=%s, isolated=%s)",
        log_transform, change_for_stop, isolated_cell
      ), {
        object <- make_object()
        connectivities <- make_graph(k = 4)
        if (isolated_cell) connectivities <- isolate_first_cell(connectivities)
        counts <- as.matrix(SeuratObject::LayerData(object, assay = "ADT", layer = "counts"))
        expected <- normalize_reference(
          counts, connectivities,
          log_transform = log_transform, change_for_stop = change_for_stop
        )
        result <- normalize_with_neighbors(
          object, connectivities,
          log_transform = log_transform, change_for_stop = change_for_stop,
          save_size_factors = TRUE, verbose = FALSE
        )
        expect_equal(t(as.matrix(normalized_data(result))), expected$x, tolerance = 1e-7, ignore_attr = TRUE)
        history <- SeuratObject::Misc(result, "knn_normalization")$size_factor_history
        expect_equal(unname(history), unname(expected$history), tolerance = 1e-7)
        expect_equal(unname(result$size_factor), unname(apply(expected$history, 1, prod)), tolerance = 1e-7)
      })
    }
  }
}

test_that("results are identical to those of the Python package for the same KNN graph", {
  fixtures <- testthat::test_path("fixtures")
  read_csv <- function(name, header = FALSE) utils::read.csv(file.path(fixtures, name), header = header)
  counts <- t(as.matrix(read_csv("counts.csv")))
  dimnames(counts) <- list(paste0("protein", seq_len(nrow(counts))), paste0("cell", seq_len(ncol(counts))))
  graph_table <- read_csv("graph.csv", header = TRUE)
  graph <- Matrix::sparseMatrix(
    i = graph_table$i, j = graph_table$j, x = graph_table$x, dims = rep(ncol(counts), 2)
  )
  object <- SeuratObject::CreateSeuratObject(methods::as(counts, "CsparseMatrix"), assay = "ADT")

  for (tag in c("log", "nolog")) {
    result <- normalize_with_neighbors(
      object, graph,
      log_transform = (tag == "log"), save_size_factors = TRUE, verbose = FALSE
    )
    expect_equal(
      unname(t(as.matrix(normalized_data(result)))),
      unname(as.matrix(read_csv(sprintf("expected_%s.csv", tag)))),
      tolerance = 1e-8
    )
    expect_equal(
      unname(SeuratObject::Misc(result, "knn_normalization")$size_factor_history),
      unname(as.matrix(read_csv(sprintf("expected_history_%s.csv", tag)))),
      tolerance = 1e-8
    )
    expect_equal(
      unname(result$size_factor),
      read_csv(sprintf("expected_size_factor_%s.csv", tag))[[1]],
      tolerance = 1e-8
    )
  }
})
