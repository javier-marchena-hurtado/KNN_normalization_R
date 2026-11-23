#' Normalize protein expression with KNN normalization
#'
#' KNN normalization is a normalization method for protein counts in CITE-seq data. It
#' learns from neighbor cells in a KNN graph in order to estimate the appropriate total
#' protein counts in each cell, accurately estimating total protein counts while
#' preserving biological information.
#'
#' @param object A [Seurat][SeuratObject::Seurat-class] object with a protein assay
#'   (and, for `calculate_neighbors_from = "rna"`, an RNA assay). The raw protein counts
#'   must be in the `counts` layer of the protein assay. The function expects count data
#'   (integers).
#' @param calculate_neighbors_from Whether to calculate neighbors from protein or RNA
#'   data. `"rna"` requires an RNA assay in `object`. The RNA counts in its `counts` layer
#'   must be **raw counts**: they are library-size normalized and log-transformed
#'   internally. `"prot"` only needs the protein assay. `"use_existing_neighbors"` uses
#'   the graph `object@graphs[[graph_name]]`, which must exist and be a square graph over
#'   all cells.
#' @param log_transform Whether to log-transform the protein data.
#' @param n_neighbors Number of neighbors. If `NULL`, calculated automatically as
#'   `max(15, min(round(n_cells / 20), 300))`. As in scanpy, a cell counts as one of its own
#'   nearest neighbors in this number, and the graph is symmetrized.
#' @param pseudocount Pseudocount to add before normalization to avoid zero-division
#'   errors. If you have no zeroes in the protein data, you don't need to use pseudocounts;
#'   in that case pass `pseudocount = 0`.
#' @param max_iterations Maximum number of iterations.
#' @param save_size_factors If `TRUE`, saves the final size factors to the metadata column
#'   `size_factor` of `object`, the size factor history (cells x iterations) to
#'   `Misc(object, "knn_normalization")$size_factor_history`, and the KNN graph to
#'   `object[["connectivities_KNN_normalization"]]`.
#' @param verbose Whether to print progress messages.
#' @param protein_assay Name of the assay with the protein data.
#' @param rna_assay Name of the assay with the RNA data. Only used if
#'   `calculate_neighbors_from = "rna"`.
#' @param graph_name Name of the graph in `object@graphs` to use if
#'   `calculate_neighbors_from = "use_existing_neighbors"`, for example `"RNA_nn"`. A graph
#'   is a cells x cells matrix, the nonzero entries of which mark the neighbors of each cell
#'   (rows). The diagonal is ignored.
#' @param n_cores Number of cores used to compute the median ratios, the slowest step of
#'   the algorithm. Defaults to the option `mc.cores`, or 1 if it is not set. Ignored on
#'   Windows.
#'
#' @return A copy of `object` in which the `data` layer of the protein assay contains the
#'   KNN-normalized (and, by default, log-transformed) protein expression. The `counts`
#'   layer is left untouched. The parameters used are saved in
#'   `Misc(object, "knn_normalization")`.
#'
#' @note Raw counts must be stored in the `counts` layer of the protein assay. In R objects
#'   are not modified in place: assign the result, e.g. `object <- knn_normalize(object)`.
#'
#' @seealso [normalize_with_neighbors()] to normalize with a precomputed KNN graph.
#'
#' @examples
#' set.seed(1)
#' n_cells <- 100
#' counts <- matrix(rpois(n_cells * 20, lambda = 10), nrow = 20,
#'                  dimnames = list(paste0("protein", 1:20), paste0("cell", 1:n_cells)))
#' object <- SeuratObject::CreateSeuratObject(counts, assay = "ADT")
#' object <- knn_normalize(object, calculate_neighbors_from = "prot",
#'                         protein_assay = "ADT", n_neighbors = 5, verbose = FALSE)
#' SeuratObject::LayerData(object, assay = "ADT", layer = "data")[1:3, 1:3]
#' @export
knn_normalize <- function(object,
                          calculate_neighbors_from = c("rna", "prot", "use_existing_neighbors"),
                          log_transform = TRUE,
                          n_neighbors = NULL,
                          pseudocount = 5,
                          max_iterations = 25,
                          save_size_factors = FALSE,
                          verbose = TRUE,
                          protein_assay = "ADT",
                          rna_assay = "RNA",
                          graph_name = NULL,
                          n_cores = getOption("mc.cores", 1L)) {
  calculate_neighbors_from <- match.arg(calculate_neighbors_from)
  # Raise an error (or warn) if the arguments or the data cannot be normalized.
  .validate_inputs(
    object,
    calculate_neighbors_from,
    log_transform,
    n_neighbors,
    pseudocount,
    max_iterations,
    save_size_factors,
    protein_assay,
    rna_assay,
    graph_name
  )
  if (is.null(n_neighbors)) n_neighbors <- .default_n_neighbors(ncol(object))
  # KNN graph (sparse matrix of connectivities) that tells which cells are neighbors.
  neighbors <- .get_neighbors(
    object,
    calculate_neighbors_from,
    n_neighbors,
    protein_assay = protein_assay,
    rna_assay = rna_assay,
    graph_name = graph_name,
    log_transform = log_transform,
    verbose = verbose
  )
  # Normalize the protein data using the neighbors. This writes the result to the data layer.
  normalize_with_neighbors(
    object,
    neighbors,
    protein_assay = protein_assay,
    log_transform = log_transform,
    pseudocount = pseudocount,
    max_iterations = max_iterations,
    verbose = verbose,
    save_size_factors = save_size_factors,
    n_cores = n_cores
  )
}

#' Apply KNN normalization given precomputed neighbors
#'
#' @param object A [Seurat][SeuratObject::Seurat-class] object with protein expression
#'   data in the `counts` layer of the assay `protein_assay`.
#' @param connectivities KNN graph containing neighbor cells: a square (cells x cells)
#'   matrix, e.g. a `dgCMatrix` or a `Graph`, in which row `i` has nonzero entries at the
#'   neighbors of cell `i`. The weights are ignored. The cells must be in the same order
#'   as in `object`.
#' @param protein_assay Name of the assay with the protein data.
#' @param log_transform If `TRUE`, log-transforms the data after normalization.
#' @param save_size_factors If `TRUE`, saves the final size factors to `object$size_factor`
#'   and the size factor history to `Misc(object, "knn_normalization")$size_factor_history`.
#' @param pseudocount Pseudocount added to avoid zero-division errors.
#' @param max_iterations Maximum number of iterations.
#' @param change_for_stop Convergence criterion. The algorithm stops when the change in
#'   size factor between iterations is smaller than this value.
#' @param verbose Whether to print progress messages.
#' @param n_cores Number of cores used to compute the median ratios. Defaults to the
#'   option `mc.cores`, or 1 if it is not set. Ignored on Windows.
#'
#' @return A copy of `object` with the normalized protein data in the `data` layer of
#'   `protein_assay`.
#'
#' @examples
#' set.seed(1)
#' n_cells <- 100
#' counts <- matrix(rpois(n_cells * 20, lambda = 10), nrow = 20,
#'                  dimnames = list(paste0("protein", 1:20), paste0("cell", 1:n_cells)))
#' object <- SeuratObject::CreateSeuratObject(counts, assay = "ADT")
#' # A toy KNN graph in which every cell is a neighbor of the next 5 cells
#' graph <- Matrix::sparseMatrix(
#'   i = rep(1:n_cells, each = 5),
#'   j = as.vector(sapply(1:n_cells, function(i) (i + 0:4) %% n_cells + 1)),
#'   x = 1, dims = c(n_cells, n_cells)
#' )
#' object <- normalize_with_neighbors(object, graph, protein_assay = "ADT", verbose = FALSE)
#' @export
normalize_with_neighbors <- function(object,
                                     connectivities,
                                     protein_assay = "ADT",
                                     log_transform = TRUE,
                                     save_size_factors = FALSE,
                                     pseudocount = 5,
                                     max_iterations = 25,
                                     change_for_stop = 0.0005,
                                     verbose = TRUE,
                                     n_cores = getOption("mc.cores", 1L)) {
  if (!inherits(object, "Seurat")) {
    stop(sprintf("object must be a Seurat object, got %s.", class(object)[1]), call. = FALSE)
  }
  if (!isTRUE(all(dim(connectivities) == ncol(object)))) {
    stop(sprintf(
      "connectivities has dimensions %s, expected %d x %d.",
      paste(dim(connectivities), collapse = " x "), ncol(object), ncol(object)
    ), call. = FALSE)
  }
  # Dense float copy of the protein data (cells x proteins), so that the original data is
  # not modified until the end.
  x <- .densify(.get_counts(object, protein_assay))
  x <- x + pseudocount
  total_counts <- sum(x)
  row_sums <- rowSums(x)
  # Every (cell, neighbor) pair of the KNN graph: pair i is (target_cells[i], neighbor_cells[i]).
  pairs <- .neighbor_pairs(connectivities)
  target_cells <- pairs$target_cells
  neighbor_cells <- pairs$neighbor_cells
  neighbors_per_cell <- pairs$neighbors_per_cell
  size_factor_history <- list()
  cumulative_factors <- rep(1, nrow(x)) # Product of the size factors of all iterations until now.

  if (verbose) message("Running KNN normalization")

  # Median over proteins of neighbor / cell expression for every pair.
  median_ratios <- .pairwise_median_ratios(x, target_cells, neighbor_cells, n_cores = n_cores)

  for (iteration in seq_len(max_iterations)) {
    # Size factor of a cell = average of its median ratios over its neighbors.
    size_factors <- .size_factors_from_ratios(median_ratios, target_cells, neighbors_per_cell)
    # Rescale so that the total counts in the matrix stay constant.
    total_after_scaling <- sum(size_factors * (row_sums * cumulative_factors))
    size_factors <- size_factors * (total_counts / total_after_scaling)
    size_factor_history[[iteration]] <- size_factors
    cumulative_factors <- cumulative_factors * size_factors

    if (verbose) message("Iteration ", iteration)
    if (iteration > 1) {
      change <- max(abs(size_factors - size_factor_history[[iteration - 1]]))
      if (verbose) message("Change wrt previous iteration: ", change)
      if (change < change_for_stop) break
    }

    # Update the median ratios for the next iteration: neighbor gets scaled by s_n, cell by s_c.
    median_ratios <- median_ratios * (size_factors[neighbor_cells] / size_factors[target_cells])
  }

  # Apply the size factors of all iterations at once.
  x <- x * cumulative_factors
  if (log_transform) x <- log(x)

  # Seurat stores features x cells.
  SeuratObject::LayerData(object, assay = protein_assay, layer = "data") <- t(x)
  info <- list(pseudocount = pseudocount, log_transform = log_transform)
  if (save_size_factors) {
    info$size_factor_history <- .size_factor_history_matrix(size_factor_history, colnames(object))
    object <- .save_size_factors(object, info$size_factor_history, connectivities, protein_assay)
  }
  # Normalizing a second time already warned in knn_normalize(); do not warn about overwriting.
  suppressWarnings(SeuratObject::Misc(object, slot = "knn_normalization") <- info)
  object
}
