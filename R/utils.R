# Internal helpers of KNN.normalization.
#
# Internally, matrices are cells x proteins (as in the Python package), while Seurat
# stores features x cells. Conversion happens only when reading from / writing to
# the Seurat object.

# Below this number of proteins, neighbors are computed directly on the protein
# data; at or above it, PCA is applied first.
PCA_PROTEIN_THRESHOLD <- 70L

NEIGHBOR_SOURCES <- c("prot", "rna", "use_existing_neighbors")
N_CELLS_TO_CHECK <- 100L # only the first cells are inspected when validating the values

# Number of principal components used to compute the KNN graph.
N_PCS <- 50L

.check_positive_int <- function(value, name) {
  ok <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= 1 && value == round(value)
  if (!ok) {
    stop(sprintf("%s must be a positive integer, got %s.", name, format(value)), call. = FALSE)
  }
}

.check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be TRUE or FALSE, got %s.", name, format(value)), call. = FALSE)
  }
}

# Check the first cells of `x` (features x cells): stop on NaN or infinite values,
# warn on negative values. Returns TRUE if all checked values are integers.
.check_values <- function(x, name) {
  head_x <- x[, seq_len(min(ncol(x), N_CELLS_TO_CHECK)), drop = FALSE]
  # sparse matrices only store their non-zeros
  values <- if (methods::is(head_x, "sparseMatrix")) {
    methods::as(head_x, "CsparseMatrix")@x
  } else {
    as.numeric(as.matrix(head_x))
  }
  if (!all(is.finite(values))) {
    stop(sprintf("The %s data contains NaN or infinite values.", name), call. = FALSE)
  }
  if (length(values) && min(values) < 0) {
    warning(sprintf(
      "The %s data contains negative values; raw (non-negative) counts are expected.", name
    ), call. = FALSE)
  }
  all(values == round(values))
}

# Raw counts of an assay as a features x cells matrix.
.get_counts <- function(object, assay) {
  SeuratObject::LayerData(object, assay = assay, layer = "counts")
}

.validate_inputs <- function(object,
                             calculate_neighbors_from,
                             log_transform,
                             n_neighbors,
                             pseudocount,
                             max_iterations,
                             save_size_factors,
                             protein_assay,
                             rna_assay,
                             graph_name) {
  # Check arguments of knn_normalize()
  if (!inherits(object, "Seurat")) {
    stop(sprintf("object must be a Seurat object, got %s.", class(object)[1]), call. = FALSE)
  }
  if (!is.character(calculate_neighbors_from) || length(calculate_neighbors_from) != 1L ||
    !calculate_neighbors_from %in% NEIGHBOR_SOURCES) {
    stop(sprintf(
      "calculate_neighbors_from must be one of %s, got %s.",
      paste(shQuote(NEIGHBOR_SOURCES, type = "cmd"), collapse = ", "),
      format(calculate_neighbors_from)
    ), call. = FALSE)
  }
  .check_flag(log_transform, "log_transform")
  .check_flag(save_size_factors, "save_size_factors")
  .check_positive_int(max_iterations, "max_iterations")
  if (!is.null(n_neighbors)) .check_positive_int(n_neighbors, "n_neighbors")
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L || is.na(pseudocount) ||
    !is.finite(pseudocount) || pseudocount < 0) {
    stop(sprintf("pseudocount must be a finite number >= 0, got %s.", format(pseudocount)), call. = FALSE)
  }

  # Check the structure of the input
  assays <- SeuratObject::Assays(object)
  if (!protein_assay %in% assays) {
    stop(sprintf(
      "The Seurat object has no protein assay called '%s'. Use the protein_assay argument to give its name.",
      protein_assay
    ), call. = FALSE)
  }
  if (calculate_neighbors_from == "rna") {
    if (!rna_assay %in% assays) {
      stop(sprintf(
        "calculate_neighbors_from = 'rna' requires an RNA assay; the object has no assay called '%s'. Use the rna_assay argument to give its name.",
        rna_assay
      ), call. = FALSE)
    }
    if (identical(rna_assay, protein_assay)) {
      stop("rna_assay and protein_assay must be different assays.", call. = FALSE)
    }
  }
  if (!"counts" %in% SeuratObject::Layers(object[[protein_assay]])) {
    stop(sprintf("The '%s' assay has no 'counts' layer.", protein_assay), call. = FALSE)
  }
  n_cells <- ncol(object)
  if (!is.null(n_neighbors) && n_neighbors >= n_cells) {
    stop(sprintf(
      "n_neighbors (%d) must be smaller than the number of cells (%d).", as.integer(n_neighbors), n_cells
    ), call. = FALSE)
  }
  if (n_cells <= .default_n_neighbors(n_cells)) {
    stop(sprintf("Too few cells (%d) to compute a KNN graph.", n_cells), call. = FALSE)
  }

  # Check the neighbor source
  if (calculate_neighbors_from == "use_existing_neighbors") {
    if (is.null(graph_name)) {
      stop(
        "calculate_neighbors_from = 'use_existing_neighbors' requires graph_name, ",
        "the name of a graph stored in the object (see names(object@graphs)). ",
        "Run Seurat::FindNeighbors() on the protein data first.",
        call. = FALSE
      )
    }
    if (!is.character(graph_name) || length(graph_name) != 1L || !graph_name %in% names(object@graphs)) {
      stop(sprintf("The object has no graph called %s.", format(graph_name)), call. = FALSE)
    }
    graph_dim <- dim(object@graphs[[graph_name]])
    if (!isTRUE(all(graph_dim == n_cells))) {
      stop(sprintf(
        "The graph '%s' has dimensions %d x %d, expected %d x %d.",
        graph_name, graph_dim[1], graph_dim[2], n_cells, n_cells
      ), call. = FALSE)
    }
  }

  # Check protein values (first cells only)
  if (!.check_values(.get_counts(object, protein_assay), "protein")) {
    warning("The protein data is not integer-valued. KNN normalization expects raw counts.", call. = FALSE)
  }
  if (!is.null(SeuratObject::Misc(object, slot = "knn_normalization"))) {
    warning(
      "This object has already been KNN-normalized; normalizing again is possibly a mistake.",
      call. = FALSE
    )
  }
  if (save_size_factors && "size_factor" %in% colnames(object[[]])) {
    warning("object$size_factor already exists and will be overwritten.", call. = FALSE)
  }

  # Check RNA values (first cells only)
  if (calculate_neighbors_from == "rna" && !.check_values(.get_counts(object, rna_assay), "RNA")) {
    warning(
      "The RNA data is not integer-valued. KNN normalization expects raw RNA counts ",
      "and will library-size normalize and log-transform them.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Dense double matrix, cells x proteins, from a features x cells matrix.
.densify <- function(x) {
  out <- t(as.matrix(x))
  storage.mode(out) <- "double"
  out
}

.default_n_neighbors <- function(n_cells) {
  max(15, min(round(n_cells / 20), 300))
}

# Principal components (scores) of a cells x features matrix after centering the features.
.pca <- function(x, n_comps = N_PCS) {
  n_comps <- min(n_comps, nrow(x) - 1L, ncol(x) - 1L)
  centers <- Matrix::colMeans(x)
  if (n_comps >= min(dim(x)) / 2) {
    # irlba is only worthwhile when few components are needed: use an exact SVD instead.
    centered <- sweep(as.matrix(x), 2L, centers)
    decomposition <- svd(centered, nu = n_comps, nv = 0L)
    scores <- sweep(decomposition$u[, seq_len(n_comps), drop = FALSE], 2L, decomposition$d[seq_len(n_comps)], "*")
  } else {
    decomposition <- irlba::irlba(x, nv = n_comps, center = centers)
    scores <- sweep(decomposition$u, 2L, decomposition$d, "*")
  }
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC_", seq_len(ncol(scores)))
  scores
}

# Drop the diagonal and make the graph symmetric (as the connectivities of scanpy are):
# two cells are neighbors if one of them is among the nearest neighbors of the other.
.symmetrize_graph <- function(graph) {
  graph <- methods::as(methods::as(graph, "CsparseMatrix"), "generalMatrix")
  graph <- Matrix::drop0(graph)
  pattern <- methods::as(graph != 0, "CsparseMatrix")
  pattern <- pattern | Matrix::t(pattern)
  pattern <- methods::as(pattern, "dMatrix")
  Matrix::diag(pattern) <- 0
  Matrix::drop0(methods::as(pattern, "CsparseMatrix"))
}

# Get the KNN graph (sparse matrix of connectivities) from the protein or RNA data, or reuse
# the existing one.
#
# For protein data, PCA is applied first if there are at least PCA_PROTEIN_THRESHOLD proteins,
# and the cosine distance is used. RNA data is library-size normalized and log-transformed,
# and always uses PCA. The graph is computed with Seurat::FindNeighbors().
.get_neighbors <- function(object,
                           calculate_neighbors_from,
                           n_neighbors,
                           protein_assay = "ADT",
                           rna_assay = "RNA",
                           graph_name = NULL,
                           log_transform = TRUE,
                           verbose = TRUE) {
  if (calculate_neighbors_from == "use_existing_neighbors") {
    return(.drop_diagonal(object@graphs[[graph_name]]))
  }

  if (calculate_neighbors_from == "prot") {
    if (verbose) message("Calculating KNN graph from protein data")
    x <- .densify(.get_counts(object, protein_assay))
    if (log_transform) x <- log1p(x)
    use_pca <- ncol(x) >= PCA_PROTEIN_THRESHOLD
    metric <- "cosine"
  } else {
    if (verbose) message("Calculating KNN graph from RNA data")
    counts <- methods::as(.get_counts(object, rna_assay), "CsparseMatrix")
    x <- Matrix::t(.log_normalize_total(counts))
    use_pca <- TRUE
    metric <- "euclidean"
  }

  if (use_pca) x <- .pca(x)
  rownames(x) <- colnames(object)
  if (is.null(colnames(x))) colnames(x) <- paste0("dim_", seq_len(ncol(x)))
  neighbors <- Seurat::FindNeighbors(
    x,
    k.param = n_neighbors,
    compute.SNN = FALSE,
    nn.method = "annoy",
    annoy.metric = metric,
    verbose = FALSE
  )
  graph <- if (is.list(neighbors)) neighbors$nn else neighbors
  .symmetrize_graph(graph)
}

# An existing graph is used as it is, except that its diagonal (a cell is not its own neighbor) is dropped.
.drop_diagonal <- function(graph) {
  graph <- methods::as(methods::as(graph, "CsparseMatrix"), "generalMatrix")
  Matrix::diag(graph) <- 0
  Matrix::drop0(graph)
}

# Library-size normalize (to the median library size, like scanpy::normalize_total) and
# log1p-transform a genes x cells count matrix.
.log_normalize_total <- function(counts) {
  library_sizes <- Matrix::colSums(counts)
  target_sum <- stats::median(library_sizes[library_sizes > 0])
  scale_factors <- ifelse(library_sizes > 0, target_sum / library_sizes, 0)
  normalized <- counts %*% Matrix::Diagonal(x = scale_factors)
  normalized <- methods::as(normalized, "CsparseMatrix")
  normalized@x <- log1p(normalized@x)
  dimnames(normalized) <- dimnames(counts)
  normalized
}

# Return list(target_cells, neighbor_cells, neighbors_per_cell) of the KNN graph (1-based).
#
# target_cells[i] and neighbor_cells[i] form the i-th (cell, neighbor) pair.
.neighbor_pairs <- function(connectivities) {
  graph <- Matrix::drop0(methods::as(methods::as(connectivities, "CsparseMatrix"), "generalMatrix"))
  triplets <- Matrix::summary(graph)
  ordering <- order(triplets$i, triplets$j)
  target_cells <- triplets$i[ordering]
  neighbor_cells <- triplets$j[ordering]
  list(
    target_cells = target_cells,
    neighbor_cells = neighbor_cells,
    neighbors_per_cell = tabulate(target_cells, nbins = nrow(graph))
  )
}

# Median over proteins of neighbor / cell expression, for every (cell, neighbor) pair.
#
# This is the expensive step of the normalization. It only has to be computed once:
# afterwards the ratios can be updated from the size factors alone (see
# normalize_with_neighbors()).
#
# The pairs are split into chunks. The chunks do not depend on each other, so several chunks
# can be processed at the same time on several cores (see n_cores).
#
# x: dense matrix of protein data (cells x proteins, pseudocount already added).
# chunk_size: number of (cell, neighbor) pairs in one chunk. Every core holds one chunk in
#   memory at a time, so a smaller chunk size means less memory use.
# Returns a numeric vector with one median ratio per (cell, neighbor) pair.
.pairwise_median_ratios <- function(x, target_cells, neighbor_cells, chunk_size = 20000L, n_cores = 1L) {
  n_pairs <- length(target_cells)
  first_pair_of_each_chunk <- seq.int(1L, n_pairs, by = chunk_size)

  median_ratios_of_chunk <- function(first_pair) {
    chunk <- first_pair:min(first_pair + chunk_size - 1L, n_pairs)
    ratios <- x[neighbor_cells[chunk], , drop = FALSE] / x[target_cells[chunk], , drop = FALSE]
    matrixStats::rowMedians(ratios)
  }

  n_cores <- as.integer(min(n_cores, length(first_pair_of_each_chunk)))
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    chunks <- parallel::mclapply(first_pair_of_each_chunk, median_ratios_of_chunk, mc.cores = n_cores)
    failed <- vapply(chunks, inherits, logical(1), what = "try-error")
    if (any(failed)) stop(as.character(chunks[[which(failed)[1]]]), call. = FALSE)
  } else {
    chunks <- lapply(first_pair_of_each_chunk, median_ratios_of_chunk)
  }
  unlist(chunks, use.names = FALSE)
}

# Size factor of each cell: the average of the median ratios of all of its neighbors. A cell
# without neighbors gets a neutral size factor of 1.
.size_factors_from_ratios <- function(median_ratios, target_cells, neighbors_per_cell) {
  n_cells <- length(neighbors_per_cell)
  has_neighbors <- neighbors_per_cell > 0

  # Add up the median ratios of all pairs that belong to the same cell.
  ratio_sum_per_cell <- numeric(n_cells)
  sums <- rowsum(median_ratios, target_cells, reorder = TRUE)
  ratio_sum_per_cell[as.integer(rownames(sums))] <- sums[, 1]

  # Average them. Cells without neighbors keep the neutral value 1.
  size_factors <- rep(1, n_cells)
  size_factors[has_neighbors] <- ratio_sum_per_cell[has_neighbors] / neighbors_per_cell[has_neighbors]
  size_factors
}

# Size factor history as a matrix (cells x iterations).
.size_factor_history_matrix <- function(size_factor_history, cell_names) {
  history <- do.call(cbind, size_factor_history)
  dimnames(history) <- list(cell_names, paste0("iteration_", seq_len(ncol(history))))
  history
}

# Save the size factors and the connectivities in the Seurat object.
.save_size_factors <- function(object, history, connectivities, protein_assay) {
  object$size_factor <- apply(history, 1L, prod) # cumulative factor over all iterations
  graph <- methods::as(connectivities, "CsparseMatrix")
  dimnames(graph) <- list(colnames(object), colnames(object))
  graph <- SeuratObject::as.Graph(graph)
  SeuratObject::DefaultAssay(graph) <- protein_assay
  object[["connectivities_KNN_normalization"]] <- graph
  object
}
