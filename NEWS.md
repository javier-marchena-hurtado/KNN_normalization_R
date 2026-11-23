# KNN.normalization 0.1.0

* First release: an R/Seurat port of the Python package
  [KNN_normalization](https://github.com/javier-marchena-hurtado/KNN_normalization) (v0.1.2).
* `knn_normalize()` normalizes the protein assay of a Seurat object, computing the KNN graph
  from the RNA data, from the protein data, or reusing an existing graph.
* `normalize_with_neighbors()` normalizes with a precomputed KNN graph. For the same graph,
  it gives the same result as the Python package.
