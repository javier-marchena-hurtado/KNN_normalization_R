# KNN.normalization

<!-- badges: start -->
[![R-CMD-check](https://github.com/javier-marchena-hurtado/KNN_normalization_R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/javier-marchena-hurtado/KNN_normalization_R/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<img src="man/figures/KNN_normalization_logo.png" alt="KNN normalization logo" width="250">

## Background and motivation

KNN normalization is a normalization method for protein counts in CITE-seq data. KNN normalization learns from neighbor cells in a KNN graph in order to estimate the appropriate total protein counts in each cell, accurately estimating total protein counts while preserving biological information.

This is the **R package** (`KNN.normalization`), which works on Seurat objects. If you work in Python with scanpy, use the **Python package** [KNN_normalization](https://github.com/javier-marchena-hurtado/KNN_normalization) instead. For the same KNN graph, both packages give the same result; they differ only in how the KNN graph is calculated (Seurat in R, scanpy in Python).

## Installation

Install the development version of `KNN.normalization` from [GitHub](https://github.com/javier-marchena-hurtado/KNN_normalization_R):

```r
# install.packages("remotes")
remotes::install_github("javier-marchena-hurtado/KNN_normalization_R")
```

To also install the tutorial vignette locally, use `remotes::install_github("javier-marchena-hurtado/KNN_normalization_R", build_vignettes = TRUE)`.

## Basic usage

```r
library(Seurat)
library(KNN.normalization)

# A Seurat object with raw counts in the "counts" layer of an "RNA" and an "ADT" (protein) assay
seurat_object <- readRDS("path/to/your/data.rds")

# Run KNN normalization. The normalized protein data is saved in the "data" layer of the ADT assay.
seurat_object <- knn_normalize(seurat_object)
```

See `?knn_normalize` for details. For a full walkthrough (neighbors from RNA or protein, downstream analysis, your own KNN graph), see the [tutorial][].

## Documentation

Please refer to the [documentation website][], in particular the [tutorial][] and the [documentation of the knn_normalize() function][].
If you installed the package with `build_vignettes = TRUE`, the tutorial is also available offline with `vignette("tutorial", package = "KNN.normalization")`.

## Release notes

See [NEWS.md](NEWS.md).

## Contact

Contributions, questions and issues are very welcome! For questions, help requests and contributions, please open an [issue](https://github.com/javier-marchena-hurtado/KNN_normalization_R/issues) on GitHub or send an [email](mailto:javier.marchena.hurtado@gmail.com).

## Citation

> t.b.a

[documentation website]: https://javier-marchena-hurtado.github.io/KNN_normalization_R/
[tutorial]: https://javier-marchena-hurtado.github.io/KNN_normalization_R/articles/tutorial.html
[documentation of the knn_normalize() function]: https://javier-marchena-hurtado.github.io/KNN_normalization_R/reference/knn_normalize.html
