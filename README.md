# KNN.normalization

<!-- badges: start -->
[![R-CMD-check](https://github.com/javier-marchena-hurtado/KNN_normalization_R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/javier-marchena-hurtado/KNN_normalization_R/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<img src="images/KNN_normalization_logo.png" width="250">

## Background and motivation

KNN normalization is a normalization method for protein counts in CITE-seq data. KNN normalization learns from neighbor cells in a KNN graph in order to estimate the appropriate total protein counts in each cell, accurately estimating total protein counts while preserving biological information.

This is the **R package** (`KNN.normalization`), which works on [Seurat](https://satijalab.org/seurat/) objects.
If you work in Python with [scanpy](https://scanpy.readthedocs.io), use the **Python package** [KNN_normalization](https://github.com/javier-marchena-hurtado/KNN_normalization) instead ([documentation](https://knn-normalization.readthedocs.io)).
For the same KNN graph, both packages give the same result; they differ only in how the KNN graph is calculated (Seurat in R, scanpy in Python).

## Installation

Install the development version of `KNN.normalization` from [GitHub](https://github.com/javier-marchena-hurtado/KNN_normalization_R):

```r
# install.packages("remotes")
remotes::install_github("javier-marchena-hurtado/KNN_normalization_R")
```

To also build the tutorial vignette, use `remotes::install_github("javier-marchena-hurtado/KNN_normalization_R", build_vignettes = TRUE)`.

## Basic usage

```r
library(Seurat)
library(KNN.normalization)

# A Seurat object with raw counts in the "counts" layer of an "RNA" and an "ADT" (protein) assay
seurat_object <- readRDS("path/to/your/data.rds")

# Run KNN normalization. The normalized protein data is saved in the "data" layer of the ADT assay.
seurat_object <- knn_normalize(seurat_object)
```

Unlike in Python, objects are not modified in place: assign the result.
The raw counts in the `counts` layer are left untouched.
See `?knn_normalize` and the tutorial (`vignette("tutorial", package = "KNN.normalization")`) for details.

## Differences with the Python package

| Python (`KNN_normalization`) | R (`KNN.normalization`) |
|---|---|
| `knn_normalize(data, ...)` on an `AnnData` or `MuData` | `knn_normalize(object, ...)` on a `Seurat` object |
| Result written to `protein.X` | Result written to the `data` layer of the protein assay (`protein_assay = "ADT"`) |
| `inplace=True/False` | Not needed: R copies on modify, so assign the result |
| `protein.obsp["connectivities"]` | A graph in `object@graphs`, selected with `graph_name` |
| `size_factor` in `.obs`, history in `.obsm` | `size_factor` in the metadata, history in `Misc(object, "knn_normalization")` |
| KNN graph with `scanpy.pp.neighbors` | KNN graph with `Seurat::FindNeighbors` |

## Release notes

See [NEWS.md](NEWS.md).

## Contact

Contributions, questions and issues are very welcome! For questions, help requests and contributions, please open an [issue](https://github.com/javier-marchena-hurtado/KNN_normalization_R/issues) on GitHub or send an [email](mailto:javier.marchena.hurtado@gmail.com).

## Citation

> t.b.a
