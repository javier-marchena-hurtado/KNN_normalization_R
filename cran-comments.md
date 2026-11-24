## Test environments

* local Linux (Debian 13), R 4.5.3, `R CMD check --as-cran`
* win-builder, R Under development (2026-11-17 r90926 ucrt), Windows Server 2022
* GitHub Actions (ubuntu, macOS, windows; R release and devel)
* TODO before submitting: win-builder (release) and R-hub, then update this file with their results

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* "KNN" is flagged as a possibly misspelled word in the DESCRIPTION. It is an acronym (k-nearest neighbors) and the name of the method implemented by the package.
* (Local only: notes about a missing `qpdf`/`tidy` and an unverifiable clock come from the local machine, not from the package.)

## Notes

The tutorial vignette is precomputed (`vignettes/tutorial.Rmd` is knitted from `vignettes/tutorial.Rmd.orig`, see `vignettes/precompile.R`), because it uses data that is not part of the package. Building the vignette therefore does not need the data or any network access.
