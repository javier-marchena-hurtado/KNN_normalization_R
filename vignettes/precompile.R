# The tutorial needs data that is not part of the package, so it is precomputed:
# run this script from the package root whenever vignettes/tutorial.Rmd.orig changes.
old_wd <- setwd("vignettes")
on.exit(setwd(old_wd))
knitr::knit("tutorial.Rmd.orig", "tutorial.Rmd")
