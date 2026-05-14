# Precompile vignettes that depend on ASReml.
# CRAN does not have ASReml-R, so vignettes must be pre-knitted to plain Rmd
# files (no executable code chunks) before submission.
#
# Usage: run this script from the package root before submitting to CRAN.
#   Rscript build_vignette.R
#
# TODO: add entries here for each new wgAim vignette as they are created.

# Example (uncomment when vignette source exists):
# knitr::knit("vignettes/wgAim_intro.Rmd_orig", "vignettes/wgAim_intro.Rmd")
# file.copy(list.files(pattern = ".png"), "vignettes/", recursive = TRUE, overwrite = TRUE)
# file.remove(list.files(pattern = ".png"))
