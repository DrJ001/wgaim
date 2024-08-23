# Precompile vignettes that depend on ASreml
# CRAN wants vignettes in raw form (i.e. Rmd files), but they don't have ASreml
# So we pre-knit the .Rmd.orig file, which outputs it as a plain .Rmd file with
# no code that needs to be run, and CRAN can just build that.

knitr::knit("vignettes/wgaim_intro.Rmd_orig", "vignettes/wgaim_intro.Rmd")

# Because the vignette was knit at the top level of the folder, we need to move the
# images generated from the top level directory to the vignettes folder
file.copy(list.files(pattern = ".png"), "vignettes/", recursive = T, overwrite = T)
file.remove(list.files(pattern = ".png"))
