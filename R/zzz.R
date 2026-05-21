# Suppress R CMD check NOTEs for asreml functions and other symbols that are
# resolved at runtime but have no visible binding at check time.
# asreml is a commercial package listed under SystemRequirements and cannot
# be declared in Imports or Suggests, so its functions must be suppressed here.
utils::globalVariables(c(
    # asreml functions used as bare calls throughout the engine
    "asreml", "asreml.options", "predict.asreml", "update.asreml",
    "asreml.is.converged", "mbf", "vm",
    # ggplot2 / ggrepel aesthetics accessed via .data$ or bare names
    "bottom", "label", "lh.pos", "name", "pos", "rh.pos", "top", "trait"
))

.onAttach <- function(libname, pkgname){
  if(!("asreml" %in% loadedNamespaces())) {
    packageStartupMessage("The analysis functions in wgAim require the R package ASReml-R to be installed.
This is currently a commercially available product with a licensing system that varies
depending on the institution.

Please visit https://vsni.co.uk/software/asreml-r/ for more information including trial
licenses and pricing.")
      }
}

