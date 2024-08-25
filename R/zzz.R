.onAttach <- function(libname, pkgname){
  if(!("asreml" %in% loadedNamespaces())) {
    packageStartupMessage("The QTL analysis functions in this package require the R package ASReml-R to be installed.
This is currently a commercially available product with a licensing system that varies
depending on the institution.

Please visit https://vsni.co.uk/software/asreml-r/ for more information including trial
licenses and pricing.")
      }
}

