# =============================================================================
# print.GPAim.R
# S3 print method for GPAim objects.
# =============================================================================

#' @describeIn GPAim Print a concise summary of the genomic prediction to the
#'   console, reporting the number of genotyped lines and markers, the
#'   computational path used (\code{"vm"} or \code{"mbf"}), estimated genetic
#'   and residual variance components, narrow-sense heritability, and the
#'   range of GEBVs.
#' @param x A \code{GPAim} object.
#' @export
print.GPAim <- function(x, ...) {
    gp <- x$GP
    cat("\nGenomic Prediction (G-BLUP)\n")
    cat("===========================\n")
    cat(sprintf("  Lines genotyped : %d\n",  nrow(gp$gebv)))
    cat(sprintf("  Markers used    : %d\n",  gp$n.markers))
    cat(sprintf("  Genotype data   : %s\n",  gp$gen.type))
    cat(sprintf("  Fitting path    : %s\n",  gp$path))
    cat("\nVariance components:\n")
    cat(sprintf("  Genetic  (Vg)   : %.4f\n", gp$var.genetic))
    cat(sprintf("  Residual (Ve)   : %.4f\n", gp$var.resid))
    cat(sprintf("  Heritability h2 : %.4f\n", gp$heritability))
    cat("\nGEBV range:\n")
    cat(sprintf("  Min : %.4f\n", min(gp$gebv$GEBV, na.rm = TRUE)))
    cat(sprintf("  Mean: %.4f\n", mean(gp$gebv$GEBV, na.rm = TRUE)))
    cat(sprintf("  Max : %.4f\n", max(gp$gebv$GEBV, na.rm = TRUE)))
    invisible(x)
}
