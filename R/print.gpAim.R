# =============================================================================
# print.gpAim.R
# S3 print method for gpAim objects.
# =============================================================================

#' Print a Fitted \code{gpAim} Object
#'
#' @description
#' Prints a concise summary of a \code{\link{gpAim}} fit to the console,
#' reporting the number of genotyped lines and markers, the computational path
#' used (\code{"vm"} or \code{"mbf"}), estimated genetic and residual variance
#' components, narrow-sense heritability, and the range of GEBVs.
#'
#' @param x A fitted object of class \code{"gpAim"}.
#' @param \dots Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{gpAim}}, \code{\link{summary.gpAim}},
#'   \code{\link{plot.gpAim}}
#'
#' @examples
#' \dontrun{
#' print(gp.fit)
#' }
#'
#' @export
print.gpAim <- function(x, ...) {
    gp    <- x$GP
    is.mv <- !is.null(gp$Trait)

    if (!is.mv) {
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
        if (!is.null(gp$gen.H2) && !is.na(gp$gen.H2))
            cat(sprintf("  Generalised H2  : %.4f\n", gp$gen.H2))
        cat("\nGEBV range:\n")
        cat(sprintf("  Min : %.4f\n", min(gp$gebv$GEBV, na.rm = TRUE)))
        cat(sprintf("  Mean: %.4f\n", mean(gp$gebv$GEBV, na.rm = TRUE)))
        cat(sprintf("  Max : %.4f\n", max(gp$gebv$GEBV, na.rm = TRUE)))
    } else {
        cat("\nGenomic Prediction (Multivariate G-BLUP)\n")
        cat("=========================================\n")
        cat(sprintf("  %-16s: %s\n", gp$Trait,
                    paste(gp$trait.levels, collapse = ", ")))
        cat(sprintf("  Lines genotyped : %d\n",
                    length(unique(gp$gebv[[gp$genetic.term]]))))
        cat(sprintf("  Markers used    : %d\n", gp$n.markers))
        cat(sprintf("  Genotype data   : %s\n", gp$gen.type))
        cat(sprintf("  Fitting path    : %s\n", gp$path))

        cat("\nVariance components (per trial):\n")
        vc.df <- data.frame(
            Vg     = round(gp$var.genetic,   4),
            Ve     = round(gp$var.resid,     4),
            h2     = round(gp$heritability,  4),
            check.names = FALSE
        )
        if (!is.null(gp$gen.H2) && !all(is.na(gp$gen.H2)))
            vc.df$gen.H2 <- round(gp$gen.H2[gp$trait.levels], 4)
        rownames(vc.df) <- gp$trait.levels
        print(vc.df)

        cat("\nGenetic correlation (Gcor):\n")
        print(round(gp$Gcor, 4))

        cat("\nGEBV range (per trial):\n")
        for (tname in gp$trait.levels) {
            tv <- gp$gebv$GEBV[gp$gebv[[gp$Trait]] == tname]
            cat(sprintf("  %-12s  Min = %8.4f   Mean = %8.4f   Max = %8.4f\n",
                        tname,
                        min(tv,  na.rm = TRUE),
                        mean(tv, na.rm = TRUE),
                        max(tv,  na.rm = TRUE)))
        }
    }
    invisible(x)
}
