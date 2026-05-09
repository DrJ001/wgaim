# =============================================================================
# print.GWASAim.R
# S3 print method for GWASAim objects.
# =============================================================================

print.GWASAim <- function(x, panelObj, ...) {
    if (missing(panelObj))
        stop("panelObj is a required argument.")
    if (!inherits(panelObj, "panel"))
        stop("panelObj must be of class \"panel\".")
    bonf.str <- if (isTRUE(x$QTL$bonferroni))
        sprintf("Bonferroni-adjusted threshold: %s (%d markers)",
                formatC(x$QTL$TypeI.eff, format = "e", digits = 3),
                x$QTL$n.markers)
    else
        sprintf("Significance threshold: %.4f (no Bonferroni correction)", x$QTL$TypeI.eff)
    cat("\nGWAS forward selection -", bonf.str, "\n")
    if (is.null(x$QTL$effects)) {
        cat("No significant markers detected.\n")
    } else {
        qtlm <- getQTL(x, panelObj)
        for (z in 1:nrow(qtlm))
            cat(sprintf("\nSignificant marker on chromosome %s: %s (%.2f cM)\n",
                        qtlm[z, 1], qtlm[z, 3], as.numeric(qtlm[z, 4])))
    }
}
