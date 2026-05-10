# =============================================================================
# print.gwasAim.R
# S3 print method for gwasAim objects.
# =============================================================================

#' @describeIn gwasAim Print a brief summary of significant markers to the
#'   console, listing the chromosome, marker name, and cM position of each
#'   detected association. Also reports the significance threshold and total
#'   number of markers tested.
#' @param x A \code{gwasAim} object.
#' @param panelObj The \code{"panel"} object passed to \code{gwasAim}.
#' @export
print.gwasAim <- function(x, panelObj, ...) {
    if (missing(panelObj))
        stop("panelObj is a required argument.")
    if (!inherits(panelObj, "panel"))
        stop("panelObj must be of class \"panel\".")
    cat(sprintf("\nGWAS forward selection  TypeI = %.4f  (%d markers in panel)\n",
                x$QTL$TypeI, x$QTL$n.markers))
    if (is.null(x$QTL$effects)) {
        cat("No significant markers detected.\n")
    } else {
        qtlm <- getQTL(x, panelObj)
        for (z in 1:nrow(qtlm))
            cat(sprintf("\nSignificant marker on chromosome %s: %s (%.2f cM)\n",
                        qtlm[z, 1], qtlm[z, 3], as.numeric(qtlm[z, 4])))
    }
}
