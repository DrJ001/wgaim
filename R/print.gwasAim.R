# =============================================================================
# print.gwasAim.R
# S3 print method for gwasAim objects.
# =============================================================================

#' Print a Fitted \code{gwasAim} Object
#'
#' @description
#' Prints a brief narrative of the significant markers detected by
#' \code{\link{gwasAim}} to the console, listing the chromosome, marker name,
#' and cM position of each association. Also reports the significance threshold
#' and the total number of markers in the panel.
#'
#' @param x A fitted object of class \code{"gwasAim"}.
#' @param genObj The \code{"wgPanel"} object used in the original
#'   \code{\link{gwasAim}} call, produced by \code{\link{primePanel}}.
#' @param \dots Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{gwasAim}}, \code{\link{summary.gwasAim}},
#'   \code{\link{getQTL}}
#'
#' @examples
#' \dontrun{
#' print(gwas.fit, genObj = panel)
#' }
#'
#' @export
print.gwasAim <- function(x, genObj, ...) {
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgPanel"))
        stop("genObj must be of class \"wgPanel\" produced by primePanel().")
    cat(sprintf("\nGWAS forward selection  TypeI = %.4f  (%d markers in panel)\n",
                x$QTL$TypeI, x$QTL$n.markers))
    if (is.null(x$QTL$effects)) {
        cat("No significant markers detected.\n")
    } else {
        is.mv <- !is.null(x$QTL$Trait)
        qtlm  <- getQTL(x, genObj)
        for (z in 1:nrow(qtlm)) {
            type <- if (is.mv)
                if (x$QTL$is.interaction[z]) " [INTERACTION]" else " [MAIN]"
            else ""
            cat(sprintf("\nSignificant marker on chromosome %s: %s (%.2f cM)%s\n",
                        qtlm[z, 1], qtlm[z, 3], as.numeric(qtlm[z, 4]), type))
        }
    }
    invisible(x)
}
