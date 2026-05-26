# =============================================================================
# print.qtlAim.R
# S3 print method for qtlAim objects.
# =============================================================================

#' Print a Fitted \code{qtlAim} Object
#'
#' @description
#' Prints a brief narrative of the significant QTL detected by
#' \code{\link{qtlAim}} to the console, listing the chromosome, interval
#' index, and flanking marker names for each QTL. Prints a message if no
#' significant QTL were found.
#'
#' @param x A fitted object of class \code{"qtlAim"}.
#' @param genObj The \code{"wgCross"} object used in the original
#'   \code{\link{qtlAim}} call, produced by \code{\link{primeCross}}.
#'   Required to resolve internal interval indices to marker names.
#' @param \dots Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{summary.qtlAim}},
#'   \code{\link{getQTL}}
#'
#' @examples
#' \dontrun{
#' print(qtl.fit, genObj = genoRxK)
#' }
#'
#' @export
print.qtlAim <- function(x, genObj, ...) {
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("genObj must be of class \"wgCross\" produced by primeCross().")
    if (is.null(x$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
    } else {
        is.mv  <- !is.null(x$QTL$Trait)
        qtlm   <- getQTL(x, genObj)
        for (z in 1:nrow(qtlm)) {
            int  <- paste(qtlm[z, 1], qtlm[z, 2], sep = ".")
            type <- if (is.mv)
                if (x$QTL$is.interaction[z]) " [INTERACTION]" else " [MAIN]"
            else ""
            if (x$QTL$type == "interval")
                cat("\nPutative QTL found on the interval", int, type,
                    "\nLeft-hand marker is", qtlm[z, 3],
                    "\nRight-hand marker is", qtlm[z, 7], "\n")
            else
                cat("\nPutative QTL found close to marker", int, type,
                    "\nMarker is", qtlm[z, 3], "\n")
        }
    }
    invisible(x)
}
