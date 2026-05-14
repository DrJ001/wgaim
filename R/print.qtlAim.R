# =============================================================================
# print.qtlAim.R
# S3 print method for qtlAim objects.
# =============================================================================

#' @describeIn qtlAim Print a brief summary of detected QTL to the console,
#'   listing the chromosome, interval number, and flanking marker names for
#'   each significant QTL. If no QTL were detected, a message is printed.
#' @param x A \code{qtlAim} object.
#' @param genObj The \code{"wgCross"} object passed to \code{qtlAim},
#'   produced by \code{\link{primeCross}}. Required to resolve internal
#'   interval indices to marker names and positions.
#' @export
print.qtlAim <- function(x, genObj, ...) {
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("genObj must be of class \"wgCross\" produced by primeCross().")
    if (is.null(x$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
    } else {
        qtlm <- getQTL(x, genObj)
        for (z in 1:nrow(qtlm)) {
            int <- paste(qtlm[z, 1], qtlm[z, 2], sep = ".")
            if (x$QTL$type == "interval")
                cat("\nPutative QTL found on the interval", int,
                    "\nLeft-hand marker is", qtlm[z, 3],
                    "\nRight-hand marker is", qtlm[z, 7], "\n")
            else
                cat("\nPutative QTL found close to marker", int,
                    "\nMarker is", qtlm[z, 3], "\n")
        }
    }
}
