# =============================================================================
# print.qtlAim.R
# S3 print method for qtlAim objects.
# =============================================================================

#' @describeIn qtlAim Print a brief summary of detected QTL to the console,
#'   listing the chromosome, interval number, and flanking marker names for
#'   each significant QTL. If no QTL were detected, a message is printed.
#' @param x A \code{qtlAim} object.
#' @param intervalObj The \code{"interval"} object passed to \code{qtlAim}.
#'   Required to resolve internal interval indices to marker names and positions.
#' @export
print.qtlAim <- function(x, intervalObj, ...) {
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj must be of class \"interval\".")
    if (is.null(x$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
    } else {
        qtlm <- getQTL(x, intervalObj)
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
