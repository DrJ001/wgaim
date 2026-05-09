# =============================================================================
# print.QTLAim.R
# S3 print method for QTLAim objects.
# =============================================================================

print.QTLAim <- function(x, intervalObj, ...) {
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
