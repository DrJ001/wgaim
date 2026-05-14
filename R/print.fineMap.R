#' Print a fineMap object
#'
#' @description
#' Prints a compact tabular summary of fine-mapping results for each QTL or
#' marker in a \code{"fineMap"} object.  For each mapped region the table
#' shows the peak position (minimum p-value), its cM location, p-value and
#' LOD score, followed by the full scan data frame.
#'
#' @param x An object of class \code{"fineMap"} produced by
#'   \code{\link{fineMap}}.
#' @param \dots Currently ignored.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{fineMap}}, \code{\link{plot.fineMap}}
#'
#' @exportS3Method
print.fineMap <- function(x, ...)
{
    type <- attr(x, "type")
    cat("\nFine-mapping results (", type, ")\n", sep = "")
    cat(strrep("=", 50), "\n", sep = "")

    for (nm in names(x)) {
        df <- x[[nm]]
        cat("\nQTL:", nm, "\n")
        if (nrow(df) == 0L) {
            cat("  No positions scanned.\n")
            next
        }
        peak_idx  <- which.max(df$LOD)
        cat(sprintf("  Peak: %s  dist = %.2f cM  p = %.4f  LOD = %.4f\n",
                    df$mark[peak_idx], df$dist[peak_idx],
                    df$pvalue[peak_idx], df$LOD[peak_idx]))
        cat("\n")
        print(df, row.names = FALSE, digits = 4)
    }
    invisible(x)
}
