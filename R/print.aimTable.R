# =============================================================================
# print.aimTable.R
# S3 print method for aimTable objects.
# =============================================================================

#' Print an \code{aimTable} Object
#'
#' @description
#' Prints a stacked multi-trait QTL or GWAS summary table to the console.
#' A brief header reports the analysis type (\code{qtlAim} or \code{gwasAim})
#' and the number of traits included, followed by the table formatted without
#' row names. Prints a \code{"(no QTL detected)"} notice for empty tables.
#'
#' @param x An \code{aimTable} object produced by \code{\link{aimTable}}.
#' @param \dots Further arguments passed to \code{\link{print.data.frame}}.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{aimTable}}, \code{\link{summary.qtlAim}},
#'   \code{\link{summary.gwasAim}}
#'
#' @examples
#' \dontrun{
#' tbl <- aimTable(qtl.yld, qtl.tgw, genObj = genoRxK)
#' print(tbl)
#' }
#'
#' @export
print.aimTable <- function(x, ...) {
    obj_class  <- attr(x, "obj.class")
    n_traits   <- if (nrow(x) == 0L) 0L else length(unique(x$Trait))
    trait_word <- if (n_traits == 1L) "trait" else "traits"
    analysis   <- if (is.null(obj_class)) "qtlAim/gwasAim" else obj_class

    cat(sprintf("\nAim table  [%s  |  %d %s]\n\n",
                analysis, n_traits, trait_word))

    if (nrow(x) == 0L) {
        cat("  (no QTL detected in any model)\n\n")
    } else {
        print.data.frame(x, row.names = FALSE, ...)
        cat("\n")
    }
    invisible(x)
}
