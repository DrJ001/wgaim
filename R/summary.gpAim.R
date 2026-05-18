# =============================================================================
# summary.gpAim.R
# S3 summary method for gpAim objects.
# Returns the GEBV data frame sorted descending by GEBV value.
# =============================================================================

#' Summary of a Fitted \code{gpAim} Model
#'
#' @description
#' Returns and prints the full genomic estimated breeding value (GEBV) table
#' from a \code{\link{gpAim}} fit, sorted from highest to lowest GEBV. A
#' header line reports the estimated narrow-sense heritability, number of
#' lines, and number of markers used.
#'
#' @param object A fitted object of class \code{"gpAim"}, as returned by
#'   \code{\link{gpAim}}.
#' @param \dots Currently unused.
#'
#' @return A \code{data.frame} with one row per genotyped line, sorted
#'   descending by GEBV, containing columns:
#' \describe{
#'   \item{(line identifier)}{The genetic line identifier (column name
#'     matches the \code{merge.by} argument used in \code{gpAim}).}
#'   \item{GEBV}{Genomic estimated breeding value, rounded to 4 decimal
#'     places.}
#'   \item{SE}{Standard error of the GEBV, rounded to 4 decimal places.}
#' }
#'
#' @seealso \code{\link{gpAim}}, \code{\link{print.gpAim}},
#'   \code{\link{plot.gpAim}}
#'
#' @examples
#' \dontrun{
#' # After running gpAim():
#' summary(gp.fit)
#' top10 <- head(summary(gp.fit), 10)
#' }
#'
#' @export
summary.gpAim <- function(object, ...) {
    gp    <- object$GP
    gebv  <- gp$gebv
    is.mv <- !is.null(gp$Trait)

    gebv$GEBV <- round(gebv$GEBV, 4)
    gebv$SE   <- round(gebv$SE,   4)

    if (!is.mv) {
        gebv <- gebv[order(gebv$GEBV, decreasing = TRUE), ]
        rownames(gebv) <- NULL
        cat(sprintf(
            "\nGenomic Prediction summary  h2 = %.4f  (%d lines, %d markers)\n",
            gp$heritability, nrow(gebv), gp$n.markers))
        cat("GEBVs ranked highest to lowest:\n\n")
    } else {
        n.lines <- length(unique(gebv[[gp$genetic.term]]))
        gebv    <- gebv[order(gebv[[gp$Trait]], -gebv$GEBV), ]
        rownames(gebv) <- NULL
        cat(sprintf(
            "\nMultivariate Genomic Prediction  (%d lines, %d markers, %d trials)\n",
            n.lines, gp$n.markers, length(gp$trait.levels)))
        cat("GEBVs ranked highest to lowest within each trial:\n\n")
    }
    gebv
}
