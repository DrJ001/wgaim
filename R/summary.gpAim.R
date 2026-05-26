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
#' @return A \code{data.frame} with one row per genotyped line (or one row per
#'   line-by-trial combination for multivariate models), sorted descending by
#'   GEBV, containing columns:
#' \describe{
#'   \item{(line identifier)}{The genetic line identifier (column name matches
#'     the \code{merge.by} argument used in \code{gpAim}).}
#'   \item{(Trait)}{Factor of environment/trial levels.  Present only for
#'     multivariate models.}
#'   \item{GEBV}{Genomic estimated breeding value, rounded to 4 decimal places.}
#'   \item{SE}{Standard error of the GEBV, rounded to 4 decimal places.}
#'   \item{Accuracy}{Prediction accuracy \eqn{\sqrt{1 - PEV/V_g}}, rounded to
#'     4 decimal places.}
#'   \item{gen.H2}{Cullis (2006) generalised \eqn{H^2}, rounded to 4 decimal
#'     places. \code{NA} for the mbf path.}
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
    if ("Accuracy" %in% names(gebv))
        gebv$Accuracy <- round(gebv$Accuracy, 4)
    if ("gen.H2" %in% names(gebv))
        gebv$gen.H2   <- round(gebv$gen.H2,   4)

    if (!is.mv) {
        gebv <- gebv[order(gebv$GEBV, decreasing = TRUE), ]
        rownames(gebv) <- NULL
        h2.str <- sprintf("h2 = %.4f", gp$heritability)
        if (!is.null(gp$gen.H2) && !is.na(gp$gen.H2))
            h2.str <- paste0(h2.str, sprintf("  gen.H2 = %.4f", gp$gen.H2))
        cat(sprintf(
            "\nGenomic Prediction summary  %s  (%d lines, %d markers)\n",
            h2.str, nrow(gebv), gp$n.markers))
        cat("GEBVs ranked highest to lowest:\n\n")
    } else {
        n.lines <- length(unique(gebv[[gp$genetic.term]]))
        gebv    <- gebv[order(gebv[[gp$Trait]], -gebv$GEBV), ]
        rownames(gebv) <- NULL
        cat(sprintf(
            "\nMultivariate Genomic Prediction  (%d lines, %d markers, %d trials)\n",
            n.lines, gp$n.markers, length(gp$trait.levels)))
        if (!is.null(gp$gen.H2) && !all(is.na(gp$gen.H2))) {
            cat("Generalised H2 per trial:\n")
            for (tname in gp$trait.levels)
                cat(sprintf("  %-12s : %.4f\n", tname, gp$gen.H2[tname]))
        }
        cat("GEBVs ranked highest to lowest within each trial:\n\n")
    }
    gebv
}
