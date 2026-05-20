# =============================================================================
# print.selIndex.R
# S3 print and summary methods for selIndex objects.
# =============================================================================

#' Print a \code{selIndex} Object
#'
#' @description
#' Prints a concise summary of a \code{\link{selIndex}} result, reporting the
#' index type, number of lines and environments, the index weights used,
#' the genetic correlation matrix (for multivariate objects), and the top and
#' bottom five lines by index value.
#'
#' @param x A \code{"selIndex"} object returned by \code{\link{selIndex}}.
#' @param n Integer. Number of lines to show at the top and bottom of the
#'   ranking. Default \code{5}.
#' @param \dots Currently unused.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{selIndex}}, \code{\link{summary.selIndex}},
#'   \code{\link{plot.selIndex}}
#'
#' @examples
#' \dontrun{
#' si <- selIndex(gp.mv, weights = c(Trial1 = 2, Trial2 = 1))
#' print(si)
#' }
#'
#' @export
print.selIndex <- function(x, n = 5L, ...) {

    is.mv <- !is.null(x$trait.levels)
    type.labels <- c(
        "weighted"      = "Weighted (direct economic weights)",
        "smith-hazel"   = "Smith-Hazel (optimal, b = P^{-1} Ga w)",
        "desired-gains" = "Pesek-Baker (desired gains, b proportional to Ga^{-1} d)"
    )

    cat("\nSelection Index\n")
    cat("===============\n")
    cat(sprintf("  Type           : %s\n",
                type.labels[x$type]))
    cat(sprintf("  Lines          : %d\n", x$n.lines))
    if (is.mv) {
        cat(sprintf("  Environments   : %d  (%s)\n",
                    x$n.environments, paste(x$trait.levels, collapse = ", ")))
        cat(sprintf("  Standardised   : %s\n",
                    if (x$standardise) "yes (unit variance per environment)" else "no"))
    }

    # Index weights
    cat("\nIndex weights (b):\n")
    if (is.mv) {
        wdf <- data.frame(
            Environment = names(x$weights),
            Weight      = round(x$weights, 4),
            row.names   = NULL
        )
        print(wdf, row.names = FALSE)
    } else {
        cat(sprintf("  %s (univariate -- index equals GEBV)\n", x$genetic.term))
    }

    # Genetic correlations (MV only)
    if (is.mv && !is.null(x$Gcor)) {
        cat("\nGenetic correlations (Gcor):\n")
        print(round(x$Gcor, 3))
    }

    # Top/bottom n lines
    n <- min(n, floor(x$n.lines / 2))
    idx <- x$index

    cat(sprintf("\nTop %d lines by index:\n", n))
    print(head(idx, n), row.names = FALSE)

    cat(sprintf("\nBottom %d lines by index:\n", n))
    print(tail(idx, n), row.names = FALSE)

    invisible(x)
}


#' Summary of a \code{selIndex} Object
#'
#' @description
#' Returns and prints the full ranked index table from a
#' \code{\link{selIndex}} result, with one row per line sorted from highest
#' to lowest index value.  The header reports the index type, combining weights,
#' number selected, and the expected genetic gain per environment at the
#' selection threshold fixed when \code{\link{selIndex}} was called.
#'
#' @param object A \code{"selIndex"} object returned by \code{\link{selIndex}}.
#' @param \dots Currently unused.
#'
#' @return A \code{data.frame} with columns for the line identifier, one GEBV
#'   column per environment (multivariate only), \code{Index}, and \code{Rank},
#'   sorted descending by \code{Index}.
#'
#' @seealso \code{\link{selIndex}}, \code{\link{print.selIndex}},
#'   \code{\link{plot.selIndex}}
#'
#' @examples
#' \dontrun{
#' summary(si)
#' top10 <- head(summary(si), 10)
#' }
#'
#' @export
summary.selIndex <- function(object, ...) {

    idx     <- object$index
    is.mv   <- !is.null(object$trait.levels)

    type.labels <- c(
        "weighted"      = "Weighted",
        "smith-hazel"   = "Smith-Hazel",
        "desired-gains" = "Pesek-Baker desired-gains"
    )

    cat(sprintf(
        "\nSelection Index (%s)  --  %d lines",
        type.labels[object$type], object$n.lines))
    if (is.mv)
        cat(sprintf(", %d environments", object$n.environments))
    cat("\n")

    w.str <- paste(
        sprintf("%s: %.3f", names(object$weights), object$weights),
        collapse = "  |  ")
    cat(sprintf("Weights : %s\n", w.str))
    if (!is.null(object$selected)) {
        cat(sprintf("Selected: %d lines  (externally supplied)\n",
                    object$n.selected))
    } else {
        cat(sprintf("Selected: %d lines (top %.0f%%)  |  threshold = %.4f\n",
                    object$n.selected,
                    100 * object$prop.select,
                    object$thr))
    }

    # Per-environment expected genetic gain
    cat("\nExpected genetic gain (original GEBV scale):\n")
    if (is.mv) {
        gain.df <- data.frame(
            Environment  = names(object$gain),
            "Delta.G"    = round(object$gain, 4),
            check.names  = FALSE,
            row.names    = NULL
        )
        names(gain.df)[2] <- "Delta G"
        print(gain.df, row.names = FALSE)
    } else {
        cat(sprintf("  %+.4f\n", object$gain))
    }
    cat("\n")

    idx
}
