# =============================================================================
# summary.GPAim.R
# S3 summary method for GPAim objects.
# Returns the GEBV data frame sorted descending by GEBV value.
# =============================================================================

summary.GPAim <- function(object, ...) {
    gp   <- object$GP
    gebv <- gp$gebv
    gebv <- gebv[order(gebv$GEBV, decreasing = TRUE), ]
    rownames(gebv) <- NULL
    gebv$GEBV <- round(gebv$GEBV, 4)
    gebv$SE   <- round(gebv$SE,   4)
    cat(sprintf("\nGenomic Prediction summary  h2 = %.4f  (%d lines, %d markers)\n",
                gp$heritability, nrow(gebv), gp$n.markers))
    cat("GEBVs ranked highest to lowest:\n\n")
    gebv
}
