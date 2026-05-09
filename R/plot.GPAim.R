# =============================================================================
# plot.GPAim.R
# S3 plot method for GPAim objects.
# Dot plot of GEBVs ranked highest to lowest, coloured by sign,
# with a vertical zero line and heritability annotation.
# =============================================================================

plot.GPAim <- function(x, pt.col = c("steelblue", "firebrick"),
                        pt.cex = 0.8, ...) {
    gp   <- x$GP
    gebv <- gp$gebv
    gebv <- gebv[order(gebv$GEBV), ]
    gebv$rank <- seq_len(nrow(gebv))
    gebv$col  <- ifelse(gebv$GEBV >= 0, pt.col[1], pt.col[2])

    gp.line  <- gp$genetic.term
    h2.label <- sprintf("h^2 == %.3f", gp$heritability)

    ggplot(gebv, aes_string(x = "rank", y = "GEBV", colour = "col")) +
        geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
        geom_point(size = pt.cex, show.legend = FALSE) +
        scale_colour_identity() +
        annotate("text", x = nrow(gebv) * 0.05,
                 y = max(gebv$GEBV, na.rm = TRUE) * 0.9,
                 label = h2.label, parse = TRUE,
                 size = 3.5, hjust = 0, colour = "grey30") +
        labs(x = "Line rank", y = "GEBV",
             title = "Genomic Estimated Breeding Values") +
        theme_scatter()
}
