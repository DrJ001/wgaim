# =============================================================================
# plot.gpAim.R
# S3 plot method for gpAim objects.
#
# type = "caterpillar" (default)
#   Ranked GEBV plot with half-HSD comparison bars. Half-HSD bars are shown
#   only for the top.n and bottom.n lines (near selection decisions); the
#   remaining lines are plotted as simple points. Non-overlapping half-HSD
#   bars indicate a statistically significant difference between two GEBVs.
#   half-HSD bar at line i = t_alpha * SE_i / sqrt(2).
#
# type = "heatmap"
#   Genomic relationship matrix heatmap, lines ordered by hierarchical
#   clustering. Reveals population structure, families, and duplicates.
#
# type = "density"
#   Density plot of the GEBV distribution with a vertical selection threshold
#   line, shaded selected region, and annotation of n selected and
#   selection differential.
# =============================================================================

#' @describeIn gpAim Produce one of three diagnostic plots for a \code{gpAim}
#'   object, selected via the \code{type} argument.
#'
#'   \describe{
#'     \item{\code{"caterpillar"} (default)}{Ranked GEBV dot plot with
#'       half-HSD comparison bars for the top and bottom \code{top.n} lines.
#'       Half-HSD bars are defined as \eqn{t_\alpha \times SE_i / \sqrt{2}};
#'       two GEBVs whose half-HSD bars do not overlap are significantly
#'       different at level \code{alpha}. Bars are shown only for the extreme
#'       lines (near selection decisions) to preserve readability with many
#'       lines. Points are coloured by selection status relative to the
#'       threshold.}
#'     \item{\code{"heatmap"}}{Tile heatmap of the genomic relationship matrix
#'       \eqn{G}, with lines ordered by Ward hierarchical clustering. Reveals
#'       population sub-structure, family groups, and potential duplicate lines.
#'       A warning is issued for datasets with more than 500 lines where
#'       rendering may be slow.}
#'     \item{\code{"density"}}{Density plot of the GEBV distribution with the
#'       selected region shaded, a threshold line, and annotation of the number
#'       selected, proportion selected, and selection differential
#'       (mean GEBV of selected minus mean of all).}
#'   }
#'
#' @param x A \code{gpAim} object.
#' @param type Character string specifying the plot type: \code{"caterpillar"}
#'   (default), \code{"heatmap"}, or \code{"density"}.
#' @param alpha Numeric significance level used to compute half-HSD bar widths
#'   in the caterpillar plot. Default \code{0.05}.
#' @param top.n Integer. Number of highest- and lowest-ranked lines for which
#'   half-HSD bars are drawn in the caterpillar plot. Default \code{20}.
#'   Reducing this value improves readability with very many lines.
#' @param threshold Numeric. Selection threshold for caterpillar and density
#'   plots. If between 0 and 1 exclusive, treated as a quantile of the GEBV
#'   distribution (e.g. \code{0.9} selects the top 10\%). If outside (0, 1),
#'   treated as a raw GEBV value. If \code{NULL} (default), \code{prop.select}
#'   is used instead.
#' @param prop.select Numeric proportion of lines to select, used when
#'   \code{threshold = NULL}. Default \code{0.10} (top 10\%).
#' @param pt.col Character vector of length 2: \code{pt.col[1]} is used for
#'   selected lines / high relatedness / positive GEBVs;
#'   \code{pt.col[2]} for unselected lines / low relatedness / negative GEBVs.
#'   Default \code{c("steelblue","firebrick")}.
#' @param pt.cex Numeric point size. Default \code{0.8}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
plot.gpAim <- function(x,
                       type        = c("caterpillar", "heatmap", "density"),
                       # --- caterpillar args ---
                       alpha       = 0.05,   # significance level for half-HSD bars
                       top.n       = 20,     # show bars for top and bottom n lines only
                       threshold   = NULL,   # selection threshold: value or proportion (0-1)
                       # --- density args ---
                       prop.select = 0.10,   # proportion to select (used if threshold=NULL)
                       # --- shared ---
                       pt.col      = c("steelblue", "firebrick"),
                       pt.cex      = 0.8,
                       ...) {

    type <- match.arg(type)
    gp   <- x$GP

    switch(type,
        caterpillar = .plot_gp_caterpillar(gp, alpha, top.n, threshold,
                                           prop.select, pt.col, pt.cex),
        heatmap     = .plot_gp_heatmap(gp, pt.col),
        density     = .plot_gp_density(gp, threshold, prop.select, pt.col)
    )
}

# =============================================================================
# Internal: caterpillar plot
# =============================================================================
.plot_gp_caterpillar <- function(gp, alpha, top.n, threshold, prop.select,
                                  pt.col, pt.cex) {
    gebv    <- gp$gebv[order(gp$gebv$GEBV), ]
    n       <- nrow(gebv)
    gebv$rank <- seq_len(n)

    # Half-HSD bars: bar_i = t_alpha * SE_i / sqrt(2)
    # Non-overlapping bars => significantly different GEBVs
    t.val      <- qnorm(1 - alpha / 2)
    gebv$half.hsd <- t.val * gebv$SE / sqrt(2)
    gebv$lower    <- gebv$GEBV - gebv$half.hsd
    gebv$upper    <- gebv$GEBV + gebv$half.hsd

    # Determine selection threshold
    thr <- .gp_threshold(gebv$GEBV, threshold, prop.select)
    n.sel   <- sum(gebv$GEBV >= thr)
    sel.dif <- mean(gebv$GEBV[gebv$GEBV >= thr]) - mean(gebv$GEBV)

    # Colour: selected (above threshold) vs not
    gebv$col <- ifelse(gebv$GEBV >= thr, pt.col[1], pt.col[2])

    # Half-HSD bars shown only for top.n and bottom.n lines
    top.n   <- min(top.n, floor(n / 2))
    bar.idx <- c(1:top.n, (n - top.n + 1):n)
    bar.df  <- gebv[bar.idx, ]

    # Annotation label
    ann.label <- sprintf(
        "Selected: %d  (%.0f%%)\nSel. differential: %+.3f\nh² = %.3f",
        n.sel, 100 * n.sel / n, sel.dif, gp$heritability)

    gp.obj <- ggplot(gebv, aes_string(x = "rank", y = "GEBV", colour = "col")) +
        geom_hline(yintercept = 0,   colour = "grey70", linewidth = 0.3) +
        geom_hline(yintercept = thr, colour = pt.col[1],
                   linewidth = 0.5, linetype = "dashed") +
        geom_point(size = pt.cex, show.legend = FALSE) +
        # Half-HSD bars for top/bottom n only
        geom_errorbar(data = bar.df,
                      aes_string(ymin = "lower", ymax = "upper"),
                      width = 0, linewidth = 0.4, alpha = 0.7,
                      show.legend = FALSE) +
        scale_colour_identity() +
        annotate("text",
                 x = n * 0.02, y = max(gebv$upper, na.rm = TRUE) * 0.95,
                 label = ann.label, hjust = 0, vjust = 1,
                 size = 2.8, colour = "grey30") +
        labs(x     = "Line rank",
             y     = "GEBV",
             title = sprintf("Genomic Estimated Breeding Values  [%s path]",
                             gp$path),
             caption = sprintf(
                 "Bars: half-HSD (\u03b1=%.2f) for top/bottom %d lines. Non-overlapping bars indicate significant difference.",
                 alpha, top.n)) +
        theme_scatter()

    gp.obj
}

# =============================================================================
# Internal: genomic relationship heatmap
# =============================================================================
.plot_gp_heatmap <- function(gp, pt.col) {
    G <- gp$rel.matrix
    n <- nrow(G)

    if (n > 500)
        warning(sprintf(
            "Heatmap with %d lines produces a large plot. Consider subsetting.", n))

    # Hierarchical clustering to order lines by genetic similarity
    G.dist  <- as.dist(1 - G / max(G, na.rm = TRUE))
    hc      <- hclust(G.dist, method = "ward.D2")
    ord     <- hc$order
    G.ord   <- G[ord, ord]
    lnames  <- rownames(G.ord)

    # Long format for ggplot2 (column-major matches as.vector)
    g.df <- data.frame(
        Line1 = factor(rep(lnames, times = n), levels = lnames),
        Line2 = factor(rep(lnames, each  = n), levels = lnames),
        value = as.vector(G.ord),
        stringsAsFactors = FALSE
    )

    # Colour limits centred on the median (expected genomic relationship)
    med.val <- median(G.ord, na.rm = TRUE)
    max.val <- max(abs(G.ord - med.val), na.rm = TRUE)
    lims    <- c(med.val - max.val, med.val + max.val)

    ggplot(g.df, aes_string(x = "Line1", y = "Line2", fill = "value")) +
        geom_tile() +
        scale_fill_gradient2(
            low      = pt.col[2],
            mid      = "white",
            high     = pt.col[1],
            midpoint = med.val,
            limits   = lims,
            name     = "Relatedness") +
        labs(title   = "Genomic Relationship Matrix",
             subtitle = sprintf("%d lines, ordered by Ward hierarchical clustering", n),
             x = NULL, y = NULL) +
        theme_scatter() +
        theme(axis.text  = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "right")
}

# =============================================================================
# Internal: GEBV density with selection threshold
# =============================================================================
.plot_gp_density <- function(gp, threshold, prop.select, pt.col) {
    gebv.vals <- gp$gebv$GEBV
    thr       <- .gp_threshold(gebv.vals, threshold, prop.select)
    n.sel     <- sum(gebv.vals >= thr)
    prop.sel  <- n.sel / length(gebv.vals)
    sel.dif   <- mean(gebv.vals[gebv.vals >= thr]) - mean(gebv.vals)

    df <- data.frame(GEBV = gebv.vals)

    # Build density for shading the selected region
    dens      <- density(gebv.vals)
    dens.df   <- data.frame(x = dens$x, y = dens$y)
    shade.df  <- dens.df[dens.df$x >= thr, ]
    # Close the shaded polygon at y = 0
    shade.df  <- rbind(
        data.frame(x = thr, y = 0),
        shade.df,
        data.frame(x = max(shade.df$x), y = 0))

    ann.label <- sprintf(
        "Selected: %d  (%.0f%%)\nSel. differential: %+.3f\nh² = %.3f",
        n.sel, round(100 * prop.sel), sel.dif, gp$heritability)

    ggplot(df, aes_string(x = "GEBV")) +
        geom_polygon(data = shade.df, aes_string(x = "x", y = "y"),
                     fill = pt.col[1], alpha = 0.3, inherit.aes = FALSE) +
        geom_line(data = dens.df, aes_string(x = "x", y = "y"),
                  colour = "grey30", linewidth = 0.7, inherit.aes = FALSE) +
        geom_vline(xintercept = thr, colour = pt.col[1],
                   linetype = "dashed", linewidth = 0.7) +
        annotate("text",
                 x     = thr + diff(range(gebv.vals)) * 0.02,
                 y     = max(dens$y) * 0.95,
                 label = ann.label,
                 hjust = 0, vjust = 1,
                 size  = 3, colour = "grey20") +
        labs(x     = "GEBV",
             y     = "Density",
             title = sprintf("GEBV Distribution  [%s path]", gp$path),
             caption = sprintf(
                 "Selection threshold: %.4f  (%s)",
                 thr,
                 if (!is.null(threshold) && threshold < 1 && threshold > 0)
                     sprintf("top %.0f%%", 100 * (1 - threshold))
                 else sprintf("top %.0f%%", 100 * prop.select))) +
        theme_scatter()
}

# =============================================================================
# Helper: resolve selection threshold from value or proportion
# =============================================================================
.gp_threshold <- function(gebv.vals, threshold, prop.select) {
    if (!is.null(threshold)) {
        if (threshold > 0 && threshold < 1)
            quantile(gebv.vals, threshold, names = FALSE)   # treat as quantile
        else
            threshold                                        # treat as raw value
    } else {
        quantile(gebv.vals, 1 - prop.select, names = FALSE)
    }
}
