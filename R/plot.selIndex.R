# =============================================================================
# plot.selIndex.R
# S3 plot method for selIndex objects.
#
# type = "index" (default)
#   Caterpillar-style ranked index plot. Lines are sorted by index value and
#   coloured by selection status (above/below threshold). Expected genetic gain
#   per environment is annotated.
#
# type = "heatmap"
#   For MV only. Lines x environments GEBV heatmap. Lines (rows) sorted by
#   index rank (rank 1 = top), environments (columns) in trait-level order.
#   Cell colour = GEBV value (steelblue = high, firebrick = low, white = mean).
#   A dashed horizontal rule separates selected from unselected lines, showing
#   at a glance which candidates are broadly vs specifically adapted.
#
# type = "weights"
#   Horizontal lollipop chart of the index weights b. Colour convention:
#   positive weight = steelblue, negative = firebrick. Useful for communicating
#   Smith-Hazel or desired-gains weights back to the breeder.
# =============================================================================

#' Plot a \code{selIndex} Object
#'
#' @description
#' Produces one of three diagnostic plots for a \code{\link{selIndex}} result.
#'
#' \describe{
#'   \item{\code{"index"} (default)}{Caterpillar-style plot of lines sorted by
#'     index value, coloured by selection status. Expected genetic gain per
#'     environment is annotated inside the plot.}
#'   \item{\code{"heatmap"}}{For multivariate objects only. A lines \eqn{\times}
#'     environments tile heatmap of raw GEBVs. Rows (lines) are sorted by index
#'     rank (rank 1 at top); columns are environments in trait-level order. Cell
#'     colour encodes GEBV magnitude: \code{pt.col[1]} (steelblue) = high,
#'     \code{pt.col[2]} (firebrick) = low, white = overall median. A dashed
#'     horizontal rule separates selected from unselected lines. Reveals broad
#'     vs specific adaptation at a glance: uniformly blue top rows indicate
#'     broadly adapted candidates; blue in some columns but not others signals
#'     environment-specific performers.}
#'   \item{\code{"weights"}}{Horizontal lollipop chart of the index weights
#'     \eqn{\mathbf{b}}, coloured steelblue (positive) or firebrick
#'     (negative). Communicates the derived combining weights for Smith-Hazel
#'     and desired-gains indices back to the breeder.}
#' }
#'
#' @param x A \code{"selIndex"} object returned by \code{\link{selIndex}}.
#' @param type Character string. One of \code{"index"} (default),
#'   \code{"heatmap"}, or \code{"weights"}.
#' @param prop.select Numeric proportion of lines to flag as selected.
#'   Default \code{0.10}.
#' @param threshold Numeric threshold applied to the index. If between 0 and
#'   1 exclusive, treated as a quantile; otherwise a raw index value. If
#'   \code{NULL} (default) \code{prop.select} is used.
#' @param pt.col Character vector of length 2: selected / high-GEBV colour
#'   (\code{pt.col[1]}) and unselected / low-GEBV colour (\code{pt.col[2]}).
#'   Default \code{c("steelblue", "firebrick")}.
#' @param pt.cex Numeric point size (used by \code{type = "index"} only).
#'   Default \code{0.9}.
#' @param \dots Currently unused.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{selIndex}}, \code{\link{print.selIndex}},
#'   \code{\link{summary.selIndex}}
#'
#' @examples
#' \dontrun{
#' si <- selIndex(gp.mv, weights = c(Trial1 = 2, Trial2 = 1))
#' plot(si)                              # ranked index (caterpillar)
#' plot(si, type = "heatmap")            # lines x environments GEBV heatmap
#' plot(si, type = "weights")            # index weight bar chart
#' }
#'
#' @export
plot.selIndex <- function(x,
                           type        = c("index", "heatmap", "weights"),
                           prop.select = 0.10,
                           threshold   = NULL,
                           pt.col      = c("steelblue", "firebrick"),
                           pt.cex      = 0.9,
                           ...) {
    type <- match.arg(type)

    if (type == "heatmap" && is.null(x$trait.levels))
        stop("type = \"heatmap\" requires a multivariate selIndex (gpAim with Trait != NULL).")
    if (type == "weights" && is.null(x$trait.levels))
        stop("type = \"weights\" requires a multivariate selIndex (gpAim with Trait != NULL).")

    switch(type,
        index   = .plot_si_index(x, prop.select, threshold, pt.col, pt.cex),
        heatmap = .plot_si_heatmap(x, prop.select, threshold, pt.col),
        weights = .plot_si_weights(x, pt.col, pt.cex)
    )
}


# =============================================================================
# Internal helper: compute per-environment expected genetic gain from the
# $index data frame at an arbitrary threshold.  Works for both MV and UV.
# MV:  ΔG_j = mean(GEBV_j[selected]) - mean(GEBV_j[all])  (one per env)
# UV:  ΔG   = mean(GEBV[selected])   - mean(GEBV[all])    (scalar)
# =============================================================================
.compute_gains <- function(si, thr) {
    idx <- si$index
    sel <- idx$Index >= thr
    if (!is.null(si$trait.levels)) {
        setNames(
            vapply(si$trait.levels, function(j)
                mean(idx[[j]][sel], na.rm = TRUE) - mean(idx[[j]], na.rm = TRUE),
                numeric(1L)),
            si$trait.levels)
    } else {
        mean(idx$GEBV[sel], na.rm = TRUE) - mean(idx$GEBV, na.rm = TRUE)
    }
}

# =============================================================================
# Internal: ranked index caterpillar
# =============================================================================
.plot_si_index <- function(x, prop.select, threshold, pt.col, pt.cex) {

    idx    <- x$index
    n      <- nrow(idx)
    vals   <- idx$Index

    # Determine selection: externally supplied lines override threshold logic
    ext.sel <- !is.null(x$selected)
    if (ext.sel) {
        n.sel <- x$n.selected
        gain  <- x$gain
        thr   <- NA_real_
    } else {
        thr   <- .si_threshold(vals, threshold, prop.select)
        n.sel <- sum(vals >= thr)
        gain  <- .compute_gains(x, thr)
    }

    # Rank already set in $index (descending), reverse for ascending plot
    plt.df <- idx
    plt.df$rank <- seq(n, 1)  # plot rank: 1 = worst (bottom), n = best (top)
    plt.df$col  <- if (ext.sel)
        ifelse(plt.df[[x$genetic.term]] %in% x$selected, pt.col[1], pt.col[2])
    else
        ifelse(plt.df$Index >= thr, pt.col[1], pt.col[2])

    # Annotation: index type, n selected, per-environment gains
    # Gains formatted 2 environments per line for readability with many trials
    type.labels <- c(
        "weighted"      = "Weighted",
        "smith-hazel"   = "Smith-Hazel",
        "desired-gains" = "Pesek-Baker"
    )
    gain.str <- if (!is.null(x$trait.levels)) {
        tl    <- x$trait.levels
        pairs <- split(tl, ceiling(seq_along(tl) / 2L))
        paste(vapply(pairs, function(g)
            paste(sprintf("%s:%+.3f", g, gain[g]), collapse = "  "),
            character(1L)), collapse = "\n")
    } else {
        sprintf("%+.4f", gain)
    }

    ann.label <- sprintf(
        "%s index  |  Selected: %d (%.0f%%)\nExpected gain:\n%s",
        type.labels[x$type], n.sel, 100 * n.sel / n, gain.str)

    y.range <- range(vals, na.rm = TRUE)

    gp.obj <- ggplot2::ggplot(plt.df,
                    ggplot2::aes(x = .data$rank, y = .data$Index,
                                 colour = .data$col)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
        ggplot2::geom_point(size = pt.cex, show.legend = FALSE) +
        ggplot2::scale_colour_identity() +
        ggplot2::annotate("text",
                          x     = n * 0.02,
                          y     = y.range[2],
                          label = ann.label,
                          hjust = 0, vjust = 1,
                          size  = 2.8, colour = "grey30") +
        ggplot2::labs(
            x       = "Line rank (by selection index)",
            y       = "Selection index",
            title   = sprintf("Selection Index  --  %d lines", n),
            caption = if (ext.sel)
                sprintf("Selection: %d lines supplied externally", n.sel)
            else
                sprintf("Threshold: %.4f  (top %.0f%%)", thr, 100 * prop.select)) +
        theme_scatter()

    # Threshold line only when selection came from the index
    if (!ext.sel)
        gp.obj <- gp.obj +
            ggplot2::geom_hline(yintercept = thr, colour = pt.col[1],
                                linewidth = 0.6, linetype = "dashed")
    gp.obj
}


# =============================================================================
# Internal: lines x environments GEBV heatmap
#
# Rows = lines, sorted by index rank (rank 1 = best, at the top).
# Cols = environments, in trait-level order.
# Fill  = raw GEBV value; colour scale diverges from the overall median.
# A dashed hline separates selected (top) from unselected (bottom) lines.
# =============================================================================
.plot_si_heatmap <- function(x, prop.select, threshold, pt.col) {

    tl    <- x$trait.levels
    gterm <- x$genetic.term
    idx   <- x$index
    n     <- x$n.lines
    nt    <- length(tl)

    # -------------------------------------------------------------------------
    # Determine selected lines and row ordering.
    # External selection: put selected lines (sorted by index) at the top,
    # then unselected lines below, so the separator is always a clean rule.
    # Index-based selection: standard top-N ordering.
    # -------------------------------------------------------------------------
    ext.sel <- !is.null(x$selected)
    if (ext.sel) {
        n.sel      <- x$n.selected
        sel.flag   <- idx[[gterm]] %in% x$selected
        # Sort: selected by index (best first), then unselected by index
        line.order <- c(idx[[gterm]][sel.flag], idx[[gterm]][!sel.flag])
    } else {
        thr        <- .si_threshold(idx$Index, threshold, prop.select)
        n.sel      <- sum(idx$Index >= thr)
        line.order <- idx[[gterm]]      # already rank-1-first
    }

    long.df <- do.call(rbind, lapply(tl, function(j) {
        data.frame(
            line  = idx[[gterm]],
            trial = j,
            gebv  = idx[[j]],
            stringsAsFactors = FALSE
        )
    }))
    long.df$line  <- factor(long.df$line,  levels = rev(line.order))
    long.df$trial <- factor(long.df$trial, levels = tl)

    # -------------------------------------------------------------------------
    # Colour midpoint: median GEBV across the whole matrix
    # -------------------------------------------------------------------------
    mid.val <- median(long.df$gebv, na.rm = TRUE)

    # -------------------------------------------------------------------------
    # Selection separator position on the y-axis.
    # With reversed factor levels: rank 1 = level n, rank n = level 1.
    # The n.sel selected lines occupy levels (n - n.sel + 1) to n (top of plot).
    # The separator sits between level (n - n.sel) and (n - n.sel + 1).
    # -------------------------------------------------------------------------
    sep.y <- n - n.sel + 0.5

    # -------------------------------------------------------------------------
    # Y-axis text: show line IDs when n is manageable; hide when very large
    # -------------------------------------------------------------------------
    y.text.size <- if (n <= 40) max(5, 9 - n %/% 8)
                   else if (n <= 80) 4
                   else 0   # element_blank() handled below
    y.axis.text <- if (n > 80)
        ggplot2::element_blank()
    else
        ggplot2::element_text(size = y.text.size)

    # -------------------------------------------------------------------------
    # Build plot
    # -------------------------------------------------------------------------
    gp.obj <- ggplot2::ggplot(long.df,
                               ggplot2::aes(x     = .data$trial,
                                            y     = .data$line,
                                            fill  = .data$gebv)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
        ggplot2::scale_fill_gradient2(
            low      = pt.col[2],
            mid      = "white",
            high     = pt.col[1],
            midpoint = mid.val,
            name     = "GEBV") +
        # Selection separator
        ggplot2::geom_hline(yintercept = sep.y,
                            colour     = pt.col[1],
                            linewidth  = 0.9,
                            linetype   = "dashed") +
        # "selected" label just above separator on left margin
        ggplot2::annotate("text",
                          x      = 0.5,
                          y      = sep.y + 0.25,
                          label  = sprintf("%d selected", n.sel),
                          hjust  = 0, vjust = 0,
                          size   = 2.6, colour = pt.col[1],
                          fontface = "bold") +
        # Trial names at top of plot
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::labs(
            x        = NULL,
            y        = "Line  (rank 1 at top)",
            title    = "GEBV Heatmap  --  lines sorted by selection index",
            subtitle = if (ext.sel)
                sprintf("%s index  |  %d environments  |  %d lines (external selection)",
                        switch(x$type, "weighted" = "Weighted",
                               "smith-hazel" = "Smith-Hazel",
                               "desired-gains" = "Pesek-Baker"),
                        nt, n.sel)
            else
                sprintf("%s index  |  %d environments  |  top %.0f%% selected (%d lines)",
                        switch(x$type, "weighted" = "Weighted",
                               "smith-hazel" = "Smith-Hazel",
                               "desired-gains" = "Pesek-Baker"),
                        nt, 100 * prop.select, n.sel)) +
        theme_scatter() +
        ggplot2::theme(
            axis.text.x      = ggplot2::element_text(
                                   size = 8, angle = 30, hjust = 0),
            axis.text.y      = y.axis.text,
            axis.ticks.y     = if (n > 80) ggplot2::element_blank()
                               else ggplot2::element_line(),
            legend.position  = "right",
            panel.grid       = ggplot2::element_blank()
        )

    gp.obj
}


# =============================================================================
# Internal: index weights lollipop chart
# =============================================================================
.plot_si_weights <- function(x, pt.col, pt.cex) {

    w     <- x$weights
    wdf   <- data.frame(
        env     = factor(names(w), levels = rev(names(w))),  # top-to-bottom
        weight  = as.numeric(w),
        col     = ifelse(as.numeric(w) >= 0, pt.col[1], pt.col[2]),
        stringsAsFactors = FALSE
    )

    type.labels <- c(
        "weighted"      = "Weighted (direct economic)",
        "smith-hazel"   = "Smith-Hazel optimal  [b = P\u207b\u00b9 Ga w]",
        "desired-gains" = "Pesek-Baker desired gains  [b \u221d Ga\u207b\u00b9 d]"
    )

    ggplot2::ggplot(wdf,
                    ggplot2::aes(x = .data$weight, y = .data$env,
                                 colour = .data$col)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
        ggplot2::geom_segment(
            ggplot2::aes(x = 0, xend = .data$weight,
                         y = .data$env, yend = .data$env),
            linewidth = 0.8, show.legend = FALSE) +
        ggplot2::geom_point(size = 3, show.legend = FALSE) +
        ggplot2::scale_colour_identity() +
        ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%.3f", .data$weight)),
            hjust  = ifelse(wdf$weight >= 0, -0.25, 1.25),
            size   = 3, colour = "grey25", show.legend = FALSE) +
        ggplot2::labs(
            x       = "Index weight (b)",
            y       = NULL,
            title   = "Selection Index Weights",
            subtitle = type.labels[x$type]) +
        ggplot2::expand_limits(
            x = range(wdf$weight) * 1.3) +
        theme_scatter()
}
