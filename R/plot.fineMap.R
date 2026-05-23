#' Plot fine-mapping results
#'
#' @description
#' Produces a ggplot2 LOD profile for each fine-mapped QTL / marker.  Each
#' QTL is shown in its own facet panel with two labelled reference points:
#' \describe{
#'   \item{Peak}{The grid position with the maximum LOD score, shown in
#'     \code{peak.col}.}
#'   \item{QTL}{The original QTL / marker position (centre of the fine-mapping
#'     window) as returned by \code{qtlAim} / \code{gwasAim}, shown in
#'     \code{qtl.col} with a vertical dashed reference line.  Suppressed when
#'     the peak coincides with the original position.}
#' }
#' Labels are placed with \code{ggrepel} to avoid overlap.
#'
#' @param x An object of class \code{"fineMap"} produced by
#'   \code{\link{fineMap}}.
#' @param col Line colour.  Default \code{"steelblue"}.
#' @param peak.col Colour for the fine-map peak point and label.  Default
#'   \code{"firebrick"}.
#' @param qtl.col Colour for the original QTL reference point, label, and
#'   vertical dashed line.  Default \code{"navy"}.
#' @param ref.col Colour of the vertical dashed reference line.  Defaults to
#'   \code{qtl.col}.
#' @param ncol Number of columns in the \code{facet_wrap} layout.  Defaults
#'   to \code{NULL}, which uses \code{min(n_qtl, 4)} columns (so up to four
#'   QTL appear on one row; five or more wrap automatically).
#' @param \dots Currently ignored.
#'
#' @return A \code{ggplot} object (invisibly) with attribute
#'   \code{suggested_height} giving an appropriate plot height in inches.
#'   The plot is also printed.
#'
#' @seealso \code{\link{fineMap}}, \code{\link{print.fineMap}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline
#'   facet_wrap labs scale_y_continuous coord_cartesian
#'   theme_bw theme element_text element_blank
#' @importFrom ggrepel geom_text_repel
#' @exportS3Method
plot.fineMap <- function(x, col = "steelblue",
                         peak.col = "firebrick", qtl.col = "navy",
                         ref.col = qtl.col, ncol = NULL, ...)
{

    # -------------------------------------------------------------------------
    # 1. Build combined long data frame
    # -------------------------------------------------------------------------
    dfs <- lapply(names(x), function(nm) {
        df <- x[[nm]]
        if (nrow(df) == 0L) return(NULL)
        df$qtl <- nm
        df
    })
    dfs <- do.call("rbind", Filter(Negate(is.null), dfs))

    if (is.null(dfs) || nrow(dfs) == 0L) {
        message("No fine-mapping data to plot.")
        return(invisible(NULL))
    }
    # Append "(window clipped)" to the facet label for any QTL whose scan
    # was truncated by a chromosome boundary -- using the flag stored by fineMap().
    qtl_labels <- stats::setNames(names(x), names(x))
    for (nm in names(x)) {
        if (isTRUE(attr(x[[nm]], "clipped")))
            qtl_labels[nm] <- paste0(nm, "  (window clipped)")
    }
    dfs$qtl <- factor(qtl_labels[as.character(dfs$qtl)],
                      levels = qtl_labels[names(x)])

    # -------------------------------------------------------------------------
    # 2. y-axis: LOD score; cap Inf values (can arise when pvalue = 0 exactly)
    # -------------------------------------------------------------------------
    if (any(is.infinite(dfs$LOD))) {
        max_fin <- max(dfs$LOD[is.finite(dfs$LOD)], na.rm = TRUE)
        dfs$LOD[is.infinite(dfs$LOD)] <- max_fin * 1.1
    }
    dfs$y_plot <- dfs$LOD
    y_lab      <- "LOD score"

    # -------------------------------------------------------------------------
    # 3. Peak point: maximum LOD per QTL (after transformation so y_plot is set)
    # -------------------------------------------------------------------------
    peaks <- do.call("rbind", lapply(split(dfs, dfs$qtl), function(d) {
        d[which.max(d$LOD), , drop = FALSE]
    }))
    peaks$label_text <- peaks$mark

    # -------------------------------------------------------------------------
    # 4. Original QTL reference point at the exact qtl_pos x position.
    #    y is linearly interpolated from the two bracketing grid points so
    #    the navy diamond sits precisely on the profile curve.
    #    The vline is drawn at the same x, guaranteeing alignment.
    # -------------------------------------------------------------------------
    ref_list <- lapply(names(x), function(nm) {
        df  <- x[[nm]]
        qp  <- attr(df, "qtl_pos")
        if (is.null(qp) || nrow(df) == 0L) return(NULL)
        sub <- dfs[dfs$qtl == nm, , drop = FALSE]

        if (nrow(sub) == 0L) return(NULL)

        # For qtlAim the fine grid is built from qtl_pos outward, so qtl_pos
        # is an exact grid point.  For gwasAim the significant marker is
        # always in the scan.  Either way, find the nearest grid point to qp
        # (within floating-point tolerance) and use it directly.
        near_idx <- which.min(abs(sub$dist - qp))
        row            <- sub[near_idx, , drop = FALSE]
        row$label_text <- sprintf("%.1f cM", qp)
        row
    })
    ref_df <- do.call("rbind", Filter(Negate(is.null), ref_list))

    # Suppress the QTL reference label when it is within 1 cM of the peak
    # (avoids duplicated annotations when the fine-map peak = original QTL).
    if (!is.null(ref_df) && nrow(ref_df) > 0L) {
        suppressions <- vapply(seq_len(nrow(ref_df)), function(i) {
            nm <- as.character(ref_df$qtl[i])
            pk <- peaks[peaks$qtl == nm, , drop = FALSE]
            abs(ref_df$dist[i] - pk$dist[1L]) < 1.0
        }, logical(1L))
        ref_df <- ref_df[!suppressions, , drop = FALSE]
    }

    # -------------------------------------------------------------------------
    # 5. Vertical reference line -- use ref_df$dist (the snapped grid position)
    #    so the line passes exactly through the navy diamond.
    #    Only drawn when the original QTL is not suppressed (i.e. not at peak).
    # -------------------------------------------------------------------------
    vline_df <- if (!is.null(ref_df) && nrow(ref_df) > 0L) {
        data.frame(
            qtl        = ref_df$qtl,
            xintercept = ref_df$dist,
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }

    # -------------------------------------------------------------------------
    # 6. Layout
    # -------------------------------------------------------------------------
    n_facets <- length(unique(dfs$qtl))
    n_cols   <- if (!is.null(ncol)) as.integer(ncol) else min(n_facets, 4L)
    n_rows   <- ceiling(n_facets / n_cols)

    # -------------------------------------------------------------------------
    # 7. Build plot
    # -------------------------------------------------------------------------
    gp <- ggplot(dfs, aes(x = .data[["dist"]], y = .data[["y_plot"]])) +
        geom_line(colour = col, linewidth = 0.7)

    # Vertical reference line at original QTL position
    if (!is.null(vline_df) && nrow(vline_df) > 0L) {
        gp <- gp + geom_vline(
            data        = vline_df,
            aes(xintercept = .data[["xintercept"]]),
            linetype    = "dashed",
            colour      = ref.col,
            linewidth   = 0.4
        )
    }

    # Combine peak and QTL reference rows into one data frame so that a single
    # geom_text_repel() call sees all labels simultaneously -- the only way
    # ggrepel can avoid clashes between the two sets.
    peaks$pt_type <- "Peak"
    peaks$pt_shape <- 16L

    if (!is.null(ref_df) && nrow(ref_df) > 0L) {
        ref_df$pt_type  <- "QTL"
        ref_df$pt_shape <- 18L
        ann_df <- rbind(peaks, ref_df)
    } else {
        ann_df <- peaks
    }

    pt_colours <- c(Peak = peak.col, QTL = qtl.col)
    pt_shapes  <- c(Peak = 16L,      QTL = 18L)

    gp <- gp +
        geom_point(
            data  = ann_df,
            aes(x     = .data[["dist"]],
                y     = .data[["y_plot"]],
                colour = .data[["pt_type"]],
                shape  = .data[["pt_type"]]),
            size  = 2.5,
            show.legend = FALSE
        ) +
        ggrepel::geom_text_repel(
            data  = ann_df,
            aes(x     = .data[["dist"]],
                y     = .data[["y_plot"]],
                label  = .data[["label_text"]],
                colour = .data[["pt_type"]]),
            size  = 3,
            segment.size = 0.3,
            min.segment.length = 0,
            show.legend = FALSE
        ) +
        scale_colour_manual(values = pt_colours) +
        scale_shape_manual( values = pt_shapes)

    gp <- gp +
        facet_wrap(~ qtl, scales = "free", ncol = n_cols) +
        labs(x = "Position (cM)", y = y_lab) +
        scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.25))) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            strip.text       = element_text(size = 9, face = "bold"),
            axis.title       = element_text(size = 9),
            panel.grid.minor = element_blank()
        )

    attr(gp, "suggested_height") <- n_rows * 3.5
    print(gp)
    invisible(gp)
}
