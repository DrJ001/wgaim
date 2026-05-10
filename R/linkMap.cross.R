# =============================================================================
# linkMap.cross.R
# S3 linkMap() method for "cross" objects.
#
# Chromosome bar geometry
# -----------------------
# Two overlapping thick vertical segments with lineend = "round" create a
# pill-shaped ideogram (BAR_BORDER_LW outer dark segment, BAR_FILL_LW inner
# white segment).
#
# Tick geometry
# -------------
# Ticks are drawn as the FIRST layer, before the bar.  The bar segments are
# then drawn on top, painting over the tick interior.  Ticks therefore emerge
# from the bar's right visual edge without any gap, regardless of the panel
# width.  TICK_X0 = 0 (bar centre); TICK_X1 controls the total tick span in
# data coordinates.
#
# Panel x-axis pinning
# --------------------
# geom_blank points at X_LEFT_LIMIT / X_RIGHT_LIMIT force every panel to the
# same fixed x-axis range.  Without them, the range is determined by ggrepel
# nudge_x (which ggplot2 includes in scale training), coupling tick proportions
# to label layout.  Matching xlim = c(...) in the ggrepel calls prevents labels
# from escaping those boundaries.
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helper: geometry constants used by linkMap.cross and linkMap.qtlAim
# -----------------------------------------------------------------------------
.lm_const <- function() {
    list(
        ## Bar linewidths (mm, device space)
        BAR_BORDER_LW  = 5.25,
        BAR_FILL_LW    = 3.75,
        BAR_BORDER_COL = "grey20",
        BAR_FILL_COL   = "white",

        ## Tick x-range in data coordinates.
        ## TICK_X0 = 0 so ticks start at bar centre; bar drawn on top hides
        ## the interior, ticks emerge from the right visual edge.
        ## TICK_X1 = 10% of panel width (X_RIGHT_LIMIT - X_LEFT_LIMIT = 0.36)
        TICK_X0 = 0,
        TICK_X1 = 0.036,

        ## Right-side (marker) label anchor / nudge
        LABEL_X       =  0.036,
        LABEL_NUDGE_X =  0.02,

        ## Left-side (trait / QTL) label anchor / nudge
        TRAIT_X       =  0.00,
        TRAIT_NUDGE_X = -0.06,

        ## Panel x-axis boundaries (pinned by geom_blank)
        X_RIGHT_LIMIT =  0.22,
        X_LEFT_LIMIT  = -0.14
    )
}

#' @title Plot a Genetic Linkage Map (Generic)
#' @description Generic function for plotting genetic linkage maps. Dispatches
#'   to \code{\link{linkMap.cross}} for \code{"cross"} objects and
#'   \code{\link{linkMap.qtlAim}} for \code{qtlAim} objects.
#' @param object An object for which a \code{linkMap} method exists.
#' @param ... Arguments passed to the appropriate method.
#' @export
linkMap <- function(object, ...) UseMethod("linkMap")

#' Plot a Genetic Linkage Map
#'
#' @description
#' Produces a clean ggplot2 visualisation of a genetic linkage map using
#' frameless facets — one column per chromosome — with a shared reversed cM
#' y-axis (0 cM at the top). Chromosome bars are rendered as narrow pill
#' shapes. Marker labels are placed with \code{ggrepel} to avoid overprinting
#' even on dense maps. Returns a \code{ggplot} object so additional layers
#' (e.g. QTL regions from \code{\link{linkMap.qtlAim}}) can be composed
#' with \code{+}.
#'
#' @param object A \code{"cross"} object produced by
#'   \code{\link[qtl]{read.cross}} or \code{\link[wgaim]{cross2int}}.
#' @param chr Optional character vector of chromosome names to display.
#'   Defaults to all chromosomes.
#' @param chr.dist Optional named list with elements \code{$start} and/or
#'   \code{$end} (scalar or per-chromosome numeric vectors of cM positions)
#'   restricting the displayed map range.
#' @param marker.names Character string controlling marker label display.
#'   \code{"markers"} (default): marker names; \code{"dist"}: cM positions;
#'   \code{"none"} or \code{NULL}: no labels (best for dense maps);
#'   \code{"flanking"}: only markers listed in
#'   \code{attr(object, "flanking")} (set automatically by
#'   \code{\link{linkMap.qtlAim}}).
#' @param tick Logical. Unused; retained for back-compatibility.
#' @param squash Logical. Unused; retained for back-compatibility.
#' @param m.cex Numeric. Text size (pt) for marker labels. Default \code{7}.
#' @param nrow Integer or \code{"auto"}.  Number of rows in the facet layout.
#'   \code{1} (default) gives the standard single-row display.  \code{"auto"}
#'   chooses \code{ceiling(n_chr / 8)} automatically.  Ignored when
#'   \code{row.chr} is supplied.
#' @param row.chr Optional named or unnamed list of character vectors, one
#'   element per row.  Each vector names the chromosomes to display in that
#'   row, in left-to-right order.  When supplied, it overrides both \code{chr}
#'   and \code{nrow}: only chromosomes present in \code{row.chr} are shown, and
#'   the facet layout is set to \code{length(row.chr)} rows by
#'   \code{max(lengths(row.chr))} columns.  Chromosomes absent from the map are
#'   silently dropped; rows that contain fewer chromosomes than the widest row
#'   are right-padded with blank panels.
#' @param ... Additional arguments passed to
#'   \code{ggrepel::\link[ggrepel]{geom_text_repel}}.
#' @return A \code{ggplot} object.
#' @importFrom ggrepel geom_text_repel
#' @rdname linkMap.cross
#' @export
linkMap.cross <- function(object, chr, chr.dist,
                          marker.names = "markers",
                          tick = FALSE, squash = TRUE,
                          m.cex = 7,
                          nrow = 1, row.chr = NULL,
                          ...) {

    K <- .lm_const()

    ## -------------------------------------------------------------------------
    ## 1.  Extract map, apply chr / row.chr selection, determine facet layout
    ## -------------------------------------------------------------------------
    map <- qtl::pull.map(object)

    if (!is.null(row.chr)) {
        ## ---- Option 2: row.chr drives both selection and layout -------------
        if (!is.list(row.chr) ||
            !all(vapply(row.chr, is.character, logical(1L))))
            stop("'row.chr' must be a list of character vectors.")

        ## Drop any chromosomes not present in the map
        row.chr <- lapply(row.chr, function(r) r[r %in% names(map)])
        row.chr <- Filter(function(r) length(r) > 0L, row.chr)
        if (length(row.chr) == 0L)
            stop("None of the chromosomes in 'row.chr' are present in the map.")

        ## Reorder map to the row.chr sequence
        ordered_chrs <- unlist(row.chr, use.names = FALSE)
        map <- map[ordered_chrs]

        facet_nrow <- length(row.chr)
        facet_ncol <- max(lengths(row.chr))

    } else {
        ## ---- Option 1 / default: chr arg + nrow ---------------------------
        if (!missing(chr)) {
            if (any(is.na(pmatch(chr, names(map)))))
                stop("Some names in 'chr' do not match chromosome names in the map.")
            map <- map[chr]
        }

        n.chr_pre <- length(map)
        facet_nrow <- if (identical(nrow, "auto"))
            ceiling(n.chr_pre / 8L)
        else
            as.integer(nrow)
        facet_ncol <- NULL
    }

    ## ---- chr.dist range filter (applies after chr / row.chr selection) ----
    if (!missing(chr.dist)) {
        if (!all(names(chr.dist) %in% c("start", "end")))
            stop("'chr.dist' must be a named list with elements 'start' and/or 'end'.")
        n        <- length(map)
        start_cm <- if (!is.null(chr.dist$start)) rep_len(chr.dist$start, n) else rep(0, n)
        end_cm   <- if (!is.null(chr.dist$end))   rep_len(chr.dist$end,   n) else
                        vapply(map, max, numeric(1))
        map <- mapply(function(m, s, e) m[m >= s & m <= e],
                      map, start_cm, end_cm, SIMPLIFY = FALSE)
        map <- Filter(function(m) length(m) > 0L, map)
    }

    n.chr <- length(map)
    if (n.chr == 0L) stop("No chromosomes remain after subsetting.")

    ## -------------------------------------------------------------------------
    ## 2.  Decide which markers to label
    ## -------------------------------------------------------------------------
    flanking_set <- attr(object, "flanking")
    show_labels  <- !is.null(marker.names) && !identical(marker.names, "none")

    label_map <- if (show_labels) {
        if (identical(marker.names, "flanking") && !is.null(flanking_set))
            lapply(map, function(m) m[names(m) %in% flanking_set])
        else
            map
    } else NULL

    ## -------------------------------------------------------------------------
    ## 3.  Build the data frames
    ## -------------------------------------------------------------------------

    ## 3a. Chromosome bar endpoints (one row per chromosome)
    bar_df <- data.frame(
        chr    = factor(names(map), levels = names(map)),
        bottom = vapply(map, max, numeric(1)),
        top    = vapply(map, min, numeric(1)),
        stringsAsFactors = FALSE
    )

    ## 3b. Marker tick positions (one row per marker)
    marker_df <- do.call(rbind, lapply(names(map), function(ch) {
        m <- map[[ch]]
        data.frame(chr  = factor(ch, levels = names(map)),
                   pos  = as.numeric(m),
                   name = names(m),
                   stringsAsFactors = FALSE)
    }))

    ## 3c. Label positions (subset of marker_df)
    label_df <- if (show_labels && !is.null(label_map)) {
        do.call(rbind, lapply(names(label_map), function(ch) {
            m <- label_map[[ch]]
            if (length(m) == 0L) return(NULL)
            data.frame(
                chr  = factor(ch, levels = names(map)),
                pos  = as.numeric(m),
                name = if (identical(marker.names, "dist"))
                           as.character(round(as.numeric(m), 1))
                       else names(m),
                stringsAsFactors = FALSE
            )
        }))
    } else NULL

    ## -------------------------------------------------------------------------
    ## 4.  Build the base plot
    ##
    ## LAYER ORDER:
    ##   0. geom_blank anchors  — pin x-axis range (invisible)
    ##   1. Marker ticks        — drawn before bar; bar hides interior
    ##   2. Bar border          — pill shape, dark outer ring
    ##   3. Bar fill            — white interior
    ## -------------------------------------------------------------------------
    p <- ggplot2::ggplot() +

        ## ------------------------------------------------------------------
        ## Layer 0: Invisible anchors that fix the panel x-axis range.
        ## ggrepel's nudge_x is included in ggplot2's scale training, so
        ## without these the panel width is slave to the nudge value.
        ggplot2::geom_blank(
            data    = bar_df,
            mapping = ggplot2::aes(x = K$X_RIGHT_LIMIT, y = top)
        ) +
        ggplot2::geom_blank(
            data    = bar_df,
            mapping = ggplot2::aes(x = K$X_LEFT_LIMIT, y = top)
        ) +

        ## ------------------------------------------------------------------
        ## Layer 1: Marker ticks (drawn BEFORE bar so bar hides interior)
        ggplot2::geom_segment(
            data      = marker_df,
            ggplot2::aes(x = K$TICK_X0, xend = K$TICK_X1, y = pos, yend = pos),
            linewidth = 0.25,
            colour    = "grey40"
        ) +

        ## ------------------------------------------------------------------
        ## Layer 2: Bar border — pill shape, dark outer ring
        ggplot2::geom_segment(
            data      = bar_df,
            ggplot2::aes(x = 0, xend = 0, y = top, yend = bottom),
            linewidth = K$BAR_BORDER_LW,
            lineend   = "round",
            colour    = K$BAR_BORDER_COL
        ) +

        ## ------------------------------------------------------------------
        ## Layer 3: Bar fill — white, slightly narrower than border
        ggplot2::geom_segment(
            data      = bar_df,
            ggplot2::aes(x = 0, xend = 0, y = top, yend = bottom),
            linewidth = K$BAR_FILL_LW,
            lineend   = "round",
            colour    = K$BAR_FILL_COL
        ) +

        ## Facet — one column per chromosome, x suppressed
        ggplot2::facet_wrap(~ chr, nrow = facet_nrow, ncol = facet_ncol,
                            scales = "free_x") +

        ## Reversed y-axis: 0 cM at top
        ggplot2::scale_y_reverse(name = "Location (cM)") +

        ## Frameless minimalistic theme
        .linkMap_theme()

    ## -------------------------------------------------------------------------
    ## 5.  Optionally add marker labels via ggrepel
    ## -------------------------------------------------------------------------
    if (show_labels && !is.null(label_df) && nrow(label_df) > 0L) {
        p <- p + ggrepel::geom_text_repel(
            data               = label_df,
            mapping            = ggplot2::aes(x = K$LABEL_X, y = pos,
                                              label = name),
            hjust              = 0,
            size               = m.cex / ggplot2::.pt,
            direction          = "y",
            nudge_x            = K$LABEL_NUDGE_X,
            xlim               = c(K$LABEL_X, K$X_RIGHT_LIMIT),
            segment.size       = 0.2,
            segment.colour     = "grey55",
            min.segment.length = 0,
            max.overlaps       = Inf,
            force              = 0.5,
            box.padding        = 0.08,
            ...
        )
    }

    p
}

# -----------------------------------------------------------------------------
# .lm_base_args — build the base_args list for do.call(linkMap, ...)
# -----------------------------------------------------------------------------
## Shared by qtlAim, gwasAim (single and multi paths).
## chr.dist = NULL means "not supplied".
## Returns a list ready for do.call(linkMap, ...).
.lm_base_args <- function(cross_obj, chr = NULL, row.chr = NULL, nrow = 1,
                            flanking_marks = NULL, marker.names = "markers",
                            m.cex = 7, chr.dist = NULL, pass_args = list()) {
    show_marks <- !is.null(marker.names) && !identical(marker.names, "none")
    mn <- if (!show_marks) "none" else
          if (!is.null(flanking_marks)) "flanking" else marker.names
    if (!is.null(flanking_marks) && show_marks)
        attr(cross_obj, "flanking") <- flanking_marks

    if (!is.null(row.chr)) {
        args <- list(object = cross_obj, row.chr = row.chr,
                     m.cex = m.cex, marker.names = mn, nrow = nrow)
    } else {
        args <- list(object = cross_obj, m.cex = m.cex,
                     marker.names = mn, nrow = nrow)
        if (!is.null(chr)) args$chr <- chr
    }
    if (!is.null(chr.dist)) args$chr.dist <- chr.dist
    c(args, pass_args)
}

# -----------------------------------------------------------------------------
# .lm_left_labels — left-side ggrepel (trait / model name labels)
# -----------------------------------------------------------------------------
## df must have columns: chr (factor), pos (numeric), label (character).
## colour may be scalar or a per-row vector.
.lm_left_labels <- function(p, df, K, m.cex, colour) {
    p + ggrepel::geom_text_repel(
        data               = df,
        mapping            = ggplot2::aes(x = K$TRAIT_X, y = pos,
                                          label = label),
        hjust              = 1,
        size               = m.cex / ggplot2::.pt,
        colour             = colour,
        direction          = "y",
        nudge_x            = K$TRAIT_NUDGE_X,
        xlim               = c(K$X_LEFT_LIMIT, K$TRAIT_X),
        segment.size       = 0.2,
        segment.colour     = "grey55",
        min.segment.length = 0,
        max.overlaps       = Inf,
        force              = 0.6,
        box.padding        = 0.08
    )
}

# -----------------------------------------------------------------------------
# .lm_right_labels — right-side ggrepel (marker name / distance labels)
# -----------------------------------------------------------------------------
.lm_right_labels <- function(p, df, K, m.cex, colour) {
    p + ggrepel::geom_text_repel(
        data               = df,
        mapping            = ggplot2::aes(x = K$LABEL_X, y = pos,
                                          label = label),
        hjust              = 0,
        size               = m.cex / ggplot2::.pt,
        colour             = colour,
        direction          = "y",
        nudge_x            = K$LABEL_NUDGE_X,
        xlim               = c(K$LABEL_X, K$X_RIGHT_LIMIT),
        segment.size       = 0.2,
        segment.colour     = "grey55",
        min.segment.length = 0,
        max.overlaps       = Inf,
        force              = 0.5,
        box.padding        = 0.08
    )
}

# -----------------------------------------------------------------------------
# linkMap.default — backward-compatible dispatcher for list-based calls
# -----------------------------------------------------------------------------
## Handles the old API:  linkMap(list(qtl1, qtl2), intervalObj = x)
## Unwraps the list and delegates to the appropriate typed method, which
## then detects the extra models via the ... inspection path.
#' @export
linkMap.default <- function(object, ...) {
    if (!is.list(object))
        stop("No linkMap method for class '",
             paste(class(object), collapse = ", "), "'.")
    types <- unique(vapply(object, function(x) class(x)[1L], character(1L)))
    if (length(types) != 1L)
        stop("All elements of 'object' must have the same class.")
    ## Unwrap: first element becomes 'object', remainder go into ... where the
    ## method's ... inspection picks them up as extra models.
    do.call(
        switch(types,
               qtlAim  = linkMap.qtlAim,
               gwasAim = linkMap.gwasAim,
               stop("No multi-model linkMap method for class '", types, "'.")),
        c(list(object[[1L]]), object[-1L], list(...))
    )
}

# -----------------------------------------------------------------------------
# Shared theme for all linkMap variants
# -----------------------------------------------------------------------------
.linkMap_theme <- function() {
    ggplot2::theme(
        panel.background  = ggplot2::element_blank(),
        panel.border      = ggplot2::element_blank(),
        panel.grid        = ggplot2::element_blank(),
        panel.spacing.x   = ggplot2::unit(0.3, "lines"),
        strip.background  = ggplot2::element_blank(),
        strip.text        = ggplot2::element_text(
                                size   = 8,
                                face   = "bold",
                                margin = ggplot2::margin(b = 3)),
        axis.text.x       = ggplot2::element_blank(),
        axis.ticks.x      = ggplot2::element_blank(),
        axis.title.x      = ggplot2::element_blank(),
        axis.line.y       = ggplot2::element_line(linewidth = 0.4),
        axis.text.y       = ggplot2::element_text(size = 8),
        axis.title.y      = ggplot2::element_text(size = 9),
        legend.position   = "right",
        legend.title      = ggplot2::element_text(size = 8),
        legend.text       = ggplot2::element_text(size = 7),
        plot.title        = ggplot2::element_text(
                                size   = 10,
                                face   = "bold",
                                hjust  = 0.5,
                                margin = ggplot2::margin(b = 6))
    )
}
