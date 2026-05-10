# =============================================================================
# linkMap.gwasAim.R
# S3 linkMap() method for gwasAim objects — single model OR multiple models.
#
# Multiple models are passed as extra positional arguments in ...:
#   linkMap(gwas1, gwas2, panelObj = my_panel)
#
# Any unnamed ... arg that inherits "gwasAim" is treated as an additional
# model.  Named ... args are passed through to linkMap.cross().
#
# The old list-based API still works via linkMap.default (backward compat):
#   linkMap(list(gwas1, gwas2), panelObj = my_panel)
#
# Shared helpers (.lm_base_args, .lm_left_labels, .lm_right_labels) live in
# linkMap.cross.R.
# =============================================================================

#' @importFrom ggrepel geom_text_repel
#' @describeIn gwasAim Plot the genetic linkage map with significant GWAS
#'   markers overlaid as coloured horizontal bands.  Accepts a single
#'   \code{gwasAim} object or multiple models as extra positional arguments —
#'   e.g. \code{linkMap(gwas1, gwas2, panelObj = my_panel)} — for a
#'   multi-trait overlay with automatic colour legend.  Returns a
#'   \code{ggplot} object.
#' @param object A \code{gwasAim} object, or the first of several.
#' @param panelObj The \code{"panel"} object used in the analysis.
#' @param chr Optional character vector of chromosomes to display.
#' @param marker.names \code{"markers"} (default), \code{"dist"},
#'   \code{"none"} or \code{NULL}.
#' @param marker.colour Colour for the significant-marker bands (single model).
#' @param label.colour Colour for right-side marker name labels.  Defaults to
#'   \code{marker.colour}.
#' @param trait.labels Optional character vector of trait names (one per model).
#' @param trait.colour Colour(s) for left-side trait labels / multi-model bands.
#'   Defaults to \code{marker.colour}.
#' @param band.lwd Band linewidth (mm).  Defaults to \code{BAR_FILL_LW}.
#' @param m.cex Text size (pt) for labels.  Default \code{7}.
#' @param nrow Integer or \code{"auto"}.  Facet row count.
#' @param row.chr Optional list of character vectors for manual row layout.
#' @param ... Additional \code{gwasAim} models (positional, unnamed) and/or
#'   arguments passed to \code{\link{linkMap.cross}}.
#' @return A \code{ggplot} object.
#' @export
linkMap.gwasAim <- function(object, ...,
                             panelObj,
                             chr,
                             marker.names  = "markers",
                             marker.colour = "firebrick",
                             label.colour  = marker.colour,
                             trait.labels  = NULL,
                             trait.colour  = marker.colour,
                             band.lwd      = NULL,
                             m.cex         = 7,
                             nrow          = 1,
                             row.chr       = NULL) {

    K <- .lm_const()
    if (is.null(band.lwd)) band.lwd <- K$BAR_FILL_LW

    ## -------------------------------------------------------------------------
    ## Separate extra gwasAim models from pass-through ... args
    ## -------------------------------------------------------------------------
    dots   <- list(...)
    dnames <- if (is.null(names(dots))) rep("", length(dots)) else
              ifelse(is.na(names(dots)), "", names(dots))
    is_extra     <- vapply(dots, inherits, logical(1L), "gwasAim") & dnames == ""
    extra_models <- dots[is_extra]
    pass_args    <- dots[!is_extra]

    ## -------------------------------------------------------------------------
    ## Guards (shared across both paths)
    ## -------------------------------------------------------------------------
    if (missing(panelObj))
        stop("'panelObj' is a required argument.")
    if (!inherits(panelObj, "interval"))
        stop("'panelObj' must be of class \"panel\" / \"interval\".")

    ## Build fake_cross once — shared by both paths
    fake_cross <- list(geno = panelObj$geno)
    class(fake_cross) <- "cross"
    panel_map <- lapply(panelObj$geno, function(ch) ch$map)

    ## Internal: resolve marker keys → (chr, name, pos) data frame
    .resolve_markers <- function(gwas_obj) {
        keys   <- gwas_obj$QTL$qtl
        if (is.null(keys) || length(keys) == 0L) return(NULL)
        parsed <- strsplit(keys, "\\.")
        sig_chr <- vapply(parsed, `[[`, character(1L), 2L)
        sig_idx <- as.integer(vapply(parsed, `[[`, character(1L), 3L))
        sig_name <- mapply(function(ch, idx) {
            m <- panel_map[[ch]]
            if (is.null(m) || idx < 1L || idx > length(m)) NA_character_
            else names(m)[idx]
        }, sig_chr, sig_idx)
        sig_pos <- mapply(function(ch, idx) {
            m <- panel_map[[ch]]
            if (is.null(m) || idx < 1L || idx > length(m)) NA_real_
            else unname(m[idx])
        }, sig_chr, sig_idx)
        valid <- !is.na(sig_pos)
        if (!any(valid)) return(NULL)
        data.frame(chr = sig_chr[valid], name = sig_name[valid],
                   pos = sig_pos[valid], stringsAsFactors = FALSE)
    }

    ## =========================================================================
    ## MULTI-MODEL PATH
    ## =========================================================================
    if (length(extra_models) > 0L) {

        all_models <- c(list(object), extra_models)
        n_models   <- length(all_models)

        ## --- Trait names -----------------------------------------------------
        if (is.null(trait.labels))
            trait_names <- vapply(all_models,
                function(m) as.character(m$call$fixed[[2L]]), character(1L))
        else
            trait_names <- rep_len(as.character(trait.labels), n_models)

        if (anyDuplicated(trait_names)) {
            cnt <- stats::ave(seq_len(n_models), trait_names,
                              FUN = function(x) seq_along(x))
            dup <- duplicated(trait_names) | duplicated(trait_names, fromLast = TRUE)
            trait_names[dup] <- paste0(trait_names[dup], cnt[dup])
        }

        ## --- Colours ---------------------------------------------------------
        if (length(trait.colour) == 1L && trait.colour == marker.colour)
            trait.colour <- grDevices::hcl.colors(n_models, palette = "Dark 2")
        band_cols  <- stats::setNames(rep_len(trait.colour, n_models), trait_names)

        ## --- Extract + combine marker data -----------------------------------
        mk_list <- lapply(seq_len(n_models), function(i) {
            df <- .resolve_markers(all_models[[i]])
            if (is.null(df)) return(NULL)
            df$trait <- trait_names[[i]]
            df
        })
        has_sig <- !vapply(mk_list, is.null, logical(1L))
        if (!any(has_sig)) {
            warning("No significant markers in any model — plotting map only.")
            return(invisible(do.call(linkMap, .lm_base_args(
                fake_cross, chr = if (missing(chr)) NULL else chr,
                row.chr = row.chr, nrow = nrow,
                marker.names = "none", m.cex = m.cex,
                pass_args = pass_args))))
        }
        if (!all(has_sig))
            warning("No significant markers in: ",
                    paste(trait_names[!has_sig], collapse = ", "), " — omitting.")

        marker_df        <- do.call(rbind, mk_list[has_sig])
        marker_df$trait  <- factor(marker_df$trait, levels = trait_names)

        ## --- Pool sig markers for annotation, determine chromosomes ----------
        flank_marks <- unique(marker_df$name)
        if (missing(chr))
            chr <- unique(marker_df$chr)[order(unique(marker_df$chr))]

        show_marks <- !is.null(marker.names) && !identical(marker.names, "none")
        fm         <- if (show_marks) flank_marks else NULL
        base_args  <- .lm_base_args(fake_cross,
                                     chr = if (is.null(row.chr)) chr else NULL,
                                     row.chr = row.chr, nrow = nrow,
                                     flanking_marks = fm,
                                     marker.names = marker.names, m.cex = m.cex,
                                     pass_args = pass_args)
        p <- do.call(linkMap, base_args)

        ## --- Restrict to displayed chromosomes --------------------------------
        oc <- marker_df$chr %in% chr
        if (!all(oc)) {
            warning("Some significant markers on chromosomes not shown — omitting.")
            marker_df <- marker_df[oc, , drop = FALSE]
        }
        if (nrow(marker_df) == 0L) return(invisible(p))
        marker_df$chr <- factor(marker_df$chr, levels = chr)

        ## --- Marker band layer (colour mapped → legend) ----------------------
        p <- p +
            ggplot2::geom_segment(
                data = marker_df,
                ggplot2::aes(x = -(K$TICK_X1 / 2), xend = (K$TICK_X1 / 2),
                             y = pos, yend = pos, colour = trait),
                linewidth = band.lwd, lineend = "butt", alpha = 0.78
            ) +
            ggplot2::scale_colour_manual(values = band_cols, name = "Trait")

        ## --- Left labels -----------------------------------------------------
        left_df <- data.frame(chr = marker_df$chr, pos = marker_df$pos,
                               label = as.character(marker_df$trait),
                               stringsAsFactors = FALSE)
        p <- .lm_left_labels(p, left_df, K, m.cex,
                              colour = band_cols[as.character(marker_df$trait)])

        ## --- Right labels (marker names) -------------------------------------
        if (show_marks) {
            right_df <- data.frame(
                chr   = marker_df$chr,
                pos   = marker_df$pos,
                label = if (identical(marker.names, "dist"))
                            as.character(round(marker_df$pos, 1))
                        else marker_df$name,
                stringsAsFactors = FALSE
            )
            p <- .lm_right_labels(p, right_df, K, m.cex,
                                   colour = band_cols[as.character(marker_df$trait)])
        }

        return(p + ggplot2::ggtitle("Genetic Map with Significant GWAS Markers"))
    }

    ## =========================================================================
    ## SINGLE-MODEL PATH
    ## =========================================================================

    if (!length(object$QTL$effects)) {
        warning("No significant GWAS markers — plotting map only.")
        base_args <- .lm_base_args(fake_cross,
                                    chr = if (missing(chr)) NULL else chr,
                                    row.chr = row.chr, nrow = nrow,
                                    marker.names = "none", m.cex = m.cex,
                                    pass_args = pass_args)
        return(invisible(do.call(linkMap, base_args)))
    }

    marker_df <- .resolve_markers(object)
    if (is.null(marker_df)) {
        warning("Could not resolve cM positions for any significant marker.")
        return(invisible(linkMap(fake_cross, ...)))
    }

    ## --- Chromosomes + base map ----------------------------------------------
    if (missing(chr))
        chr <- unique(marker_df$chr)[order(unique(marker_df$chr))]

    show_marks <- !is.null(marker.names) && !identical(marker.names, "none")
    fm         <- if (show_marks) unique(marker_df$name) else NULL
    base_args  <- .lm_base_args(fake_cross,
                                  chr = if (is.null(row.chr)) chr else NULL,
                                  row.chr = row.chr, nrow = nrow,
                                  flanking_marks = fm,
                                  marker.names = marker.names, m.cex = m.cex,
                                  pass_args = pass_args)
    p <- do.call(linkMap, base_args)

    ## --- Restrict to displayed chromosomes -----------------------------------
    oc <- marker_df$chr %in% chr
    if (!all(oc)) {
        warning("Some significant markers on chromosomes not shown — omitting.")
        marker_df <- marker_df[oc, , drop = FALSE]
    }
    if (nrow(marker_df) == 0L) return(invisible(p))
    marker_df$chr <- factor(marker_df$chr, levels = chr)

    ## --- Marker band layer ---------------------------------------------------
    p <- p +
        ggplot2::geom_segment(
            data      = marker_df,
            ggplot2::aes(x = -(K$TICK_X1 / 2), xend = (K$TICK_X1 / 2),
                         y = pos, yend = pos),
            linewidth = band.lwd, colour = marker.colour, lineend = "butt"
        )

    ## --- Left labels (trait name) --------------------------------------------
    trait_label <- if (is.null(trait.labels))
        as.character(object$call$fixed[[2L]])
    else
        trait.labels[[1L]]

    left_df <- data.frame(chr = marker_df$chr, pos = marker_df$pos,
                           label = trait_label, stringsAsFactors = FALSE)
    p <- .lm_left_labels(p, left_df, K, m.cex, colour = trait.colour)

    ## --- Right labels (marker names / distances) -----------------------------
    if (show_marks) {
        right_df <- data.frame(
            chr   = marker_df$chr,
            pos   = marker_df$pos,
            label = if (identical(marker.names, "dist"))
                        as.character(round(marker_df$pos, 1))
                    else marker_df$name,
            stringsAsFactors = FALSE
        )
        p <- .lm_right_labels(p, right_df, K, m.cex, colour = label.colour)
    }

    p + ggplot2::ggtitle("Genetic Map with Significant GWAS Markers")
}
