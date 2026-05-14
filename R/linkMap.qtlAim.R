# =============================================================================
# linkMap.qtlAim.R
# S3 linkMap() method for qtlAim objects — single model OR multiple models.
#
# Multiple models are passed as extra positional arguments in ...:
#   linkMap(qtl_yld, qtl_tgw, genObj = my_int)
#
# Any unnamed ... arg that inherits "qtlAim" is treated as an additional model.
# Named ... args (e.g. ggrepel options) are passed through to linkMap.cross().
#
# The old list-based API still works via linkMap.default (backward compat):
#   linkMap(list(qtl_yld, qtl_tgw), genObj = my_int)
#
# genObj must be a "wgCross" object produced by primeCross().
#
# Shared helpers (.lm_base_args, .lm_left_labels, .lm_right_labels) live in
# linkMap.cross.R.
# =============================================================================

#' @importFrom ggrepel geom_text_repel
#' @describeIn qtlAim Plot the genetic linkage map with detected QTL overlaid.
#'   Accepts a single \code{qtlAim} object or multiple models passed as extra
#'   positional arguments — e.g.
#'   \code{linkMap(qtl_yld, qtl_tgw, genObj = my_int)} — for a
#'   multi-trait overlay.  QTL intervals are drawn as coloured fills inside
#'   the chromosome bar; flanking markers are annotated on the right; trait
#'   labels sit to the left.  Returns a \code{ggplot} object.
#' @param object A \code{qtlAim} object, or the first of several when multiple
#'   models are passed.
#' @param genObj The \code{"wgCross"} object produced by \code{primeCross()}
#'   used in the analysis.
#' @param chr Optional character vector of chromosome names to display.
#' @param chr.dist Optional named list with \code{$start} / \code{$end}
#'   restricting the displayed cM range.
#' @param marker.names \code{"markers"} (default), \code{"dist"},
#'   \code{"none"} or \code{NULL}.
#' @param flanking Logical; if \code{TRUE} (default) only QTL flanking markers
#'   are annotated.
#' @param qtl.colour Colour(s) for QTL fills.  Recycled across traits.
#' @param marker.colour Colour for flanking marker labels.
#' @param trait.colour Colour(s) for left-side trait labels.  Defaults to
#'   \code{qtl.colour}.
#' @param m.cex Text size (pt) for labels.  Default \code{7}.
#' @param trait.labels Optional character vector of trait names (one per
#'   model).
#' @param tick Logical.  Unused; retained for back-compatibility.
#' @param nrow Integer or \code{"auto"}.  Facet row count.
#' @param row.chr Optional list of character vectors for manual row layout.
#' @param \dots Additional \code{qtlAim} models (positional, unnamed) passed
#'   for multi-trait display, and/or further arguments to
#'   \code{\link{linkMap.cross}}.
#' @return A \code{ggplot} object.
#' @export
linkMap.qtlAim <- function(object, ...,
                            genObj,
                            chr, chr.dist,
                            marker.names  = "markers",
                            flanking      = TRUE,
                            qtl.colour    = "steelblue",
                            marker.colour = "firebrick",
                            trait.colour  = qtl.colour,
                            m.cex         = 7,
                            trait.labels  = NULL,
                            tick          = FALSE,
                            nrow          = 1,
                            row.chr       = NULL) {

    K <- .lm_const()

    ## -------------------------------------------------------------------------
    ## Separate extra qtlAim models from pass-through ... args
    ## -------------------------------------------------------------------------
    dots  <- list(...)
    dnames <- if (is.null(names(dots))) rep("", length(dots)) else
              ifelse(is.na(names(dots)), "", names(dots))
    is_extra    <- vapply(dots, inherits, logical(1L), "qtlAim") & dnames == ""
    extra_models <- dots[is_extra]
    pass_args    <- dots[!is_extra]

    ## Convert chr.dist missing-ness to NULL for helpers
    chr.dist_val <- if (missing(chr.dist)) NULL else chr.dist

    ## =========================================================================
    ## MULTI-MODEL PATH
    ## =========================================================================
    if (length(extra_models) > 0L) {

        all_models <- c(list(object), extra_models)
        n_models   <- length(all_models)

        ## --- Validate --------------------------------------------------------
        if (missing(genObj))
            stop("'genObj' is a required argument.")
        if (!inherits(genObj, "wgCross"))
            stop("'genObj' must be of class \"wgCross\".")

        types <- vapply(all_models, function(m) m$QTL$type, character(1L))
        if (length(unique(types)) > 1L)
            stop("All models must have the same QTL type.")
        qtl_type <- types[[1L]]

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
        qtl_cols   <- stats::setNames(rep_len(qtl.colour,   n_models), trait_names)
        trait_cols <- stats::setNames(rep_len(trait.colour, n_models), trait_names)

        ## --- Extract + combine QTL data --------------------------------------
        qtl_list <- lapply(seq_len(n_models), function(i) {
            mod <- all_models[[i]]
            if (!length(mod$QTL$effects)) return(NULL)
            qtlm <- getQTL(mod, genObj)
            wchr <- qtlm[, 1L]
            if (qtl_type == "interval") {
                data.frame(chr = wchr, lh.mark = qtlm[, 3L],
                           lh.pos = as.numeric(qtlm[, 4L]),
                           rh.mark = qtlm[, 7L],
                           rh.pos  = as.numeric(qtlm[, 8L]),
                           trait = trait_names[[i]],
                           stringsAsFactors = FALSE)
            } else {
                data.frame(chr = wchr, lh.mark = qtlm[, 3L],
                           lh.pos = as.numeric(qtlm[, 4L]),
                           rh.mark = qtlm[, 3L],
                           rh.pos  = as.numeric(qtlm[, 4L]),
                           trait = trait_names[[i]],
                           stringsAsFactors = FALSE)
            }
        })

        has_qtl <- !vapply(qtl_list, is.null, logical(1L))
        if (!any(has_qtl)) {
            warning("No QTL found in any model — plotting map only.")
            return(invisible(linkMap(genObj, ...)))
        }
        if (!any(has_qtl))
            warning("No QTL in: ",
                    paste(trait_names[!has_qtl], collapse = ", "), " — omitting.")

        qtl_df      <- do.call(rbind, qtl_list[has_qtl])
        qtl_df$trait <- factor(qtl_df$trait, levels = trait_names)

        ## --- Flanking markers, chromosomes -----------------------------------
        flank_marks <- if (qtl_type == "interval")
            unique(c(qtl_df$lh.mark, qtl_df$rh.mark))
        else
            unique(qtl_df$lh.mark)

        if (missing(chr))
            chr <- unique(qtl_df$chr)[order(unique(qtl_df$chr))]

        ## --- Base map --------------------------------------------------------
        show_marks <- !is.null(marker.names) && !identical(marker.names, "none")
        fm         <- if (flanking && show_marks) flank_marks else NULL
        base_args  <- .lm_base_args(genObj, chr = if (is.null(row.chr)) chr else NULL,
                                     row.chr = row.chr, nrow = nrow,
                                     flanking_marks = fm,
                                     marker.names = marker.names, m.cex = m.cex,
                                     chr.dist = chr.dist_val,
                                     pass_args = pass_args)
        p <- do.call(linkMap, base_args)

        ## --- Restrict to displayed chromosomes -------------------------------
        oc <- qtl_df$chr %in% chr
        if (!all(oc)) {
            warning("Some QTL on chromosomes not shown — omitting.")
            qtl_df <- qtl_df[oc, , drop = FALSE]
        }
        if (nrow(qtl_df) == 0L) return(invisible(p))
        qtl_df$chr   <- factor(qtl_df$chr, levels = chr)
        qtl_df$fill_col <- qtl_cols[as.character(qtl_df$trait)]

        ## --- QTL layer (colour mapped → legend) ------------------------------
        if (qtl_type == "interval") {
            p <- p +
                ggplot2::geom_segment(
                    data = qtl_df,
                    ggplot2::aes(x = 0, xend = 0, y = lh.pos, yend = rh.pos,
                                 colour = trait),
                    linewidth = K$BAR_FILL_LW, lineend = "butt", alpha = 0.78
                )
        } else {
            p <- p +
                ggplot2::geom_segment(
                    data = qtl_df,
                    ggplot2::aes(x = -(K$TICK_X1 / 2), xend = (K$TICK_X1 / 2),
                                 y = lh.pos, yend = lh.pos, colour = trait),
                    linewidth = 1.5, alpha = 0.78
                )
        }
        p <- p + ggplot2::scale_colour_manual(values = qtl_cols, name = "Trait")

        ## --- Left labels -----------------------------------------------------
        left_df <- data.frame(
            chr   = qtl_df$chr,
            pos   = (qtl_df$lh.pos + qtl_df$rh.pos) / 2,
            label = as.character(qtl_df$trait),
            stringsAsFactors = FALSE
        )
        p <- .lm_left_labels(p, left_df, K, m.cex,
                              colour = qtl_cols[as.character(qtl_df$trait)])

        return(p + ggplot2::ggtitle("Genetic Map with QTL"))
    }

    ## =========================================================================
    ## SINGLE-MODEL PATH
    ## =========================================================================

    if (missing(genObj))
        stop("'genObj' is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("'genObj' must be of class \"wgCross\".")
    if (!length(object$QTL$effects)) {
        warning("No significant QTL detected — plotting map only.")
        return(invisible(linkMap(genObj, ...)))
    }

    ## --- Extract QTL positions -----------------------------------------------
    qtlm <- getQTL(object, genObj)
    wchr <- qtlm[, 1L]

    if (object$QTL$type == "interval") {
        qtl_df <- data.frame(
            chr     = wchr, lh.mark = qtlm[, 3L],
            lh.pos  = as.numeric(qtlm[, 4L]), rh.mark = qtlm[, 7L],
            rh.pos  = as.numeric(qtlm[, 8L]), stringsAsFactors = FALSE
        )
        flank_marks <- unique(c(qtlm[, 3L], qtlm[, 7L]))
    } else {
        qtl_df <- data.frame(
            chr     = wchr, lh.mark = qtlm[, 3L],
            lh.pos  = as.numeric(qtlm[, 4L]), rh.mark = qtlm[, 3L],
            rh.pos  = as.numeric(qtlm[, 4L]), stringsAsFactors = FALSE
        )
        flank_marks <- unique(qtlm[, 3L])
    }

    ## --- Trait name + colours ------------------------------------------------
    n_qtl     <- nrow(qtl_df)
    trait_vec <- if (is.null(trait.labels))
        rep(as.character(object$call$fixed[[2L]]), n_qtl)
    else if (length(trait.labels) == 1L)
        rep(trait.labels, n_qtl)
    else
        trait.labels

    qtl_df$trait   <- factor(trait_vec, levels = unique(trait_vec))
    qtrait         <- levels(qtl_df$trait)
    n_traits       <- length(qtrait)
    qtl_cols       <- stats::setNames(rep_len(qtl.colour,   n_traits), qtrait)
    trait_cols     <- stats::setNames(rep_len(trait.colour, n_traits), qtrait)

    ## --- Chromosomes + base map ----------------------------------------------
    if (missing(chr))
        chr <- unique(wchr)[order(unique(wchr))]

    show_marks <- !is.null(marker.names) && !identical(marker.names, "none")
    fm         <- if (flanking && show_marks) flank_marks else NULL
    base_args  <- .lm_base_args(genObj,
                                  chr = if (is.null(row.chr)) chr else NULL,
                                  row.chr = row.chr, nrow = nrow,
                                  flanking_marks = fm,
                                  marker.names = marker.names, m.cex = m.cex,
                                  chr.dist = chr.dist_val,
                                  pass_args = pass_args)
    p <- do.call(linkMap, base_args)

    ## --- Restrict QTL to displayed chromosomes --------------------------------
    oc <- qtl_df$chr %in% chr
    if (!all(oc)) {
        warning("Some QTL on chromosomes not shown — omitting.")
        qtl_df <- qtl_df[oc, , drop = FALSE]
    }
    if (nrow(qtl_df) == 0L) return(invisible(p))
    qtl_df$chr      <- factor(qtl_df$chr, levels = chr)
    qtl_df$fill_col <- qtl_cols[as.character(qtl_df$trait)]

    ## --- QTL layer -----------------------------------------------------------
    if (object$QTL$type == "interval") {
        p <- p +
            ggplot2::geom_segment(
                data = qtl_df,
                ggplot2::aes(x = 0, xend = 0, y = lh.pos, yend = rh.pos),
                linewidth = K$BAR_FILL_LW, lineend = "butt",
                colour = qtl_df$fill_col, alpha = 0.85
            )
    } else {
        p <- p +
            ggplot2::geom_segment(
                data = qtl_df,
                ggplot2::aes(x = -(K$TICK_X1 / 2), xend = (K$TICK_X1 / 2),
                             y = lh.pos, yend = lh.pos),
                linewidth = 1.5, colour = qtl_df$fill_col
            )
    }

    ## --- Left labels (trait names) -------------------------------------------
    left_df <- data.frame(
        chr   = qtl_df$chr,
        pos   = (qtl_df$lh.pos + qtl_df$rh.pos) / 2,
        label = as.character(qtl_df$trait),
        stringsAsFactors = FALSE
    )
    p <- .lm_left_labels(p, left_df, K, m.cex,
                          colour = trait_cols[as.character(qtl_df$trait)])

    p + ggplot2::ggtitle("Genetic Map with QTL")
}
