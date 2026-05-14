# =============================================================================
# plot.gwasAim.R
# S3 plot method for gwasAim objects.
#
# type = "manhattan" (default)
#   Classic Manhattan plot — points coloured by chromosome with alternating
#   shaded backgrounds. Significant markers highlighted and labelled.
#
# type = "outlier"
#   Line plot of genome-wide outlier statistics per iteration.
#   Uses the shared engine from plot.qtlAim.R.
#
# type = "blups"
#   Line plot of genome-wide scaled BLUPs per iteration.
#   Uses the shared engine from plot.qtlAim.R.
# =============================================================================

#' @describeIn gwasAim Produce a diagnostic or results plot for a
#'   \code{gwasAim} object. Five plot types are available via \code{type}.
#'
#'   \describe{
#'     \item{\code{"manhattan"} (default)}{Classic Manhattan plot with
#'       alternating shaded chromosome backgrounds. Points are coloured by
#'       chromosome; significant markers are highlighted in \code{sig.col}
#'       and labelled. One facet per requested iteration.}
#'     \item{\code{"outlier"}}{Line plot of genome-wide marker outlier
#'       statistics \eqn{\tilde{q}_i^2/\tilde{v}_i} per iteration.}
#'     \item{\code{"blups"}}{Line plot of genome-wide scaled BLUPs
#'       \eqn{\tilde{q}_i/\sqrt{|\tilde{v}_i|}} per iteration.}
#'     \item{\code{"effects"}}{Lollipop plot of estimated additive effect sizes
#'       for each significant marker, with \eqn{\pm} 1 SE error bars and
#'       percentage of phenotypic variance annotated. Positive effect = alt
#'       allele (dosage 2) favoured; negative = ref allele (dosage 0) favoured.}
#'     \item{\code{"contrast"}}{Allele contrast plot showing total line genetic
#'       values split by 0/1/2 dosage class at each significant marker. Requires
#'       \code{data}.}
#'     \item{\code{"heatmap"}}{Genome × iteration heatmap: tile fill encodes
#'       the outlier statistic at every marker position across all
#'       forward-selection iterations. Excluded positions are shown in grey.
#'       Detected significant markers are annotated with a filled diamond.
#'       Iteration 1 is at the top.}
#'   }
#'
#' @param x A \code{gwasAim} object.
#' @param genObj The \code{"wgPanel"} object passed to \code{gwasAim}, produced
#'   by \code{\link{primePanel}}.
#' @param type Character string: \code{"manhattan"} (default), \code{"outlier"},
#'   \code{"blups"}, \code{"effects"}, or \code{"contrast"}.
#' @param data Data frame. The phenotypic data used in the analysis. Required
#'   only when \code{type = "contrast"}.
#' @param iter Integer vector of iterations to display. Default is all stored
#'   iterations. Not used for \code{"effects"} or \code{"contrast"}.
#' @param chr Optional character vector of chromosome names to display. Not
#'   used for \code{"effects"} or \code{"contrast"}.
#' @param chr.lines Logical. Draw vertical lines at chromosome boundaries in
#'   \code{"outlier"} and \code{"blups"} types. Default \code{FALSE}.
#' @param sig.col Colour for significant marker highlighting and labels.
#'   Default \code{"red"}.
#' @param pt.col Character vector of length 2. Alternating point/line colours
#'   for chromosomes. Default \code{c("steelblue","grey60")}.
#' @param band.col Character vector of length 2. Alternating background band
#'   colours used in the Manhattan plot. Default \code{c("grey92","white")}.
#' @param pt.cex Numeric point/line size. Default \code{0.6}.
#' @param ncol Integer. Number of columns in the \code{facet_wrap} layout for
#'   \code{type = "contrast"}. Default \code{1} (vertical stack). Set to
#'   \code{NULL} for automatic layout or any positive integer for a grid.
#'   Not used for other plot types.
#' @param cap Numeric scalar. Upper limit for the colour scale in
#'   \code{type = "heatmap"}. Outlier statistics above this value are clamped
#'   to \code{cap} and rendered in the maximum colour. Default \code{NULL}
#'   (no capping; full data range used). Not used for other plot types.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
plot.gwasAim <- function(x, genObj,
                          type      = c("manhattan", "outlier", "blups",
                                        "effects", "contrast", "heatmap"),
                          data      = NULL,
                          ncol      = 1L,
                          iter      = NULL,
                          chr       = NULL,
                          chr.lines = FALSE,
                          sig.col   = "red",
                          pt.col    = c("steelblue", "grey60"),
                          band.col  = c("grey92", "white"),
                          pt.cex    = 0.6,
                          cap       = NULL,
                          ...) {

    type <- match.arg(type)

    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgPanel"))
        stop("genObj must be of class \"wgPanel\" produced by primePanel().")

    if (is.null(x$QTL$effects) && type %in% c("effects", "contrast"))
        stop("No significant markers found in object.")

    # -------------------------------------------------------------------------
    # effects / contrast: delegate to shared helpers in plot.qtlAim.R
    # -------------------------------------------------------------------------
    if (type == "effects") {
        edf <- .build_effects_df(x, genObj)
        return(.plot_effects(edf))
    }

    if (type == "contrast") {
        if (is.null(data))
            stop("data is required for type = \"contrast\".\n",
                 "Pass the phenotypic data frame used in the analysis ",
                 "(e.g. the <response>.data object).")
        cdf <- .build_contrast_df(x, genObj, data)
        return(.plot_contrast(cdf, ncol = ncol))
    }

    # -------------------------------------------------------------------------
    # Resolve iterations
    # -------------------------------------------------------------------------
    n.stored <- length(x$QTL$diag$oint)
    if (is.null(iter))
        iter <- seq_len(n.stored)
    else if (any(iter > n.stored))
        stop("iter contains values greater than number of stored iterations (",
             n.stored, ").")

    # -------------------------------------------------------------------------
    # Chromosome subset
    # -------------------------------------------------------------------------
    all.chr <- names(genObj$geno)
    if (!is.null(chr)) {
        bad <- chr[!chr %in% all.chr]
        if (length(bad))
            stop("Chromosome(s) not found in genObj: ",
                 paste(bad, collapse = ", "), ".")
    } else {
        chr <- all.chr
    }

    # -------------------------------------------------------------------------
    # Dispatch
    # -------------------------------------------------------------------------

    # heatmap: genome x iteration outlier heatmap (shared engine)
    if (type == "heatmap")
        return(.plot_heatmap(x, genObj, iter, chr, sig.col, cap))

    cp <- .build_cumpos(genObj, "marker", chr)

    if (type == "manhattan")
        return(.plot_gwas_manhattan(x, genObj, iter, chr, cp,
                                    sig.col, pt.col, band.col, pt.cex))

    # outlier / blups: shared line-plot engine
    stat_slot <- if (type == "outlier") "oint" else "blups"
    y.lab     <- if (type == "outlier") "Outlier Statistic" else "Scaled BLUPs"

    stat_df <- .build_stat_df(x, stat_slot, iter, chr, cp)
    pt.col.map        <- stats::setNames(
        rep_len(pt.col, length(chr)), chr)
    stat_df$col <- pt.col.map[stat_df$chr]

    gp <- ggplot2::ggplot(stat_df,
               ggplot2::aes(x = .data$dist, y = .data$values,
                            colour = .data$col, group = .data$chr)) +
        ggplot2::facet_wrap(~iteration, ncol = 1, scales = "free_y") +
        ggplot2::geom_line(linewidth = pt.cex) +
        ggplot2::geom_rug(sides = "b",
                          length = ggplot2::unit(0.02, "npc"),
                          linewidth = 0.3) +
        ggplot2::scale_colour_identity() +
        ggplot2::scale_x_continuous(
            breaks = cp$chr.mid, labels = names(cp$chr.mid)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
        ggplot2::ylab(y.lab) +
        ggplot2::xlab("") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter()

    if (chr.lines && length(chr) > 1)
        gp <- gp + ggplot2::geom_vline(
            xintercept = cp$chr.end[-length(cp$chr.end)],
            colour = "grey80", linewidth = 0.5, linetype = "dashed")

    .add_sig_labels(gp, x, stat_slot, iter, chr, cp, sig.col)
}

# =============================================================================
# Internal: Manhattan plot with alternating chromosome backgrounds
# =============================================================================
.plot_gwas_manhattan <- function(object, genObj, iter, chr, cp,
                                  sig.col, pt.col, band.col, pt.cex) {

    stat_df <- .build_stat_df(object, "oint", iter, chr, cp)

    # Alternating point colours by chromosome
    chr.idx        <- match(stat_df$chr, chr)
    stat_df$col    <- pt.col[(chr.idx %% 2) + 1L]

    # -------------------------------------------------------------------------
    # Chromosome background bands
    # -------------------------------------------------------------------------
    band.df <- data.frame(
        xmin  = c(0, cp$chr.end[-length(cp$chr.end)] + 5),
        xmax  = cp$chr.end,
        fill  = band.col[(seq_along(chr) %% 2) + 1L],
        stringsAsFactors = FALSE
    )

    gp <- ggplot2::ggplot(stat_df,
               ggplot2::aes(x = .data$dist, y = .data$values)) +
        ggplot2::facet_wrap(~iteration, ncol = 1, scales = "free_y") +
        # Background bands — drawn first, behind points
        ggplot2::geom_rect(
            data        = band.df,
            ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                         ymin = -Inf,        ymax = Inf,
                         fill = .data$fill),
            inherit.aes = FALSE,
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_identity() +
        # Points
        ggplot2::geom_point(
            ggplot2::aes(colour = .data$col),
            size = pt.cex, show.legend = FALSE
        ) +
        ggplot2::scale_colour_identity() +
        ggplot2::scale_x_continuous(
            breaks = cp$chr.mid, labels = names(cp$chr.mid)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
        ggplot2::ylab("Outlier Statistic") +
        ggplot2::xlab("Chromosome") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter() +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.border     = ggplot2::element_rect(fill = NA,      colour = "grey70")
        )

    # -------------------------------------------------------------------------
    # Highlight significant markers per iteration facet
    # -------------------------------------------------------------------------
    qtl.keys <- object$QTL$qtl
    if (!is.null(qtl.keys) && length(qtl.keys) > 0L) {
        iter.sig <- iter[iter <= length(qtl.keys)]
        sig.rows <- lapply(iter.sig, function(it) {
            key  <- qtl.keys[it]
            parts <- strsplit(key, "\\.")[[1L]]
            qchr  <- parts[2L]
            if (!qchr %in% chr) return(NULL)
            idx   <- as.integer(parts[3L])
            mkr   <- names(genObj$geno[[qchr]]$map)[idx]
            val   <- object$QTL$diag$oint[[it]][key]
            if (is.null(val) || is.na(val)) return(NULL)
            data.frame(
                dist      = cp$pos_lookup[key],
                values    = as.numeric(val),
                iteration = factor(paste0("Iteration: ", it),
                                   levels = levels(stat_df$iteration)),
                label     = mkr,
                stringsAsFactors = FALSE
            )
        })
        sig.df <- do.call(rbind, sig.rows[!sapply(sig.rows, is.null)])
        if (!is.null(sig.df) && nrow(sig.df) > 0L) {
            gp <- gp +
                ggplot2::geom_point(
                    data        = sig.df,
                    ggplot2::aes(x = .data$dist, y = .data$values),
                    colour      = sig.col, size = pt.cex * 5,
                    shape       = 18, inherit.aes = FALSE
                ) +
                ggplot2::geom_text(
                    data        = sig.df,
                    ggplot2::aes(x = .data$dist, y = .data$values,
                                 label = .data$label),
                    colour      = sig.col, size = 2.8,
                    vjust       = -0.8, hjust = 0.5,
                    inherit.aes = FALSE
                )
        }
    }
    gp
}
