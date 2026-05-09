# =============================================================================
# manhattan.R
# S3 generic and GWASAim method for Manhattan plots.
#
# Plots per-iteration outlier statistics across the genome, faceted by
# iteration. The marker selected as significant at iteration k is labelled
# in facet k only. No threshold line is drawn — significant markers have
# already been determined by the forward selection algorithm.
# =============================================================================

#' Manhattan Plot for GWASAim Objects
#'
#' @description
#' Produces a genome-wide Manhattan plot of outlier statistics from a
#' \code{\link{GWASAim}} forward-selection analysis, with one facet per
#' requested iteration. Each facet shows the outlier statistic for all active
#' markers at that iteration, with the marker selected as significant at that
#' iteration highlighted and labelled. No significance threshold line is drawn:
#' significant markers have already been determined by the forward-selection
#' likelihood ratio test.
#'
#' @param object A \code{GWASAim} object.
#' @param \dots Additional arguments passed to methods.
#'
#' @name manhattan
#' @export
manhattan <- function(object, ...)
    UseMethod("manhattan")

#' @describeIn manhattan Manhattan plot for a \code{GWASAim} object.
#' @param object A \code{GWASAim} object.
#' @param panelObj The \code{"panel"} object passed to \code{GWASAim}, used
#'   to look up marker cM positions for the x-axis.
#' @param iter Integer vector of iterations to display. Default is all stored
#'   iterations (\code{NULL}). Each iteration produces one facet panel.
#' @param chr Optional character vector of chromosome names to display. Default
#'   is all chromosomes.
#' @param pt.col Character vector of length 2 giving alternating chromosome
#'   point colours. Default is \code{c("steelblue","grey50")}.
#' @param sig.col Colour used to highlight and label the significant marker
#'   in each facet. Default is \code{"red"}.
#' @param pt.cex Numeric point size. Default is \code{0.6}.
#' @return A \code{\link[ggplot2]{ggplot}} object. The y-axis shows the raw
#'   outlier statistic \eqn{\tilde{q}_i^2 / \tilde{v}_i} (not a transformed
#'   p-value). The selected marker in each iteration is shown as a larger
#'   coloured diamond with its marker name as a label.
#' @export
manhattan.GWASAim <- function(object, panelObj,
                               iter    = NULL,
                               chr     = NULL,
                               pt.col  = c("steelblue", "grey50"),
                               sig.col = "red",
                               pt.cex  = 0.6, ...) {

    if (missing(panelObj))
        stop("panelObj is a required argument.")
    if (!inherits(panelObj, "panel"))
        stop("panelObj must be of class \"panel\".")

    # Default: all stored iterations
    n.stored <- length(object$QTL$diag$oint)
    if (is.null(iter))
        iter <- seq_len(n.stored)
    if (any(iter > n.stored))
        stop("iter contains values exceeding the number of stored iterations (", n.stored, ").")

    # -------------------------------------------------------------------------
    # Chromosome subset
    # -------------------------------------------------------------------------
    chr.names <- names(panelObj$geno)
    if (!is.null(chr)) {
        if (!all(chr %in% chr.names))
            stop("Some chromosome names not found in panelObj.")
        chr.names <- chr
    }

    # -------------------------------------------------------------------------
    # Build cumulative position lookup: internal name -> cumulative cM
    # Internal names are "Chr.{chrName}.{markerIndex}" (from .buildGenoData)
    # -------------------------------------------------------------------------
    pos.lookup  <- c()
    chr.mid.pos <- c()
    cum.offset  <- 0
    gap.cM      <- 5    # gap between chromosomes on x-axis

    for (ch in chr.names) {
        map.pos <- as.numeric(panelObj$geno[[ch]]$map)
        n.mk    <- length(map.pos)
        span    <- map.pos[n.mk] - map.pos[1]
        for (k in seq_len(n.mk)) {
            key            <- paste("Chr", ch, k, sep = ".")
            pos.lookup[key] <- cum.offset + (map.pos[k] - map.pos[1])
        }
        chr.mid.pos[ch] <- cum.offset + span / 2
        cum.offset      <- cum.offset + span + gap.cM
    }

    # -------------------------------------------------------------------------
    # Build long-format data frame across all requested iterations
    # -------------------------------------------------------------------------
    c.iter <- paste("Iteration:", iter)

    oint.list <- object$QTL$diag$oint[iter]
    names(oint.list) <- c.iter

    ointl <- lapply(oint.list, function(el) {
        # Filter to requested chromosomes and non-zero entries
        echr <- sapply(strsplit(names(el), "\\."), "[", 2)
        whc  <- echr %in% chr.names & el > 0
        data.frame(
            values = as.numeric(el[whc]),
            chr    = echr[whc],
            iname  = names(el)[whc],          # internal "Chr.chrName.idx"
            dist   = pos.lookup[names(el)[whc]],
            stringsAsFactors = FALSE
        )
    })

    char.iter <- factor(
        rep(c.iter, times = sapply(ointl, nrow)),
        levels = unique(c.iter)
    )

    ointd      <- do.call("rbind.data.frame", ointl)
    ointd$iteration <- char.iter
    # Alternating chromosome colours
    ointd$col  <- pt.col[(match(ointd$chr, chr.names) %% 2) + 1]

    # -------------------------------------------------------------------------
    # Significant marker per iteration
    # Each facet k labels ONLY the marker selected at that iteration.
    # object$QTL$qtl[k] is the internal name of the marker selected at iter k.
    # The last iteration has no selection (it failed the significance test).
    # -------------------------------------------------------------------------
    qtl.names <- object$QTL$qtl   # may be NULL if nothing found

    sig.df <- NULL
    if (!is.null(qtl.names) && length(qtl.names) > 0) {
        # Only iterations that actually produced a significant selection
        iter.with.sig <- iter[iter <= length(qtl.names)]

        sig.rows <- lapply(iter.with.sig, function(it) {
            qn   <- qtl.names[it]
            echr <- strsplit(qn, "\\.")[[1]][2]
            if (!echr %in% chr.names) return(NULL)
            idx  <- as.integer(strsplit(qn, "\\.")[[1]][3])
            mkr  <- names(panelObj$geno[[echr]]$map)[idx]
            val  <- object$QTL$diag$oint[[it]][qn]
            if (is.null(val) || is.na(val)) return(NULL)
            data.frame(
                values    = as.numeric(val),
                chr       = echr,
                iname     = qn,
                dist      = pos.lookup[qn],
                iteration = factor(paste("Iteration:", it), levels = levels(char.iter)),
                label     = mkr,
                stringsAsFactors = FALSE
            )
        })
        sig.df <- do.call("rbind", sig.rows[!sapply(sig.rows, is.null)])
    }

    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------
    gp <- ggplot(ointd, aes_string(x = "dist", y = "values", colour = "col")) +
        facet_wrap(~iteration, ncol = 1, scales = "free_y") +
        geom_point(size = pt.cex, show.legend = FALSE) +
        scale_colour_identity() +
        scale_x_continuous(breaks = chr.mid.pos, labels = chr.names) +
        labs(x = "Chromosome", y = "Outlier Statistic") +
        theme_scatter()

    # Overlay significant marker point + label in the matching facet only
    if (!is.null(sig.df) && nrow(sig.df) > 0) {
        gp <- gp +
            geom_point(data = sig.df,
                       aes_string(x = "dist", y = "values"),
                       colour = sig.col, size = pt.cex * 4,
                       shape = 18, inherit.aes = FALSE) +
            geom_text(data = sig.df,
                      aes_string(x = "dist", y = "values", label = "label"),
                      colour = sig.col, size = 2.8,
                      vjust = -1.0, hjust = 0.5, inherit.aes = FALSE)
    }
    gp
}
