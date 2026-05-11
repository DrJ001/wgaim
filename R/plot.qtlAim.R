# =============================================================================
# plot.qtlAim.R
# S3 plot method for qtlAim objects.
#
# Plots the forward-selection outlier statistics or scaled BLUPs across the
# genome for each requested iteration, with detected QTL labelled.
#
# Shared internal engine functions used also by plot.gwasAim:
#   .build_cumpos()    — cumulative cM position lookup from a genObj
#   .build_stat_df()   — long-format data frame from diag oint/blups vectors
#   .add_sig_labels()  — annotate selected QTL/markers onto a ggplot
#   .plot_chr_bars()   — chromosome-level bar chart (plot.chr mode)
# =============================================================================

#' @describeIn qtlAim Plot genome-wide outlier statistics or scaled BLUPs from
#'   the forward-selection analysis, with one facet per iteration. Detected QTL
#'   are annotated in the facet where they were selected.
#'
#' @param x A \code{qtlAim} object.
#' @param genObj The \code{"wgCross"} object passed to \code{qtlAim}, produced
#'   by \code{\link{primeCross}}.
#' @param type Character string. \code{"outlier"} (default) plots the
#'   interval/marker outlier statistics \eqn{\tilde{q}_i^2/\tilde{v}_i}.
#'   \code{"blups"} plots the scaled BLUPs
#'   \eqn{\tilde{q}_i / \sqrt{|\tilde{v}_i|}}.
#'   \code{"chr"} plots chromosome-level aggregated outlier statistics as a
#'   bar chart (only available when \code{qtlAim} was run with
#'   \code{selection = "chromosome"}).
#' @param iter Integer vector of iterations to display. Default is all stored
#'   iterations.
#' @param chr Optional character vector of chromosome names to display.
#' @param chr.lines Logical. If \code{TRUE}, vertical lines are drawn at
#'   chromosome boundaries. Default \code{FALSE}.
#' @param sig.col Colour used to highlight and label the selected QTL in each
#'   facet. Default \code{"firebrick"}.
#' @param pt.col Character vector of length equal to the number of displayed
#'   chromosomes, or length 2 for alternating colours. Colours the line/rug
#'   traces by chromosome. Default uses \code{hcl.colors}.
#' @param pt.cex Numeric line width. Default \code{0.6}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
plot.qtlAim <- function(x, genObj,
                         type      = c("outlier", "blups", "chr"),
                         iter      = NULL,
                         chr       = NULL,
                         chr.lines = FALSE,
                         sig.col   = "firebrick",
                         pt.col    = NULL,
                         pt.cex    = 0.6,
                         ...) {

    type <- match.arg(type)

    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("genObj must be of class \"wgCross\" produced by primeCross().")

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
    # chr type: chromosome-level bar chart
    # -------------------------------------------------------------------------
    if (type == "chr") {
        if (is.null(x$QTL$diag$ochr[[1]]))
            stop("No chromosome outlier statistics found. ",
                 "Re-run qtlAim with selection = \"chromosome\".")
        return(.plot_chr_bars(x, iter, chr, sig.col))
    }

    # -------------------------------------------------------------------------
    # outlier / blups: line plot
    # -------------------------------------------------------------------------
    gen.type <- x$QTL$type
    cp       <- .build_cumpos(genObj, gen.type, chr)
    stat_slot <- if (type == "outlier") "oint" else "blups"
    y.lab     <- if (type == "outlier") "Outlier Statistic" else "Scaled BLUPs"

    stat_df  <- .build_stat_df(x, stat_slot, iter, chr, cp)
    if (is.null(pt.col))
        pt.col <- grDevices::hcl.colors(length(chr), "Dark 2")
    pt.col.map <- stats::setNames(rep_len(pt.col, length(chr)), chr)
    stat_df$col <- pt.col.map[stat_df$chr]

    c.iter <- paste0("Iteration: ", iter)
    char.iter <- factor(rep(c.iter, times = tapply(stat_df$iteration,
                                                    stat_df$iteration, length)[c.iter]),
                        levels = c.iter)

    gp <- ggplot2::ggplot(stat_df,
               ggplot2::aes(x = .data$dist, y = .data$values,
                            colour = .data$col, group = .data$chr)) +
        ggplot2::facet_wrap(~iteration, ncol = 1, scales = "free_y") +
        ggplot2::geom_line(linewidth = pt.cex) +
        ggplot2::geom_rug(sides = "b",
                          length = ggplot2::unit(0.02, "npc"),
                          linewidth = 0.3) +
        ggplot2::scale_colour_identity() +
        ggplot2::scale_x_continuous(breaks = cp$chr.mid, labels = names(cp$chr.mid)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
        ggplot2::ylab(y.lab) +
        ggplot2::xlab("") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter()

    # Chromosome boundary lines
    if (chr.lines && length(chr) > 1)
        gp <- gp + ggplot2::geom_vline(
            xintercept = cp$chr.end[-length(cp$chr.end)],
            colour = "grey80", linewidth = 0.5, linetype = "dashed")

    # Annotate selected QTL
    gp <- .add_sig_labels(gp, x, stat_slot, iter, chr, cp, sig.col)
    gp
}

# =============================================================================
# Shared engine: .build_cumpos
# Build cumulative cM position lookup for a genObj.
#
# Returns a list:
#   $pos_lookup  — named vector: internal key "Chr.CHR.IDX" -> cumulative cM
#   $chr.mid     — named vector: chr name -> mid-point cumulative cM
#   $chr.end     — named vector: chr name -> end cumulative cM (for chr.lines)
# =============================================================================
.build_cumpos <- function(genObj, gen.type, chr) {
    gap.cM      <- 5
    pos.lookup  <- c()
    chr.mid     <- c()
    chr.end     <- c()
    cum.offset  <- 0

    for (ch in chr) {
        if (gen.type == "interval" && !is.null(genObj$geno[[ch]]$inferred.map))
            map.pos <- genObj$geno[[ch]]$inferred.map
        else
            map.pos <- genObj$geno[[ch]]$map

        n.mk <- length(map.pos)
        if (n.mk == 0) next
        span <- map.pos[n.mk] - map.pos[1]

        for (k in seq_len(n.mk)) {
            key             <- paste("Chr", ch, k, sep = ".")
            pos.lookup[key] <- cum.offset + (map.pos[k] - map.pos[1])
        }
        chr.mid[ch] <- cum.offset + span / 2
        chr.end[ch] <- cum.offset + span
        cum.offset  <- cum.offset + span + gap.cM
    }
    list(pos_lookup = pos.lookup, chr.mid = chr.mid, chr.end = chr.end)
}

# =============================================================================
# Shared engine: .build_stat_df
# Build a long-format data frame from the oint or blups diagnostic vectors.
# =============================================================================
.build_stat_df <- function(object, stat_slot, iter, chr, cp) {
    diag_data <- object$QTL$diag[[stat_slot]][iter]
    c.iter    <- paste0("Iteration: ", iter)
    names(diag_data) <- c.iter

    rows <- lapply(seq_along(diag_data), function(i) {
        el    <- diag_data[[i]]
        echr  <- sapply(strsplit(names(el), "\\."), "[", 2)
        whc   <- echr %in% chr
        if (!any(whc)) return(NULL)
        data.frame(
            values    = as.numeric(el[whc]),
            chr       = echr[whc],
            iname     = names(el)[whc],
            dist      = cp$pos_lookup[names(el)[whc]],
            iteration = c.iter[i],
            stringsAsFactors = FALSE
        )
    })
    rows <- rows[!sapply(rows, is.null)]
    if (length(rows) == 0L)
        stop("No statistic values found for the requested iterations/chromosomes.")

    df            <- do.call(rbind, rows)
    all.iters     <- paste0("Iteration: ", iter)
    df$iteration  <- factor(df$iteration, levels = all.iters)
    df
}

# =============================================================================
# Shared engine: .add_sig_labels
# Annotate the selected QTL/marker in its corresponding iteration facet.
# =============================================================================
.add_sig_labels <- function(gp, object, stat_slot, iter, chr, cp, sig.col) {
    qtl.keys <- object$QTL$qtl
    if (is.null(qtl.keys) || length(qtl.keys) == 0L) return(gp)

    iter.with.sig <- iter[iter <= length(qtl.keys)]
    diag_data     <- object$QTL$diag[[stat_slot]]

    sig.rows <- lapply(iter.with.sig, function(it) {
        key  <- qtl.keys[it]
        parts <- strsplit(key, "\\.")[[1L]]
        qchr  <- parts[2L]
        if (!qchr %in% chr) return(NULL)
        val  <- diag_data[[it]][key]
        if (is.null(val) || is.na(val) || val == 0) return(NULL)
        dist <- cp$pos_lookup[key]
        if (is.na(dist)) return(NULL)
        label <- sub("Chr\\.", "", key)   # "CHRNAME.IDX" for display
        data.frame(
            dist      = dist,
            values    = as.numeric(val),
            iteration = factor(paste0("Iteration: ", it),
                               levels = levels(gp$data$iteration)),
            label     = label,
            stringsAsFactors = FALSE
        )
    })
    sig.df <- do.call(rbind, sig.rows[!sapply(sig.rows, is.null)])
    if (is.null(sig.df) || nrow(sig.df) == 0L) return(gp)

    gp +
        ggplot2::geom_point(data = sig.df,
                            ggplot2::aes(x = .data$dist, y = .data$values),
                            colour = sig.col, size = 2.5,
                            shape = 18, inherit.aes = FALSE) +
        ggplot2::geom_text(data = sig.df,
                           ggplot2::aes(x = .data$dist, y = .data$values,
                                        label = .data$label),
                           colour = sig.col, size = 2.8,
                           vjust = -0.8, hjust = 0.5,
                           inherit.aes = FALSE)
}

# =============================================================================
# Shared engine: .plot_chr_bars
# Chromosome-level bar chart for selection = "chromosome" analyses.
# =============================================================================
.plot_chr_bars <- function(object, iter, chr, sig.col) {
    c.iter   <- paste0("Iteration: ", iter)
    ochr_raw <- object$QTL$diag$ochr[iter]
    names(ochr_raw) <- c.iter

    rows <- lapply(seq_along(ochr_raw), function(i) {
        el  <- ochr_raw[[i]]
        whc <- names(el) %in% chr
        data.frame(values    = as.numeric(el[whc]),
                   chr       = names(el)[whc],
                   iteration = c.iter[i],
                   stringsAsFactors = FALSE)
    })
    ochrd <- do.call(rbind, rows)
    ochrd$iteration <- factor(ochrd$iteration, levels = c.iter)

    gp <- ggplot2::ggplot(ochrd,
               ggplot2::aes(x = .data$chr, y = .data$values)) +
        ggplot2::facet_wrap(~iteration, ncol = 1) +
        ggplot2::geom_col(fill = "lightblue", colour = "grey50") +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::ylab("Outlier Statistic") +
        ggplot2::xlab("Chromosome") +
        theme_scatter()

    # Annotate selected chromosome
    qtl.keys <- object$QTL$qtl
    if (!is.null(qtl.keys)) {
        iter.sig <- iter[iter <= length(qtl.keys)]
        sig.rows <- lapply(iter.sig, function(it) {
            key  <- qtl.keys[it]
            qchr <- strsplit(key, "\\.")[[1L]][2L]
            if (!qchr %in% chr) return(NULL)
            val  <- object$QTL$diag$ochr[[it]][qchr]
            if (is.null(val) || is.na(val)) return(NULL)
            data.frame(
                chr       = qchr,
                values    = as.numeric(val),
                iteration = factor(paste0("Iteration: ", it), levels = c.iter),
                stringsAsFactors = FALSE
            )
        })
        sig.df <- do.call(rbind, sig.rows[!sapply(sig.rows, is.null)])
        if (!is.null(sig.df) && nrow(sig.df) > 0L)
            gp <- gp +
                ggplot2::geom_col(data = sig.df,
                                  ggplot2::aes(x = .data$chr, y = .data$values),
                                  fill = sig.col, colour = sig.col,
                                  inherit.aes = FALSE)
    }
    gp
}
