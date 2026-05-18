# =============================================================================
# plot.qtlAim.R
# S3 plot method for qtlAim objects.
#
# Plots the forward-selection outlier statistics or scaled BLUPs across the
# genome for each requested iteration, with detected QTL labelled.
#
# Shared internal engine functions used also by plot.gwasAim:
#   .build_cumpos()      -- cumulative cM position lookup from a genObj
#   .build_stat_df()     -- long-format data frame from diag oint/blups vectors
#   .add_sig_labels()    -- annotate selected QTL/markers onto a ggplot
#   .plot_chr_bars()     -- chromosome-level bar chart (plot.chr mode)
#   .build_effects_df()  -- data frame of effect sizes, SEs, % var for lollipop
#   .plot_effects()      -- lollipop plot of QTL/marker effects
#   .build_contrast_df() -- data frame of total genetic values by allele class
#   .plot_contrast()     -- violin/box allele contrast plot
# =============================================================================

#' Plot a Fitted \code{qtlAim} Object
#'
#' @description
#' Produces one of six diagnostic or results plots for a \code{qtlAim} object.
#'
#' @param x A \code{qtlAim} object.
#' @param genObj The \code{"wgCross"} object passed to \code{qtlAim}, produced
#'   by \code{\link{primeCross}}.
#' @param type Character string selecting the plot type:
#'   \describe{
#'     \item{\code{"outlier"} (default)}{Line plot of genome-wide interval/marker
#'       outlier statistics \eqn{\tilde{q}_i^2/\tilde{v}_i}, one facet per
#'       iteration. Detected QTL are annotated.}
#'     \item{\code{"blups"}}{Line plot of genome-wide scaled BLUPs
#'       \eqn{\tilde{q}_i / \sqrt{|\tilde{v}_i|}}, one facet per iteration.}
#'     \item{\code{"chr"}}{Chromosome-level aggregated outlier statistics as a
#'       bar chart (only when \code{qtlAim} was run with
#'       \code{selection = "chromosome"}).}
#'     \item{\code{"effects"}}{Lollipop plot of estimated additive effect sizes
#'       for each detected QTL, with \eqn{\pm} 1 SE error bars, percentage of
#'       phenotypic variance explained annotated on each head, and allele
#'       direction (A / B; or 0/1/2 for GWAS panels) colour-coded.}
#'     \item{\code{"contrast"}}{Allele contrast plot: for each detected QTL one
#'       facet shows the distribution of total line genetic values (additive
#'       genomic BLUP + residual genetic BLUP + QTL fixed contributions) split
#'       by allele class. Requires \code{data} (the phenotypic data frame used
#'       in the analysis). Biparental populations use A / AB / B classes; GWAS
#'       panels use 0 / 1 / 2 dosage classes.}
#'     \item{\code{"heatmap"}}{Genome x iteration heatmap: tile fill encodes
#'       the outlier statistic \eqn{\tilde{q}_i^2/\tilde{v}_i} at every
#'       interval/marker position across all forward-selection iterations.
#'       Positions excluded by the exclusion window are shown in grey.
#'       Detected QTL are marked with a filled diamond at the (position,
#'       iteration) cell where they were selected. Iteration 1 is at the top.
#'       Useful for seeing where genome-wide signal builds, dissipates, or
#'       shifts as each QTL is absorbed into the model.}
#'   }
#' @param data Data frame. The phenotypic data used in the analysis (the
#'   \code{<response>.data} object created by \code{qtlAim} / \code{gwasAim}
#'   in the calling environment). Required only when \code{type = "contrast"}.
#' @param iter Integer vector of iterations to display. Default is all stored
#'   iterations. Not used for \code{"effects"} or \code{"contrast"}.
#' @param chr Optional character vector of chromosome names to display. Not
#'   used for \code{"effects"} or \code{"contrast"}.
#' @param chr.lines Logical. If \code{TRUE}, vertical lines are drawn at
#'   chromosome boundaries (\code{"outlier"} and \code{"blups"} only).
#'   Default \code{FALSE}.
#' @param sig.col Colour used to highlight and label the selected QTL in each
#'   iteration facet (\code{"outlier"} / \code{"blups"} / \code{"chr"} types).
#'   Default \code{"firebrick"}.
#' @param pt.col Character vector of length equal to the number of displayed
#'   chromosomes, or length 2 for alternating colours. Colours the line/rug
#'   traces by chromosome. Default uses \code{hcl.colors}.
##' @param pt.cex Numeric line/point width. Default \code{0.6}.
#' @param ncol Integer. Number of columns in the \code{facet_wrap} layout for
#'   \code{type = "contrast"}. Default \code{1} (facets stacked vertically,
#'   consistent with the \code{"outlier"} and \code{"blups"} layout). Set to
#'   \code{NULL} to let ggplot2 choose automatically, or any positive integer
#'   for a multi-column grid. Not used for other plot types.
#' @param cap Numeric scalar. Upper limit for the colour scale in
#'   \code{type = "heatmap"}. Outlier statistics above this value are clamped
#'   to \code{cap} and rendered in the maximum colour, stretching the gradient
#'   across the lower range where most signal variation occurs. Default
#'   \code{NULL} (no capping; full data range used). Not used for other plot
#'   types.
#' @param \dots Currently unused.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{summary.qtlAim}},
#'   \code{\link{aimTrace}}, \code{\link{linkMap.qtlAim}}
#'
#' @examples
#' \dontrun{
#' plot(qtl.fit, genObj = genoRxK)
#' plot(qtl.fit, genObj = genoRxK, type = "effects")
#' plot(qtl.fit, genObj = genoRxK, type = "heatmap")
#' plot(qtl.fit, genObj = genoRxK, type = "contrast", data = yld.data)
#' }
#'
#' @export
plot.qtlAim <- function(x, genObj,
                         type      = c("outlier", "blups", "chr",
                                       "effects", "contrast", "heatmap"),
                         data      = NULL,
                         ncol      = 1L,
                         iter      = NULL,
                         chr       = NULL,
                         chr.lines = FALSE,
                         sig.col   = "firebrick",
                         pt.col    = NULL,
                         pt.cex    = 0.6,
                         cap       = NULL,
                         ...) {

    type <- match.arg(type)

    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("genObj must be of class \"wgCross\" produced by primeCross().")

    if (is.null(x$QTL$effects))
        stop("No significant QTL found in object.")

    # -------------------------------------------------------------------------
    # effects: lollipop plot -- no iter/chr needed
    # -------------------------------------------------------------------------
    if (type == "effects") {
        if (!is.null(x$QTL$Trait)) {
            edf <- .build_mv_effects_df(x, genObj)
            return(.plot_mv_effects(edf, x))
        }
        edf <- .build_effects_df(x, genObj)
        return(.plot_effects(edf))
    }

    # -------------------------------------------------------------------------
    # contrast: allele contrast plot -- requires data
    # -------------------------------------------------------------------------
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
    # heatmap: genome x iteration outlier heatmap
    # -------------------------------------------------------------------------
    if (type == "heatmap")
        return(.plot_heatmap(x, genObj, iter, chr, sig.col, cap))

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

    # Multivariate blups: one line per trial with flag annotations
    if (type == "blups" && !is.null(x$QTL$Trait)) {
        mv_df        <- .build_mv_blups_df(x, iter, chr, cp)
        trial.levels <- x$QTL$trait.levels
        trial.cols   <- grDevices::hcl.colors(length(trial.levels), "Dark 2")
        names(trial.cols) <- trial.levels
        mv_df$trial  <- factor(mv_df$trial, levels = trial.levels)

        gp <- ggplot2::ggplot(mv_df,
                   ggplot2::aes(x      = .data$dist,
                                y      = .data$values,
                                colour = .data$trial,
                                group  = interaction(.data$chr, .data$trial))) +
            ggplot2::facet_wrap(~iteration, ncol = 1, scales = "free_y") +
            ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                                linewidth = 0.3, linetype = "dashed") +
            ggplot2::geom_line(linewidth = pt.cex) +
            ggplot2::scale_colour_manual(values = trial.cols,
                                         name   = x$QTL$Trait) +
            ggplot2::scale_x_continuous(
                breaks = cp$chr.mid, labels = names(cp$chr.mid)) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.05, 0.25))) +
            ggplot2::ylab("Scaled BLUPs") +
            ggplot2::xlab("") +
            ggplot2::coord_cartesian(clip = "off") +
            theme_scatter()
        return(.add_mv_blups_flags(gp, x, iter, chr, cp, sig.col))
    }

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
            expand = ggplot2::expansion(mult = c(0.02, 0.25))) +
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
#   $pos_lookup  -- named vector: internal key "Chr.CHR.IDX" -> cumulative cM
#   $chr.mid     -- named vector: chr name -> mid-point cumulative cM
#   $chr.end     -- named vector: chr name -> end cumulative cM (for chr.lines)
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
        ggrepel::geom_text_repel(
            data        = sig.df,
            ggplot2::aes(x = .data$dist, y = .data$values,
                         label = .data$label),
            colour      = sig.col,
            size        = 2.8,
            nudge_y     = diff(range(sig.df$values, na.rm = TRUE)) * 0.08 + 0.5,
            direction   = "x",
            segment.size    = 0.3,
            segment.colour  = sig.col,
            segment.alpha   = 0.6,
            box.padding     = 0.2,
            point.padding   = 0.3,
            force           = 1,
            max.overlaps    = Inf,
            inherit.aes     = FALSE)
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

# =============================================================================
# Shared engine: .match_perc_var
#
# Extracts the "Perc.Var" column from a summary() table and matches each
# element back to the corresponding entry in 'effects' (which is in ASReml
# detection order, possibly reversed by rev() in engine_effects.R).
#
# Matching key: chromosome name + rounded effect size (4 dp).  Within a
# chromosome, two QTL with identical rounded effects are astronomically
# unlikely, but if they occur the first match is used and a warning is issued.
# =============================================================================
.match_perc_var <- function(summ, object, genObj, effects, gen.type) {
    qtlm     <- getQTL(object, genObj)
    chr_summ <- as.character(summ[["Chromosome"]])

    # Match by cM position, not effect size.  Two detected QTL on the same
    # chromosome are always >= exclusion.window cM apart, so position is an
    # unambiguous key.  Effect sizes can be coincidentally close on crowded
    # chromosomes, making the old size-based match unreliable.
    #
    # Column index for cM position in the summary data.frame:
    #   interval type : col 5 = inferred midpoint cM
    #                   (summary layout: Chr | LMark | dist | IMark | dist | RMark | dist | Size | Prob | Perc.Var)
    #   marker/gwas   : col 3 = marker cM
    #                   (summary layout: Chr | Marker | dist | Size | Prob | Perc.Var)
    pos_col  <- if (gen.type == "interval") 5L else 3L
    pos_summ <- as.numeric(summ[, pos_col])

    # Corresponding position from getQTL():
    #   interval: col 6 = round(imark, 2) = inferred midpoint cM
    #   marker  : col 4 = round(lhmark, 2) = marker cM
    pos_qtlm_col <- if (gen.type == "interval") 6L else 4L

    vapply(seq_along(effects), function(i) {
        chr_i <- as.character(qtlm[i, 1L])
        pos_i <- as.numeric(qtlm[i, pos_qtlm_col])
        j     <- which(chr_summ == chr_i & abs(pos_summ - pos_i) < 0.1)
        if (length(j) == 0L) {
            warning("Could not match QTL ", names(effects)[i],
                    " in summary table -- Perc.Var set to NA.")
            return(NA_real_)
        }
        summ[["Perc.Var"]][j[1L]]
    }, numeric(1L))
}

# =============================================================================
# Shared engine: .build_effects_df
#
# Builds a data frame with one row per detected QTL/marker containing:
#   qtl_key  -- internal "Chr.CHR.IDX" key
#   label    -- display label "CHR . POS cM"
#   effect   -- estimated additive effect (beta)
#   se       -- standard error  sqrt(veffect x sigma2)
#   lo / hi  -- effect +/- 1 SE
#   perc_var -- % phenotypic variance explained
#   direction -- "A" (effect > 0) or "B" (effect < 0), or "0","1","2" for GWAS
#   col      -- fill colour for lollipop head
#   x_order  -- genome-order integer for x-axis sorting
#
# Works for both qtlAim (interval/marker type) and gwasAim (marker type).
# =============================================================================
.build_effects_df <- function(object, genObj) {

    effects  <- object$QTL$effects
    veffects <- object$QTL$veffects
    sigma2   <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4L)
        sigma2 <- 1

    gen.type <- object$QTL$type

    # % variance explained -- taken directly from summary() so the values
    # shown on the plot are guaranteed to match the summary table output.
    summ      <- summary(object, genObj, LOD = FALSE)
    perc_var  <- .match_perc_var(summ, object, genObj, effects, gen.type)

    # SE and confidence limits
    se <- sqrt(veffects * sigma2)

    # Position labels from getQTL
    qtlm <- getQTL(object, genObj)

    rows <- lapply(seq_along(effects), function(i) {
        eff  <- effects[i]
        nm   <- names(effects)[i]    # "X.CHR.IDX"
        chr  <- qtlm[i, 1L]
        pos_cM <- if (gen.type == "interval") as.numeric(qtlm[i, 6L])
                  else                        as.numeric(qtlm[i, 4L])
        lbl <- paste0(chr, " \u00b7 ", round(pos_cM, 1), " cM")

        dir <- if (eff >= 0) "A" else "B"
        col <- if (eff >= 0) "firebrick" else "steelblue"

        chr_order <- match(chr, names(genObj$geno))

        # % var label sits just outside the SE bar endpoint, not the head
        label_y     <- if (eff >= 0) eff + se[i] else eff - se[i]
        label_vjust <- if (eff >= 0) -0.7 else 1.7

        data.frame(
            qtl_key     = nm,
            label       = lbl,
            chr         = chr,
            pos_cM      = pos_cM,
            effect      = eff,
            se          = se[i],
            lo          = eff - se[i],
            hi          = eff + se[i],
            perc_var    = perc_var[i],
            label_y     = label_y,
            label_vjust = label_vjust,
            direction   = dir,
            col         = col,
            x_order     = chr_order * 1e6 + pos_cM,
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)
    df[order(df$x_order), ]
}

# =============================================================================
# Shared engine: .plot_effects
#
# Lollipop plot of additive QTL/marker effect sizes.
#   -- stem from 0 to the effect estimate
#   -- dot sized to emphasise position, coloured by direction
#   -- +/-1 SE error bars
#   -- % variance explained annotated above/below each head
#   -- allele axis annotation:  "<- B allele favoured" / "A allele favoured ->"
# =============================================================================
.plot_effects <- function(edf) {

    # Ordered factor for x-axis (genome order, left to right)
    edf$label    <- factor(edf$label, levels = edf$label)
    edf$pv_label <- paste0(edf$perc_var, "%")

    # Axis subtitle with allele encoding reminder
    x_sub <- "QTL position  \u2014  firebrick = A allele favoured (+1),  steelblue = B allele favoured (\u22121)"

    gp <- ggplot2::ggplot(edf,
               ggplot2::aes(x = .data$label, y = .data$effect)) +

        # Zero reference line
        ggplot2::geom_hline(yintercept = 0, colour = "grey50",
                            linewidth = 0.5, linetype = "dashed") +

        # Stems (from 0 to estimate)
        ggplot2::geom_segment(
            ggplot2::aes(x    = .data$label, xend = .data$label,
                         y    = 0,           yend = .data$effect,
                         colour = .data$col),
            linewidth = 0.8, show.legend = FALSE) +

        # +/-1 SE error bars
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$lo, ymax = .data$hi,
                         colour = .data$col),
            width = 0.25, linewidth = 0.7, show.legend = FALSE) +

        # Lollipop heads
        ggplot2::geom_point(
            ggplot2::aes(colour = .data$col),
            size = 3.5, show.legend = FALSE) +

        # % variance explained -- anchored at the SE bar tip, not the head
        ggplot2::geom_text(
            ggplot2::aes(y     = .data$label_y,
                         label = .data$pv_label,
                         vjust = .data$label_vjust),
            size = 3, colour = "grey30") +

        ggplot2::scale_colour_identity() +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.18, 0.18))) +
        ggplot2::ylab("Additive Effect") +
        ggplot2::xlab(x_sub) +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter() +
        ggplot2::theme(
            # angle = 45, hjust = 1, vjust = 1 aligns the top-right corner
            # of the label to the tick mark -- the standard "flush" appearance
            axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1,
                                                  vjust = 1, size = 8),
            axis.title.x = ggplot2::element_text(size = 7.5,
                                                  colour = "grey40",
                                                  margin = ggplot2::margin(t = 8))
        )
    gp
}

# =============================================================================
# Shared engine: .build_contrast_df
#
# Builds the data frame needed for the allele contrast plot.
#
# Total genetic value for each line:
#   total_i  =  g_combined_i  +  qtl_contrib_i
#
# where:
#   g_combined_i  =  vm(gterm, covObj) BLUP  +  gterm/Gsave BLUP
#                    obtained from a single predict() call with both terms
#                    in the 'only' argument.
#
#   qtl_contrib_i  =  Sigma_k ( beta_k x x_ik )
#
# Allele scoring strategy (differs by object class):
#
#   qtlAim (biparental) -- uses the same data matrix qtlAim.asreml() used:
#     interval.data[, idx_k] when gen.type == "interval" (idx_k indexes the
#     inferred.map / interval columns directly), or imputed.data[, idx_k]
#     when gen.type == "marker" (idx_k indexes the marker map columns).
#     Both are post-imputation, so values cluster near +/-1 with no spurious
#     zeros from missing raw genotypes.  No AB class: biparental (BC/DH/RI)
#     populations have no true heterozygotes.
#     x-axis in the plot is the continuous score; labels pinned at +/-1.
#
#   gwasAim (diversity panel) -- uses the X.<chr>.<idx> column in data
#     (integer dosage 0/1/2).  Three discrete classes.
#
# Returns a list (one data frame per QTL) with attribute "is_gwas" so
# .plot_contrast() can choose the correct geom strategy.
# =============================================================================
.build_contrast_df <- function(object, genObj, data) {

    effects  <- object$QTL$effects
    veffects <- object$QTL$veffects
    sigma2   <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4L)
        sigma2 <- 1

    gen.type <- object$QTL$type
    is_gwas  <- inherits(object, "gwasAim")
    gterm    <- object$QTL$diag$genetic.term

    # ----------------------------------------------------------------
    # 1.  Predict combined additive + residual genetic BLUPs per line
    # ----------------------------------------------------------------
    rterms    <- attr(stats::terms(object$call$random), "term.labels")
    vm_label  <- grep("vm\\(|mbf\\(", rterms, value = TRUE)
    res_label <- if ("Gsave" %in% colnames(object$mf)) "Gsave" else gterm

    only_terms <- unique(c(vm_label, res_label))
    pv <- predict(object,
                  classify = res_label,
                  only     = only_terms,
                  data     = data)
    g_df <- pv$pvals[, c(res_label, "predicted.value")]
    colnames(g_df) <- c("line", "g_combined")
    g_df$line <- as.character(g_df$line)

    # ----------------------------------------------------------------
    # 2.  Compute QTL contributions  Sigma_k beta_k x x_ik  per line
    # ----------------------------------------------------------------
    qtl_cols         <- gsub("Chr\\.", "X.", names(effects))
    qtl_cols_present <- qtl_cols[qtl_cols %in% names(data)]
    if (length(qtl_cols_present) < length(qtl_cols))
        warning("Some QTL genotype columns not found in data -- ",
                "QTL contributions may be incomplete.")

    id_col      <- object$QTL$diag$genetic.term
    line_ids    <- as.character(data[[id_col]])
    qtl_mat     <- as.matrix(data[, qtl_cols_present, drop = FALSE])
    eff_vec     <- effects[match(qtl_cols_present,
                                 gsub("Chr\\.", "X.", names(effects)))]
    qtl_contrib <- as.vector(qtl_mat %*% eff_vec)

    contrib_df <- stats::aggregate(
        data.frame(qtl_contrib = qtl_contrib),
        by  = list(line = line_ids),
        FUN = mean)

    # ----------------------------------------------------------------
    # 3.  Merge and compute total genetic value
    # ----------------------------------------------------------------
    merged <- merge(g_df, contrib_df, by = "line", all.x = TRUE)
    merged$total_genetic <- merged$g_combined +
        ifelse(is.na(merged$qtl_contrib), 0, merged$qtl_contrib)

    # ----------------------------------------------------------------
    # 4.  % var explained and p-values
    # ----------------------------------------------------------------
    qtlm <- getQTL(object, genObj)

    # Use summary() for perc_var -- guaranteed to match the summary table.
    summ     <- summary(object, genObj, LOD = FALSE)
    perc_var <- .match_perc_var(summ, object, genObj, effects, gen.type)

    zrat   <- effects / sqrt(veffects * sigma2)
    pvalue <- if (object$QTL$method == "random")
        round((1 - pchisq(zrat^2, 1)) / 2, 4)
    else
        round(2 * (1 - pnorm(abs(zrat))), 4)
    pvalue[as.numeric(pvalue) < 0.0001] <- "<0.0001"

    # ----------------------------------------------------------------
    # 5.  Per-QTL allele scoring and facet data assembly
    # ----------------------------------------------------------------
    result_list <- lapply(seq_along(effects), function(i) {

        chr_k <- qtlm[i, 1L]
        idx_k <- as.integer(qtlm[i, 2L])

        if (is_gwas) {
            # ----- GWAS: use X.<chr>.<idx> column from data ------------------
            # The column holds the +/-1-centred encoding (0/1/2 -> -1/0/+1)
            # that gwasAim.asreml() uses internally.  Map back to 0/1/2 labels
            # for display: -1 -> "0 (ref)", 0 -> "1 (het)", +1 -> "2 (alt)".
            col_x <- qtl_cols[i]
            if (!col_x %in% names(data)) {
                warning("Column ", col_x, " not found in data -- skipping QTL ", i)
                return(NULL)
            }
            raw_score  <- data[[col_x]]
            line_ids_i <- as.character(data[[id_col]])
            score_df   <- stats::aggregate(raw_score ~ line_ids_i, FUN = mean)
            colnames(score_df) <- c("line", "score")

            score_df$allele_class <- factor(
                as.character(round(score_df$score)),
                levels = c("-1", "0", "1"),
                labels = c("0", "1", "2"))
            allele_cols <- c("0" = "steelblue", "1" = "grey60", "2" = "firebrick")

        } else {
            # ----- Biparental: use the same data matrix the model used -------
            # gen.type == "interval": model used interval.data; idx_k indexes
            #   the inferred.map / interval columns directly.
            # gen.type == "marker":  model used imputed.data; idx_k indexes
            #   the marker map / imputed columns directly.
            # Both are post-imputation so values cluster near +/-1 -- no spurious
            # zeros from missing raw genotypes.  No AB class for biparental
            # (BC/DH/RI) populations.
            geno_mat <- if (gen.type == "interval")
                genObj$geno[[chr_k]]$interval.data
            else
                genObj$geno[[chr_k]]$imputed.data

            imp_vec  <- geno_mat[, idx_k]
            score_df <- data.frame(
                line  = rownames(geno_mat),
                score = as.numeric(imp_vec),
                stringsAsFactors = FALSE
            )
            score_df$allele_class <- factor(
                ifelse(score_df$score < 0, "B (\u22121)", "A (+1)"),
                levels = c("B (\u22121)", "A (+1)"))
            allele_cols <- c("B (\u22121)" = "steelblue",
                             "A (+1)"     = "firebrick")
        }

        # Merge total genetic values with allele scores
        out <- merge(merged[, c("line", "total_genetic")],
                     score_df[, c("line", "score", "allele_class")],
                     by = "line")

        # Facet label and effect annotation
        pos_cM  <- if (gen.type == "interval") as.numeric(qtlm[i, 6L])
                   else                        as.numeric(qtlm[i, 4L])
        se_i    <- sqrt(veffects[i] * sigma2)
        eff_txt <- sprintf("Effect: %+.3f \u00b1 %.3f  (%s%%  var)   p = %s",
                           effects[i], se_i, perc_var[i], pvalue[i])

        out$facet_label     <- paste0("Chr ", chr_k, "  \u00b7  ",
                                      round(pos_cM, 1), " cM")
        out$effect_txt      <- eff_txt
        out$allele_cols_map <- list(allele_cols)
        out
    })

    result_list <- result_list[!sapply(result_list, is.null)]
    attr(result_list, "is_gwas") <- is_gwas
    result_list
}

# =============================================================================
# Shared engine: .plot_contrast
#
# Allele contrast plot.  One facet per QTL.  Two visual strategies:
#
#   Biparental (qtlAim) -- scatter / jitter plot.
#     x = imputed genotype score (continuous, post-imputation ~= +/-1).
#     Points jittered in x only (width = 0.03) to reveal density.
#     Colour: steelblue (B, score < 0) / firebrick (A, score > 0).
#     x-axis pinned to [-1.2, 1.2] with breaks and labels at -1 / +1.
#     Group means +/- SE drawn at x = -1 and x = +1, connected by a line.
#
#   GWAS (gwasAim) -- violin / box / jitter plot.
#     x = dosage class (discrete: 0 / 1 / 2).
#     Colour: steelblue (0) / grey60 (1) / firebrick (2).
#
# Both paths:
#   y = total genetic value (vm BLUP + gterm BLUP + QTL contributions).
#   Effect +/- SE, % var explained, p-value annotated inside each facet strip.
# =============================================================================
.plot_contrast <- function(cdf_list, ncol = 1L) {

    is_gwas <- isTRUE(attr(cdf_list, "is_gwas"))

    # Combine all per-QTL data frames
    all_df <- do.call(rbind, lapply(seq_along(cdf_list), function(i) {
        df <- cdf_list[[i]]
        df$qtl_idx <- i
        df
    }))

    facet_levels       <- unique(all_df$facet_label)
    all_df$facet_label <- factor(all_df$facet_label, levels = facet_levels)

    allele_cols <- cdf_list[[1L]]$allele_cols_map[[1L]]
    all_df$allele_class <- factor(all_df$allele_class, levels = names(allele_cols))

    sub_df <- unique(all_df[, c("facet_label", "effect_txt")])

    # ------------------------------------------------------------------
    # Shared mean_df: per-facet group means (no SE needed)
    # ------------------------------------------------------------------
    mean_df <- do.call(rbind, lapply(levels(all_df$facet_label), function(fl) {
        sub <- all_df[all_df$facet_label == fl, ]
        agg <- stats::aggregate(
            total_genetic ~ allele_class + facet_label,
            data = sub,
            FUN  = mean, na.rm = TRUE)
        agg$facet_label <- factor(agg$facet_label, levels = facet_levels)
        agg
    }))

    # ------------------------------------------------------------------
    # Biparental path: scatter + jitter on continuous imputed score
    # ------------------------------------------------------------------
    if (!is_gwas) {

        # Mean positions pinned to the canonical +/-1 encoding values
        mean_df$x_pos <- ifelse(mean_df$allele_class == "B (\u22121)", -1, 1)

        gp <- ggplot2::ggplot(all_df,
                   ggplot2::aes(x      = .data$score,
                                y      = .data$total_genetic,
                                colour = .data$allele_class)) +
            ggplot2::facet_wrap(~ facet_label, ncol = ncol,
                                scales = "free_y") +

            # Jittered points -- horizontal jitter only; y is meaningful
            ggplot2::geom_jitter(width = 0.05, height = 0,
                                 size = 1.2, alpha = 0.55,
                                 show.legend = FALSE) +

            # Mean diamond at fixed x = +/-1
            ggplot2::geom_point(
                data   = mean_df,
                ggplot2::aes(x = .data$x_pos, y = .data$total_genetic),
                colour = "grey15", size = 2.8,
                shape  = 18, inherit.aes = FALSE) +

            # Line connecting the two group means
            ggplot2::geom_line(
                data   = mean_df,
                ggplot2::aes(x     = .data$x_pos,
                             y     = .data$total_genetic,
                             group = .data$facet_label),
                colour = "grey15", linewidth = 0.7,
                inherit.aes = FALSE) +

            ggplot2::scale_colour_manual(values = allele_cols, guide = "none") +
            ggplot2::scale_x_continuous(
                limits = c(-1.2, 1.2),
                breaks = c(-1, 1),
                labels = c("B (\u22121)", "A (+1)")) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
            ggplot2::xlab("Allele  (imputed genotype score)") +
            ggplot2::ylab("Total Genetic Value")

    } else {
    # ------------------------------------------------------------------
    # GWAS path: jitter scatter on discrete dosage classes (0 / 1 / 2)
    # ------------------------------------------------------------------
        gp <- ggplot2::ggplot(all_df,
                   ggplot2::aes(x      = .data$allele_class,
                                y      = .data$total_genetic,
                                colour = .data$allele_class)) +
            ggplot2::facet_wrap(~ facet_label, ncol = ncol,
                                scales = "free_y") +

            # Jittered points around each discrete dosage class --
            # small width matches the tight vertical clouds of the QTL plot
            ggplot2::geom_jitter(width = 0.05, height = 0,
                                 size = 1.2, alpha = 0.55,
                                 show.legend = FALSE) +

            # Mean diamond per dosage class
            ggplot2::geom_point(
                data   = mean_df,
                ggplot2::aes(x = .data$allele_class,
                             y = .data$total_genetic),
                colour = "grey15", size = 2.8,
                shape  = 18, inherit.aes = FALSE) +

            # Line connecting the three group means (additive trend)
            ggplot2::geom_line(
                data   = mean_df,
                ggplot2::aes(x     = .data$allele_class,
                             y     = .data$total_genetic,
                             group = .data$facet_label),
                colour = "grey15", linewidth = 0.7,
                inherit.aes = FALSE) +

            ggplot2::scale_colour_manual(values = allele_cols, guide = "none") +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
            ggplot2::xlab("Dosage class  (0 = ref hom,  1 = het,  2 = alt hom)") +
            ggplot2::ylab("Total Genetic Value")
    }

    # Effect annotation inside each facet -- shared by both paths
    gp <- gp +
        ggplot2::geom_text(
            data        = sub_df,
            ggplot2::aes(x = -Inf, y = Inf, label = .data$effect_txt),
            inherit.aes = FALSE,
            hjust = -0.04, vjust = 1.6,
            size = 2.8, colour = "grey30") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter()
    gp
}

# =============================================================================
# Shared engine: .build_heatmap_df
#
# Builds the data frame for the genome x iteration heatmap.
#
# Returns a list with two elements:
#   $heat_df  -- long-format data frame, one row per (position x iteration):
#       dist      cumulative cM position
#       iter_num  integer iteration index
#       fill_val  outlier statistic (zero -> NA so excluded regions show grey)
#       tile_w    per-position tile width (Voronoi half-distance to neighbours
#                 within the chromosome; clipped at chromosome start/end)
#       chr       chromosome name
#
#   $sig_df   -- data frame of selected QTL positions (dist, iter_num), or
#               NULL if no QTL were selected
#
# Zero/negative -> NA strategy: positions excluded via the exclusion window are
# set to 0 in the oint vector by .qtlSelect().  Tiny floating-point negatives
# can also arise from numerical precision.  Both are converted to NA so they
# render as neutral grey in the colour scale, making active signal regions
# immediately visible.  The fill is the raw outlier statistic on a linear scale.
# =============================================================================
.build_heatmap_df <- function(object, iter, chr, cp) {

    # ---- per-position Voronoi tile widths ------------------------------------
    # Compute half-distance to neighbours within each chromosome block so
    # tiles fill the space without gaps or overlaps.
    all_keys <- names(cp$pos_lookup)
    tile_w   <- numeric(length(all_keys))
    names(tile_w) <- all_keys

    for (ch in chr) {
        ch_keys <- grep(paste0("^Chr\\.", ch, "\\."), all_keys, value = TRUE)
        if (length(ch_keys) == 0L) next
        ch_pos  <- sort(cp$pos_lookup[ch_keys])     # sorted cumulative cM
        n       <- length(ch_pos)
        if (n == 1L) {
            tile_w[names(ch_pos)] <- 2              # arbitrary 2 cM for singleton
            next
        }
        left_bounds  <- c(ch_pos[1L],
                          (ch_pos[-n]  + ch_pos[-1L]) / 2)
        right_bounds <- c((ch_pos[-n]  + ch_pos[-1L]) / 2,
                          ch_pos[n])
        w <- right_bounds - left_bounds
        names(w) <- names(ch_pos)
        tile_w[names(w)] <- w
    }

    # ---- build long-format heat data frame -----------------------------------
    rows <- lapply(iter, function(it) {
        el    <- object$QTL$diag$oint[[it]]
        echr  <- sapply(strsplit(names(el), "\\."), "[", 2L)
        whc   <- echr %in% chr
        if (!any(whc)) return(NULL)

        keys  <- names(el)[whc]
        vals  <- as.numeric(el[whc])
        vals[is.finite(vals) & vals <= 0] <- NA  # excluded/fp-zero -> NA (grey)

        data.frame(
            dist     = cp$pos_lookup[keys],
            iter_num = it,
            fill_val = vals,
            tile_w   = tile_w[keys],
            chr      = echr[whc],
            iname    = keys,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    })
    rows    <- rows[!sapply(rows, is.null)]
    heat_df <- do.call(rbind, rows)

    # ---- selected QTL positions ----------------------------------------------
    qtl_keys <- object$QTL$qtl
    if (!is.null(qtl_keys) && length(qtl_keys) > 0L) {
        iter_sig <- iter[iter <= length(qtl_keys)]
        sig_rows <- lapply(iter_sig, function(it) {
            key   <- qtl_keys[it]
            parts <- strsplit(key, "\\.")[[1L]]
            qchr  <- parts[2L]
            if (!qchr %in% chr) return(NULL)
            dist  <- cp$pos_lookup[key]
            if (is.na(dist)) return(NULL)
            data.frame(dist = dist, iter_num = it,
                       stringsAsFactors = FALSE)
        })
        sig_df <- do.call(rbind, sig_rows[!sapply(sig_rows, is.null)])
        if (nrow(sig_df) == 0L) sig_df <- NULL
    } else {
        sig_df <- NULL
    }

    list(heat_df = heat_df, sig_df = sig_df)
}

# =============================================================================
# Shared engine: .plot_heatmap
#
# Genome x iteration outlier-statistic heatmap.
#
# Visual conventions:
#   -- x-axis   : cumulative cM, chromosome names at midpoints
#   -- y-axis   : iteration number, reversed (iteration 1 at top) so the
#                sequence reads top-to-bottom, matching the faceted line plots
#   -- fill     : outlier statistic, sqrt-transformed for skewness.
#                NA (excluded) positions shown in grey ("grey92").
#                Colour ramp: light cream -> teal -> dark navy (YlGnBu-style,
#                using base-R hcl.colors so no extra dependencies).
#   -- borders  : thin white vertical lines at chromosome boundaries
#   -- QTL mark : open diamond (shape 23) with sig.col fill, white border,
#                at the (position, iteration) cell where each QTL was selected
#
# Works for both qtlAim (gen.type = "interval" or "marker") and gwasAim
# (always "marker") because .build_cumpos() abstracts the difference.
# =============================================================================
.plot_heatmap <- function(object, genObj, iter, chr, sig.col, cap = NULL) {

    gen.type <- object$QTL$type
    cp       <- .build_cumpos(genObj, gen.type, chr)
    hd       <- .build_heatmap_df(object, iter, chr, cp)
    heat_df  <- hd$heat_df
    sig_df   <- hd$sig_df

    # Apply cap: clamp fill_val so the gradient spans 0 -> cap and anything
    # above cap maps to the maximum colour.  This stretches colour variation
    # into the lower range without transforming the data in any other way.
    if (!is.null(cap) && is.finite(cap) && cap > 0)
        heat_df$fill_val <- pmin(heat_df$fill_val, cap)

    # YlGnBu-style ramp via base-R hcl.colors (no viridis dependency)
    heat_cols <- grDevices::hcl.colors(256, palette = "YlGnBu", rev = TRUE)

    gp <- ggplot2::ggplot(heat_df,
               ggplot2::aes(x      = .data$dist,
                            y      = .data$iter_num,
                            width  = .data$tile_w,
                            fill   = .data$fill_val)) +
        ggplot2::geom_tile(height = 0.85, colour = NA) +
        ggplot2::scale_fill_gradientn(
            colours  = heat_cols,
            na.value = "grey92",
            name     = "Outlier",
            guide    = ggplot2::guide_colourbar(
                barwidth  = 0.8,
                barheight = 8,
                ticks     = TRUE)) +
        ggplot2::scale_x_continuous(
            breaks = cp$chr.mid,
            labels = names(cp$chr.mid),
            expand = ggplot2::expansion(mult = 0.01)) +
        ggplot2::scale_y_reverse(
            breaks = iter,
            labels = paste0("Iter. ", iter),
            expand = ggplot2::expansion(add = 0.5)) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position  = "right")

    # White chromosome boundary lines
    if (length(chr) > 1L)
        gp <- gp + ggplot2::geom_vline(
            xintercept = cp$chr.end[-length(cp$chr.end)],
            colour     = "white",
            linewidth  = 0.5)

    # Selected QTL diamonds
    if (!is.null(sig_df) && nrow(sig_df) > 0L)
        gp <- gp + ggplot2::geom_point(
            data        = sig_df,
            ggplot2::aes(x = .data$dist, y = .data$iter_num),
            shape       = 23,
            fill        = sig.col,
            colour      = "white",
            size        = 3,
            stroke      = 0.8,
            inherit.aes = FALSE)

    gp
}

# =============================================================================
# Shared engine: .build_mv_blups_df
#
# Builds a long-format data frame for the multivariate blups plot.
# The blups slot for ntrait > 1 is a list of matrices (one per iteration),
# each of dimension (n_all_markers x ntrait) with rownames = marker keys
# and colnames = trait levels.
#
# Returns a data frame with columns:
#   values    -- scaled BLUP z-score
#   chr       -- chromosome name
#   iname     -- marker key "Chr.CHR.IDX"
#   dist      -- cumulative cM position
#   iteration -- factor "Iteration: i"
#   trial     -- trial/environment name (column from the blups matrix)
# =============================================================================
.build_mv_blups_df <- function(object, iter, chr, cp) {
    # Guard: blups slot must be matrices (ntrait > 1 format from current engine).
    # Old in-memory objects (fitted before this engine version) store named
    # vectors instead.  Re-running the analysis produces the correct format.
    first <- object$QTL$diag$blups[[iter[1L]]]
    if (!is.matrix(first))
        stop("The blups slot contains univariate-format vectors, not per-trial ",
             "matrices.\nPlease re-run qtlAim() / gwasAim() to regenerate the ",
             "object with the updated engine.")

    diag_data <- object$QTL$diag$blups[iter]
    c.iter    <- paste0("Iteration: ", iter)
    names(diag_data) <- c.iter

    rows <- lapply(seq_along(diag_data), function(i) {
        mat  <- diag_data[[i]]           # nall_markers x ntrait matrix
        echr <- sapply(strsplit(rownames(mat), "\\."), "[", 2L)
        whc  <- echr %in% chr  # keep all markers including exclusion-window zeros
        if (!any(whc)) return(NULL)

        mat_sub <- mat[whc, , drop = FALSE]
        keys    <- rownames(mat_sub)

        do.call(rbind, lapply(seq_len(ncol(mat_sub)), function(j) {
            data.frame(
                values    = mat_sub[, j],
                chr       = echr[whc],
                iname     = keys,
                dist      = cp$pos_lookup[keys],
                iteration = c.iter[i],
                trial     = colnames(mat_sub)[j],
                stringsAsFactors = FALSE
            )
        }))
    })
    rows <- rows[!sapply(rows, is.null)]
    if (length(rows) == 0L)
        stop("No multivariate BLUP values found for the requested ",
             "iterations/chromosomes.")

    df           <- do.call(rbind, rows)
    df$iteration <- factor(df$iteration, levels = paste0("Iteration: ", iter))
    df
}

# =============================================================================
# Shared engine: .add_mv_blups_flags
#
# Annotates selected QTL on the multivariate blups plot using a "flag" style:
#   -- a faint dashed vertical line at the selected QTL's cumulative cM
#      position, drawn only in its corresponding iteration facet
#   -- a text label fixed at y = Inf (top of each panel) showing:
#        line 1: "CHR.IDX"
#        line 2: "[MAIN]" or "[INT]" when x$QTL$Trait is non-NULL
#
# This approach is independent of the spread of trial BLUP values at the
# selected position, so it works equally well for MAIN QTL (all trials peak
# together) and INTERACTION/crossover QTL (trials diverge in sign).
# =============================================================================
.add_mv_blups_flags <- function(gp, object, iter, chr, cp, sig.col) {
    qtl.keys <- object$QTL$qtl
    if (is.null(qtl.keys) || length(qtl.keys) == 0L) return(gp)

    is.mv          <- !is.null(object$QTL$Trait)
    iter.with.sig  <- iter[iter <= length(qtl.keys)]
    iter.levels    <- levels(gp$data$iteration)

    sig.rows <- lapply(iter.with.sig, function(it) {
        key   <- qtl.keys[it]
        parts <- strsplit(key, "\\.")[[1L]]
        qchr  <- parts[2L]
        if (!qchr %in% chr) return(NULL)
        dist <- cp$pos_lookup[key]
        if (is.na(dist)) return(NULL)

        base_label <- sub("Chr\\.", "", key)   # "CHRNAME.IDX"
        if (is.mv && !is.null(object$QTL$is.interaction)) {
            type_tag <- if (object$QTL$is.interaction[it]) "[INT]" else "[MAIN]"
            label    <- paste0(base_label, "\n", type_tag)
        } else {
            label <- base_label
        }

        data.frame(
            dist      = dist,
            iteration = factor(paste0("Iteration: ", it), levels = iter.levels),
            label     = label,
            stringsAsFactors = FALSE
        )
    })
    sig.df <- do.call(rbind, sig.rows[!sapply(sig.rows, is.null)])
    if (is.null(sig.df) || nrow(sig.df) == 0L) return(gp)

    gp +
        ggplot2::geom_vline(
            data        = sig.df,
            ggplot2::aes(xintercept = .data$dist),
            colour      = sig.col,
            linewidth   = 0.4,
            linetype    = "dashed",
            alpha       = 0.55,
            inherit.aes = FALSE) +
        ggplot2::geom_text(
            data        = sig.df,
            ggplot2::aes(x = .data$dist, label = .data$label),
            y           = Inf,
            colour      = sig.col,
            size        = 2.5,
            vjust       = 1.2,
            hjust       = 0.5,
            lineheight  = 0.85,
            inherit.aes = FALSE)
}

# =============================================================================
# Shared engine: .build_mv_effects_df
#
# Builds the data frame for the multivariate effects lollipop plot (Option A).
# One row per (QTL x trial combination):
#   MAIN QTL  -- one row,  trial = "all trials"
#   INT QTL   -- one row per trait level, trial = trait level name
#
# Columns returned:
#   qtl_key   -- "X.CHR.IDX" key
#   qtl_label -- facet strip label "CHR . POS cM  [MAIN|INT]"
#   y_row     -- y-axis factor: "all trials" or trial name
#   chr / pos_cM -- genomic position
#   effect / se / lo / hi -- estimate and +/- 1 SE
#   perc_var  -- % variance (same value for all rows of same QTL)
#   trial     -- "all trials" (MAIN) or trait level name (INT)
#   is_int    -- logical: TRUE = interaction QTL
#   show_pv   -- logical: TRUE on the row with largest |effect| per QTL
#                (this row gets the % var annotation)
# =============================================================================
.build_mv_effects_df <- function(object, genObj) {
    sigma2   <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4L)
        sigma2 <- 1
    gen.type     <- object$QTL$type
    trait.levels <- object$QTL$trait.levels

    # Coefficients from final pruned model (identical extraction to summary())
    all.fc <- object$coefficients$fixed
    all.vc <- object$vcoeff$fixed
    zind   <- grep("X\\.", rownames(all.fc))
    qtle   <- setNames(rev(all.fc[zind, 1L]), rev(rownames(all.fc)[zind]))
    veff   <- rev(all.vc[zind])
    se     <- sqrt(veff * sigma2)

    # Map each coefficient to its QTL key and trial label
    enams     <- names(qtle)
    qtl.x     <- sub("^.*:(X\\..*)$", "\\1", enams)
    qtl.x[!grepl(":", enams)] <- enams[!grepl(":", enams)]
    trait.lab <- rep("MAIN", length(enams))
    int.rows  <- grepl(":", enams)
    trait.lab[int.rows] <- sub("^(.*):X\\..*$", "\\1", enams[int.rows])
    prefix    <- paste0(object$QTL$Trait, "_")
    trait.lab <- gsub(prefix, "", trait.lab, fixed = TRUE)

    # % variance and positions from summary()
    summ    <- summary(object, genObj, LOD = FALSE)
    # With Trait column prepended: interval pos = col 6, marker pos = col 4
    pos_col <- if (gen.type == "interval") 6L else 4L
    chr_summ <- as.character(summ[["Chromosome"]])
    pos_summ <- as.numeric(summ[, pos_col])

    # Position info from getQTL()
    orig.qtl       <- object$QTL$qtl
    object$QTL$qtl <- sub("^X\\.", "Chr.", qtl.x)
    qtlm           <- getQTL(object, genObj)
    object$QTL$qtl <- orig.qtl

    # Unique QTL in genome order
    unique.qtl.x <- unique(qtl.x)
    chr_vals <- vapply(unique.qtl.x, function(uq) {
        i <- which(qtl.x == uq)[1L]
        match(as.character(qtlm[i, 1L]), names(genObj$geno))
    }, integer(1L))
    pos_vals <- vapply(unique.qtl.x, function(uq) {
        i <- which(qtl.x == uq)[1L]
        as.numeric(if (gen.type == "interval") qtlm[i, 6L] else qtlm[i, 4L])
    }, numeric(1L))
    uq_ordered <- unique.qtl.x[order(chr_vals * 1e6 + pos_vals)]

    rows <- lapply(uq_ordered, function(uq) {
        idx    <- which(qtl.x == uq)
        is.int <- any(trait.lab[idx] != "MAIN")

        i_first <- idx[1L]
        chr     <- as.character(qtlm[i_first, 1L])
        pos_cM  <- as.numeric(if (gen.type == "interval") qtlm[i_first, 6L]
                              else                         qtlm[i_first, 4L])
        tag       <- if (is.int) "[INT]" else "[MAIN]"
        qtl_label <- paste0(chr, " \u00b7 ", round(pos_cM, 1), " cM  ", tag)

        # % variance: first matching row by chr + cM position
        j  <- which(chr_summ == chr & abs(pos_summ - pos_cM) < 0.1)
        pv <- if (length(j) > 0L) summ[["Perc.Var"]][j[1L]] else NA_real_

        if (!is.int) {
            # MAIN: one row
            data.frame(
                qtl_key   = uq,
                qtl_label = qtl_label,
                y_row     = "all trials",
                chr       = chr,
                pos_cM    = pos_cM,
                effect    = qtle[idx[1L]],
                se        = se[idx[1L]],
                lo        = qtle[idx[1L]] - se[idx[1L]],
                hi        = qtle[idx[1L]] + se[idx[1L]],
                perc_var  = pv,
                trial     = "all trials",
                is_int    = FALSE,
                stringsAsFactors = FALSE
            )
        } else {
            # INT: one row per trial in trait.levels order
            do.call(rbind, lapply(trait.levels, function(tname) {
                i_t <- idx[trait.lab[idx] == tname]
                if (length(i_t) == 0L) return(NULL)
                eff <- qtle[i_t]; se_t <- se[i_t]
                data.frame(
                    qtl_key   = uq,
                    qtl_label = qtl_label,
                    y_row     = tname,
                    chr       = chr,
                    pos_cM    = pos_cM,
                    effect    = eff,
                    se        = se_t,
                    lo        = eff - se_t,
                    hi        = eff + se_t,
                    perc_var  = pv,
                    trial     = tname,
                    is_int    = TRUE,
                    stringsAsFactors = FALSE
                )
            }))
        }
    })

    df <- do.call(rbind, rows)

    # y_row factor: "all trials" first (bottom), then trait levels top to bottom
    df$y_row <- factor(df$y_row, levels = c(rev(trait.levels), "all trials"))

    # qtl_label factor: genome order (controls facet top-to-bottom order)
    qtl_labels_ord <- vapply(uq_ordered,
                             function(uq) df$qtl_label[df$qtl_key == uq][1L],
                             character(1L))
    df$qtl_label <- factor(df$qtl_label, levels = qtl_labels_ord)

    # Flag the row with the largest |effect| per QTL for the % var annotation
    df$show_pv <- FALSE
    for (uq in uq_ordered) {
        i_grp          <- which(df$qtl_key == uq)
        i_max          <- i_grp[which.max(abs(df$effect[i_grp]))]
        df$show_pv[i_max] <- TRUE
    }
    df
}

# =============================================================================
# Shared engine: .plot_mv_effects
#
# Horizontal lollipop plot for multivariate QTL/marker effects (Option A).
#
# Layout:
#   facet_grid(qtl_label ~ ., scales="free_y", space="free_y") -- one panel
#     per QTL; panel height proportional to number of trial rows
#   y-axis  : trial names ("all trials" for MAIN, level names for INT)
#   x-axis  : additive effect estimate
#   Colour  : MAIN rows in grey; INT rows coloured by trial (Dark 2 palette)
#             Legend shows trial names only (breaks = trait.levels)
#   % var   : annotated at the tip of the row with the largest |effect| per QTL
#   Facet strips (right): "CHR . POS cM [MAIN]" / "[INT]" labels, angle=0
# =============================================================================
.plot_mv_effects <- function(edf, object) {
    # Direction colour: firebrick = A allele favoured (+1),
    #                   steelblue = B allele favoured (-1)
    # Matches the univariate effects plot convention.
    edf$dir_col    <- ifelse(edf$effect >= 0, "firebrick", "steelblue")

    # % variance annotation: one label per QTL on the max-|effect| row
    pv_df             <- edf[edf$show_pv, ]
    pv_df$label_x     <- ifelse(pv_df$effect >= 0, pv_df$hi, pv_df$lo)
    pv_df$label_hjust <- ifelse(pv_df$effect >= 0, -0.35, 1.35)
    pv_df$pv_label    <- paste0(pv_df$perc_var, "%")

    # QTL label data frame: one row per facet, for top-right annotation
    qtl_label_df <- data.frame(
        qtl_label = levels(edf$qtl_label),
        stringsAsFactors = FALSE
    )
    qtl_label_df$qtl_label <- factor(qtl_label_df$qtl_label,
                                      levels = levels(edf$qtl_label))

    x_sub <- paste0(
        "Additive Effect  \u2014  ",
        "firebrick = A allele favoured (+1),  ",
        "steelblue = B allele favoured (\u22121)")

    ggplot2::ggplot(edf,
        ggplot2::aes(y = .data$y_row, x = .data$effect,
                     colour = .data$dir_col)) +

        ggplot2::facet_grid(qtl_label ~ .,
                            scales = "free_y", space = "free_y") +

        # Zero reference line
        ggplot2::geom_vline(xintercept = 0, colour = "grey50",
                            linewidth = 0.5, linetype = "dashed") +

        # Stems (0 to estimate)
        ggplot2::geom_segment(
            ggplot2::aes(y    = .data$y_row, yend = .data$y_row,
                         x    = 0,           xend = .data$effect),
            linewidth = 0.8, show.legend = FALSE) +

        # +/- 1 SE error bars
        ggplot2::geom_errorbar(
            ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
            width = 0.3, linewidth = 0.7, show.legend = FALSE) +

        # Lollipop heads
        ggplot2::geom_point(size = 3.5, show.legend = FALSE) +

        # % variance label at the tip of the largest-|effect| row per QTL
        ggplot2::geom_text(
            data        = pv_df,
            ggplot2::aes(x     = .data$label_x,
                         y     = .data$y_row,
                         label = .data$pv_label,
                         hjust = .data$label_hjust),
            size        = 3,
            colour      = "grey30",
            show.legend = FALSE,
            inherit.aes = FALSE) +

        # QTL label: top-right corner of each panel, no strip box
        ggplot2::geom_text(
            data        = qtl_label_df,
            ggplot2::aes(label = .data$qtl_label),
            x           = Inf,
            y           = Inf,
            hjust       = 1.05,
            vjust       = 1.2,
            size        = 2.8,
            colour      = "grey20",
            fontface    = "bold",
            inherit.aes = FALSE) +

        ggplot2::scale_colour_identity() +

        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.2, 0.2))) +

        ggplot2::xlab(x_sub) +
        ggplot2::ylab("") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter() +
        ggplot2::theme(
            strip.background = ggplot2::element_blank(),
            strip.text       = ggplot2::element_blank(),
            panel.spacing    = ggplot2::unit(0.8, "lines"),
            axis.text.y      = ggplot2::element_text(size = 8),
            axis.title.x     = ggplot2::element_text(size = 7.5,
                                                      colour = "grey40",
                                                      margin = ggplot2::margin(t = 8)),
            legend.position  = "none"
        )
}
