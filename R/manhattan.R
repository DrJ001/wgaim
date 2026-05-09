# =============================================================================
# manhattan.R
# S3 generic and GWASAim method for Manhattan plots.
# Uses first-iteration outlier statistics to derive approximate single-marker
# p-values from the genome-wide conditional scan.
# =============================================================================

manhattan <- function(object, ...)
    UseMethod("manhattan")

manhattan.GWASAim <- function(object, panelObj, iter = 1,
                               chr = NULL, sig.col = "red",
                               pt.col  = c("steelblue", "grey50"),
                               pt.cex  = 0.6,
                               labels  = TRUE, ...) {
    if (missing(panelObj))
        stop("panelObj is a required argument.")
    if (!inherits(panelObj, "panel"))
        stop("panelObj must be of class \"panel\".")
    if (iter > length(object$QTL$diag$oint))
        stop("iter exceeds the number of analysis iterations stored.")

    # -------------------------------------------------------------------------
    # Build position data from panelObj
    # -------------------------------------------------------------------------
    chr.names <- names(panelObj$geno)
    if (!is.null(chr)) {
        if (!all(chr %in% chr.names))
            stop("Some chromosome names not found in panelObj.")
        chr.names <- chr
    }

    # Cumulative positions across chromosomes
    pos.list <- lapply(panelObj$geno[chr.names], function(el) el$map)
    chr.lens  <- sapply(pos.list, function(p) max(p) - min(p))
    chr.offsets <- c(0, cumsum(chr.lens + 5))   # 5 cM gap between chromosomes
    names(chr.offsets) <- c(chr.names, "end")

    marker.df <- do.call(rbind, lapply(seq_along(chr.names), function(i) {
        ch  <- chr.names[i]
        pos <- pos.list[[ch]]
        data.frame(
            marker   = names(pos),
            chr      = ch,
            pos.cM   = as.numeric(pos),
            cum.pos  = as.numeric(pos) - min(as.numeric(pos)) + chr.offsets[i],
            chr.idx  = i,
            stringsAsFactors = FALSE
        )
    }))

    # -------------------------------------------------------------------------
    # Extract outlier statistics and convert to -log10(p)
    # -------------------------------------------------------------------------
    oint <- object$QTL$diag$oint[[iter]]

    # oint is over all markers (with state=0 for excluded/selected)
    # Convert: oint ~ chi-sq(1), so p = (1 - pchisq(oint, 1)) / 2
    pvals           <- (1 - pchisq(oint, 1)) / 2
    pvals[pvals == 0] <- .Machine$double.xmin   # avoid -Inf
    log10p          <- -log10(pvals)

    # Match to marker.df by name
    mnames         <- paste0("Chr.", marker.df$chr, ".",
                              match(marker.df$marker,
                                    names(panelObj$geno[[marker.df$chr[1]]]$map)))
    # Build lookup from internal names to log10p
    stat.names     <- names(oint)
    marker.df$oint <- NA_real_
    for (k in seq_len(nrow(marker.df))) {
        ch  <- marker.df$chr[k]
        mkr <- marker.df$marker[k]
        idx <- which(names(panelObj$geno[[ch]]$map) == mkr)
        if (length(idx)) {
            key <- paste("Chr", ch, idx, sep = ".")
            if (key %in% stat.names)
                marker.df$oint[k] <- oint[key]
        }
    }
    marker.df$log10p <- ifelse(is.na(marker.df$oint), NA,
                               -log10(pmax((1 - pchisq(marker.df$oint, 1)) / 2,
                                           .Machine$double.xmin)))

    # Alternating chromosome colours
    marker.df$col <- pt.col[(marker.df$chr.idx %% 2) + 1]

    # -------------------------------------------------------------------------
    # Significance threshold line: -log10(TypeI.eff)
    # -------------------------------------------------------------------------
    thresh.log10p <- -log10(object$QTL$TypeI.eff)
    sig.markers   <- if (!is.null(object$QTL$qtl) && iter <= length(object$QTL$qtl))
        object$QTL$qtl[1:iter] else character(0)

    # -------------------------------------------------------------------------
    # ggplot2 Manhattan plot
    # -------------------------------------------------------------------------
    chr.mids <- sapply(chr.names, function(ch) {
        rows <- marker.df[marker.df$chr == ch, ]
        mean(range(rows$cum.pos, na.rm = TRUE))
    })

    gp <- ggplot(marker.df[!is.na(marker.df$log10p), ],
                 aes_string(x = "cum.pos", y = "log10p")) +
        geom_point(aes_string(colour = "col"), size = pt.cex, show.legend = FALSE) +
        scale_colour_identity() +
        geom_hline(yintercept = thresh.log10p, colour = sig.col,
                   linetype = "dashed", linewidth = 0.6) +
        scale_x_continuous(breaks = chr.mids, labels = chr.names) +
        labs(x = "Chromosome", y = expression(-log[10](p)),
             title = sprintf("Manhattan Plot (Iteration %d)", iter),
             caption = sprintf("Dashed line: Bonferroni threshold (p = %s)",
                               formatC(object$QTL$TypeI.eff, format = "e", digits = 2))) +
        theme_scatter()

    # Label significant markers from this iteration
    if (labels && length(sig.markers)) {
        sig.df <- do.call(rbind, lapply(sig.markers, function(qtl.name) {
            parts <- strsplit(qtl.name, "\\.")[[1]]
            ch    <- parts[2]; idx <- as.integer(parts[3])
            mkr   <- names(panelObj$geno[[ch]]$map)[idx]
            rows  <- marker.df[marker.df$chr == ch & marker.df$marker == mkr, ]
            if (nrow(rows)) rows[1, , drop = FALSE] else NULL
        }))
        if (!is.null(sig.df) && nrow(sig.df))
            gp <- gp + geom_point(data = sig.df,
                                   aes_string(x = "cum.pos", y = "log10p"),
                                   colour = sig.col, size = pt.cex * 2.5) +
                geom_text(data = sig.df,
                          aes_string(x = "cum.pos", y = "log10p", label = "marker"),
                          colour = sig.col, size = 2.5, vjust = -0.8, hjust = 0.5)
    }
    gp
}
