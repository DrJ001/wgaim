# =============================================================================
# checkPanel.R
# Diagnostic check of a raw marker panel before filtering and preparation.
#
# checkPanel() is purely diagnostic — it reports data quality issues without
# modifying the data. Use filterPanel() to act on the findings, then
# primePanel() to prepare the filtered data for analysis.
# =============================================================================

#' Check a Marker Panel for Data Quality Issues
#'
#' @description
#' Performs a comprehensive diagnostic check of a raw marker genotype matrix
#' and its associated genetic map, reporting data quality issues without
#' modifying the data.  The results guide the choice of thresholds for the
#' subsequent \code{\link{filterPanel}} step.
#'
#' The following checks are performed:
#' \describe{
#'   \item{Encoding validation}{Confirms that genotype values lie within the
#'     expected range for the stated \code{encoding}.}
#'   \item{Map consistency}{Identifies markers present in \code{geno} but
#'     absent from \code{map} (unplaceable) and vice versa.}
#'   \item{Duplicate lines}{Lines with identical genotype profiles.}
#'   \item{Duplicate markers}{Markers with identical genotype profiles across
#'     all lines.}
#'   \item{Marker missingness}{Per-marker proportion of missing genotypes.}
#'   \item{Line missingness}{Per-line proportion of missing genotypes.}
#'   \item{MAF distribution}{Per-marker minor allele frequency; summary
#'     statistics and a table showing how many markers would be removed at
#'     standard thresholds (0.01, 0.02, 0.05, 0.10).}
#'   \item{Monomorphic markers}{Markers with MAF = 0 that carry no
#'     information.}
#'   \item{Heterozygosity}{Per-line proportion of heterozygous calls.  Excess
#'     heterozygosity may indicate mixed or mislabelled samples.}
#'   \item{Chromosome coverage}{Number of markers and cM range per
#'     chromosome.}
#' }
#'
#' @param geno A numeric matrix (\code{lines x markers}) with row names
#'   identifying lines, or a \code{data.frame} with a column named by
#'   \code{id} identifying lines and remaining columns containing marker
#'   genotypes.
#' @param map A \code{data.frame} containing the genetic map with at minimum
#'   three columns: marker names (\code{map.id}), chromosome identifiers
#'   (\code{map.chr}), and cM positions (\code{map.pos}).
#' @param id Character string naming the line identifier column when
#'   \code{geno} is a \code{data.frame}.  Ignored when \code{geno} is a
#'   matrix.  Default \code{"id"}.
#' @param map.id Character string naming the marker column in \code{map}.
#'   Default \code{"marker"}.
#' @param map.chr Character string naming the chromosome column in \code{map}.
#'   Default \code{"chr"}.
#' @param map.pos Character string naming the cM position column in
#'   \code{map}.  Default \code{"pos"}.
#' @param encoding Character string specifying the genotype encoding:
#'   \code{"012"} (default, allele counts 0/1/2 or dosage in [0,2]) or
#'   \code{"pm1"} (already in additive \eqn{\pm 1} coding).
#'
#' @return An object of class \code{"checkPanel"} — a list containing:
#' \describe{
#'   \item{\code{$geno}}{The original (unmodified) genotype matrix.}
#'   \item{\code{$map}}{The original (unmodified) map data frame.}
#'   \item{\code{$encoding}, \code{$id}, \code{$map.id}, \code{$map.chr},
#'     \code{$map.pos}}{Arguments as supplied, carried through for
#'     \code{\link{filterPanel}}.}
#'   \item{\code{$n.lines}, \code{$n.markers}, \code{$n.chr}}{Dimensions of
#'     the panel.}
#'   \item{\code{$encoding.ok}}{Logical; \code{FALSE} if values lie outside
#'     the expected range.}
#'   \item{\code{$consistency}}{List with \code{$not_in_map} and
#'     \code{$not_in_geno}: character vectors of marker names.}
#'   \item{\code{$dup.lines}}{Character vector of duplicate line IDs (the
#'     second and subsequent copies).}
#'   \item{\code{$dup.markers}}{Character vector of duplicate marker names
#'     (the second and subsequent copies).}
#'   \item{\code{$miss.marker}}{Named numeric vector of per-marker missing
#'     rates.}
#'   \item{\code{$miss.line}}{Named numeric vector of per-line missing
#'     rates.}
#'   \item{\code{$maf}}{Named numeric vector of per-marker minor allele
#'     frequencies (computed on markers present in both \code{geno} and
#'     \code{map}).}
#'   \item{\code{$maf.table}}{Data frame showing the number and percentage of
#'     markers removed at standard MAF thresholds.}
#'   \item{\code{$monomorphic}}{Character vector of monomorphic marker
#'     names (MAF = 0).}
#'   \item{\code{$het}}{Named numeric vector of per-line heterozygosity
#'     rates.}
#'   \item{\code{$chr.coverage}}{Data frame with columns \code{chr},
#'     \code{n.markers}, \code{cM.min}, \code{cM.max}.}
#' }
#'
#' @seealso \code{\link{filterPanel}}, \code{\link{primePanel}}
#' @export
checkPanel <- function(geno, map,
                       id      = "id",
                       map.id  = "marker",
                       map.chr = "chr",
                       map.pos = "pos",
                       encoding = "012") {

    encoding <- match.arg(encoding, c("012", "pm1"))

    # -------------------------------------------------------------------------
    # 1.  Extract raw matrix and line IDs
    # -------------------------------------------------------------------------
    if (is.data.frame(geno)) {
        if (!id %in% names(geno))
            stop("Column '", id, "' not found in geno data frame.")
        ids      <- as.character(geno[[id]])
        geno.mat <- as.matrix(geno[, !names(geno) %in% id, drop = FALSE])
        rownames(geno.mat) <- ids
    } else {
        geno.mat <- as.matrix(geno)
        ids      <- rownames(geno.mat)
        if (is.null(ids))
            stop("geno matrix must have row names identifying lines.")
    }
    storage.mode(geno.mat) <- "numeric"

    for (col in c(map.id, map.chr, map.pos))
        if (!col %in% names(map))
            stop("Column '", col, "' not found in map.")

    n.lines   <- nrow(geno.mat)
    n.markers <- ncol(geno.mat)

    # -------------------------------------------------------------------------
    # 2.  Encoding validation (on raw values before any transformation)
    # -------------------------------------------------------------------------
    rng <- range(geno.mat, na.rm = TRUE)
    encoding.ok <- if (encoding == "012")
        rng[1] >= -0.05 && rng[2] <= 2.05   # allow small floating-point slack
    else
        rng[1] >= -1.05 && rng[2] <= 1.05
    encoding.warn <- if (!encoding.ok)
        sprintf("Values outside expected range for encoding='%s': [%.3f, %.3f].",
                encoding, rng[1], rng[2])
    else
        NULL

    # -------------------------------------------------------------------------
    # 3.  Map consistency
    # -------------------------------------------------------------------------
    map.markers  <- as.character(map[[map.id]])
    geno.markers <- colnames(geno.mat)
    if (is.null(geno.markers))
        stop("geno must have column names matching marker names in map.")

    not_in_map  <- setdiff(geno.markers, map.markers)
    not_in_geno <- setdiff(map.markers, geno.markers)

    # Restrict to common markers for all subsequent checks
    common   <- intersect(geno.markers, map.markers)
    geno.chk <- geno.mat[, common, drop = FALSE]
    map.chk  <- map[map[[map.id]] %in% common, , drop = FALSE]

    # -------------------------------------------------------------------------
    # 4.  Duplicate lines (string-signature approach; NA treated as distinct)
    # -------------------------------------------------------------------------
    line_sigs <- apply(geno.chk, 1L, function(x)
        paste(ifelse(is.na(x), "NA", as.character(x)), collapse = "|"))
    dup.lines <- ids[duplicated(line_sigs)]

    # -------------------------------------------------------------------------
    # 5.  Duplicate markers
    # -------------------------------------------------------------------------
    mark_sigs <- apply(geno.chk, 2L, function(x)
        paste(ifelse(is.na(x), "NA", as.character(x)), collapse = "|"))
    dup.markers <- common[duplicated(mark_sigs)]

    # -------------------------------------------------------------------------
    # 6.  Marker missingness
    # -------------------------------------------------------------------------
    miss.marker <- colMeans(is.na(geno.chk))

    # -------------------------------------------------------------------------
    # 7.  Line missingness
    # -------------------------------------------------------------------------
    miss.line <- rowMeans(is.na(geno.chk))

    # -------------------------------------------------------------------------
    # 8.  MAF  (computed on raw values before encoding transformation)
    # -------------------------------------------------------------------------
    if (encoding == "012") {
        # values in [0, 2]: p = colMean / 2
        col_means <- colMeans(geno.chk, na.rm = TRUE)
        p         <- col_means / 2
    } else {
        # values in [-1, 1]: p = (colMean + 1) / 2
        col_means <- colMeans(geno.chk, na.rm = TRUE)
        p         <- (col_means + 1) / 2
    }
    maf_vec     <- pmin(p, 1 - p)
    monomorphic <- names(maf_vec)[!is.na(maf_vec) & maf_vec == 0]

    maf_thresh  <- c(0.01, 0.02, 0.05, 0.10)
    maf.table   <- data.frame(
        threshold   = maf_thresh,
        n.removed   = sapply(maf_thresh, function(t)
            sum(!is.na(maf_vec) & maf_vec < t)),
        pct.removed = sapply(maf_thresh, function(t)
            round(100 * mean(!is.na(maf_vec) & maf_vec < t), 1)),
        stringsAsFactors = FALSE
    )

    # -------------------------------------------------------------------------
    # 9.  Heterozygosity per line
    # -------------------------------------------------------------------------
    het <- if (encoding == "012") {
        # heterozygous call = 1 in 0/1/2 coding
        rowMeans(geno.chk == 1L, na.rm = TRUE)
    } else {
        # heterozygous call = 0 in pm1 coding
        rowMeans(geno.chk == 0L, na.rm = TRUE)
    }
    names(het) <- ids

    # -------------------------------------------------------------------------
    # 10.  Chromosome coverage
    # -------------------------------------------------------------------------
    chr_vals <- as.character(map.chk[[map.chr]])
    chrs     <- unique(chr_vals)
    chr.coverage <- do.call("rbind", lapply(chrs, function(ch) {
        idx  <- chr_vals == ch
        mpos <- as.numeric(map.chk[[map.pos]][idx])
        data.frame(
            chr       = ch,
            n.markers = sum(idx),
            cM.min    = min(mpos, na.rm = TRUE),
            cM.max    = max(mpos, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }))

    # -------------------------------------------------------------------------
    # 11.  Assemble and return
    # -------------------------------------------------------------------------
    structure(
        list(
            geno         = geno,
            map          = map,
            encoding     = encoding,
            id           = id,
            map.id       = map.id,
            map.chr      = map.chr,
            map.pos      = map.pos,
            n.lines      = n.lines,
            n.markers    = n.markers,
            n.chr        = length(chrs),
            encoding.ok  = encoding.ok,
            encoding.warn = encoding.warn,
            consistency  = list(not_in_map  = not_in_map,
                                not_in_geno = not_in_geno),
            dup.lines    = dup.lines,
            dup.markers  = dup.markers,
            miss.marker  = miss.marker,
            miss.line    = miss.line,
            maf          = maf_vec,
            maf.table    = maf.table,
            monomorphic  = monomorphic,
            het          = het,
            chr.coverage = chr.coverage
        ),
        class = "checkPanel"
    )
}

# =============================================================================
# print.checkPanel
# =============================================================================

#' Print a checkPanel summary
#'
#' @param x A \code{"checkPanel"} object from \code{\link{checkPanel}}.
#' @param miss.line.thresh   Threshold for flagging lines in the missingness
#'   report.  Default \code{0.20}.
#' @param miss.marker.thresh Threshold for flagging markers in the missingness
#'   report.  Default \code{0.20}.
#' @param het.thresh Threshold for flagging lines in the heterozygosity
#'   report.  Default \code{0.10}.
#' @param \dots Currently ignored.
#' @return \code{x} invisibly.
#' @exportS3Method
print.checkPanel <- function(x,
                              miss.line.thresh   = 0.20,
                              miss.marker.thresh = 0.20,
                              het.thresh         = 0.10,
                              ...) {
    .sep  <- function(n = 52) cat(strrep("=", n), "\n", sep = "")
    .hdr  <- function(s)      cat("\n", s, "\n", sep = "")
    .flag <- function(n, what) if (n > 0L) cat("  *** ", n, what, "\n")

    cat("\nPanel Check Summary\n")
    .sep()
    cat(sprintf("  Lines: %d  |  Markers: %d  |  Chromosomes: %d\n",
                x$n.lines, x$n.markers, x$n.chr))
    cat(sprintf("  Encoding: %s\n", x$encoding))

    # Encoding
    .hdr("Encoding validation:")
    if (x$encoding.ok)
        cat("  OK -- values within expected range.\n")
    else
        cat("  WARNING:", x$encoding.warn, "\n")

    # Map consistency
    .hdr("Map consistency:")
    cat(sprintf("  Markers in geno, not in map : %d\n",
                length(x$consistency$not_in_map)))
    cat(sprintf("  Markers in map,  not in geno: %d\n",
                length(x$consistency$not_in_geno)))

    # Duplicates
    .hdr("Duplicates:")
    cat(sprintf("  Duplicate lines  : %d\n", length(x$dup.lines)))
    cat(sprintf("  Duplicate markers: %d\n", length(x$dup.markers)))

    # Missingness
    .hdr("Missingness:")
    ml <- x$miss.line;   mm <- x$miss.marker
    cat(sprintf("  Lines   -- Mean: %.1f%%  Max: %.1f%%",
                100 * mean(ml), 100 * max(ml)))
    n_fl <- sum(ml > miss.line.thresh)
    if (n_fl > 0L) cat(sprintf("  (%d above %.0f%% threshold)", n_fl,
                                100 * miss.line.thresh))
    cat("\n")
    cat(sprintf("  Markers -- Mean: %.1f%%  Max: %.1f%%",
                100 * mean(mm), 100 * max(mm)))
    n_fm <- sum(mm > miss.marker.thresh)
    if (n_fm > 0L) cat(sprintf("  (%d above %.0f%% threshold)", n_fm,
                                100 * miss.marker.thresh))
    cat("\n")

    # MAF
    .hdr("MAF distribution:")
    maf <- x$maf[!is.na(x$maf)]
    cat(sprintf("  Min: %.3f  Median: %.3f  Mean: %.3f  Max: %.3f\n",
                min(maf), median(maf), mean(maf), max(maf)))
    cat(sprintf("  Monomorphic markers (MAF = 0): %d\n",
                length(x$monomorphic)))
    cat("  Markers removed at standard thresholds:\n")
    for (i in seq_len(nrow(x$maf.table))) {
        r <- x$maf.table[i, ]
        cat(sprintf("    MAF < %.2f : %4d markers (%4.1f%%)\n",
                    r$threshold, r$n.removed, r$pct.removed))
    }

    # Heterozygosity
    .hdr("Heterozygosity (per line):")
    ht <- x$het
    cat(sprintf("  Mean: %.3f  Max: %.3f",  mean(ht), max(ht)))
    n_fh <- sum(ht > het.thresh)
    if (n_fh > 0L) cat(sprintf("  (%d lines above %.0f%% threshold)",
                                n_fh, 100 * het.thresh))
    cat("\n")

    # Chromosome coverage
    .hdr("Chromosome coverage:")
    cc <- x$chr.coverage
    cat(sprintf("  %-12s  %8s  %s\n", "Chromosome", "Markers", "Range (cM)"))
    cat(sprintf("  %-12s  %8s  %s\n", strrep("-", 12),
                strrep("-", 7), strrep("-", 18)))
    for (i in seq_len(nrow(cc))) {
        cat(sprintf("  %-12s  %8d  %.1f \u2013 %.1f cM\n",
                    cc$chr[i], cc$n.markers[i], cc$cM.min[i], cc$cM.max[i]))
    }
    cat("\n")
    invisible(x)
}
