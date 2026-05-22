# =============================================================================
# filterPanel.R
# Filter a raw marker panel based on data quality thresholds.
#
# filterPanel() applies filters in the statistically correct order:
#   1. Map consistency  -- drop unplaceable markers
#   2. Duplicate lines  -- drop exact genotype duplicates
#   3. Duplicate markers-- drop redundant marker columns
#   4. Marker missingness -- drop high-missing markers before assessing lines
#   5. Line missingness   -- drop high-missing lines (on cleaned marker set)
#   6. MAF               -- drop rare/monomorphic markers on cleaned line set
#   7. Heterozygosity     -- flag/drop lines with excess heterozygosity
#
# The function is an S3 generic:
#   filterPanel.default()     -- accepts raw geno + map
#   filterPanel.checkPanel()  -- accepts a checkPanel object (reuses stored args)
# =============================================================================

#' Filter a Marker Panel by Data Quality Thresholds
#'
#' @description
#' Applies a sequential set of data quality filters to a raw marker genotype
#' matrix and its associated genetic map.  Filters are applied in a
#' statistically principled order so that each step operates on data already
#' cleaned by the preceding steps:
#'
#' \enumerate{
#'   \item \strong{Map consistency} — markers absent from \code{map} are
#'     dropped; these cannot be placed on the genome.
#'   \item \strong{Duplicate lines} — lines with identical genotype profiles
#'     are dropped (second and subsequent copies).  Performed before
#'     population-level statistics so allele frequencies are not inflated.
#'   \item \strong{Duplicate markers} — markers with identical genotype
#'     profiles are dropped (second and subsequent copies).  Performed before
#'     missingness calculations.
#'   \item \strong{Marker missingness} — markers with a missing rate above
#'     \code{miss.marker} are removed \emph{before} line missingness is
#'     assessed, because a failed genotyping assay can make many lines appear
#'     more missing than they really are.
#'   \item \strong{Line missingness} — lines with a missing rate above
#'     \code{miss.line} are removed, calculated on the cleaned marker set.
#'   \item \strong{MAF} — markers with minor allele frequency below
#'     \code{maf} are removed, calculated on the cleaned line set so that
#'     allele frequencies reflect the retained population.
#'   \item \strong{Heterozygosity} — lines with a heterozygosity rate above
#'     \code{het} are removed.  Performed last on the MAF-filtered data.
#'     Skipped when \code{het = NULL} (default).
#' }
#'
#' @param x For the \code{checkPanel} method: a \code{"checkPanel"} object
#'   from \code{\link{checkPanel}}.  The stored \code{geno}, \code{map},
#'   \code{encoding}, and column-name arguments are used automatically.
#' @param geno A numeric matrix (\code{lines x markers}) with row names
#'   identifying lines, or a \code{data.frame} with a line identifier column.
#'   (Default method only.)
#' @param map A \code{data.frame} containing the genetic map.
#'   (Default method only.)
#' @param id Character string naming the line identifier column when
#'   \code{geno} is a \code{data.frame}.  Default \code{"id"}.
#' @param map.id Character string naming the marker column in \code{map}.
#'   Default \code{"marker"}.
#' @param map.chr Character string naming the chromosome column in \code{map}.
#'   Default \code{"chr"}.
#' @param map.pos Character string naming the cM position column in
#'   \code{map}.  Default \code{"pos"}.
#' @param encoding Character string specifying the genotype encoding:
#'   \code{"012"} (default) or \code{"pm1"}.
#' @param miss.line Numeric scalar.  Lines with a missing rate above this
#'   threshold are removed.  Default \code{0.20}.  Set to \code{NULL} to
#'   skip.
#' @param miss.marker Numeric scalar.  Markers with a missing rate above this
#'   threshold are removed.  Default \code{0.20}.  Set to \code{NULL} to
#'   skip.
#' @param maf Numeric scalar.  Markers with minor allele frequency below this
#'   threshold are removed.  Default \code{0.05}.  Set to \code{NULL} to
#'   skip.
#' @param dup.lines Logical.  If \code{TRUE} (default), duplicate lines are
#'   removed.
#' @param dup.markers Logical.  If \code{FALSE} (default), duplicate markers
#'   are \emph{not} removed.  Set to \code{TRUE} to drop them.
#' @param het Numeric scalar.  Lines with a heterozygosity rate above this
#'   threshold are removed.  \code{NULL} (default) skips this filter.
#' @param \dots Currently ignored.
#'
#' @return An object of class \code{"filteredPanel"} — a list containing:
#' \describe{
#'   \item{\code{$geno}}{Filtered genotype matrix.}
#'   \item{\code{$map}}{Filtered map data frame.}
#'   \item{\code{$encoding}, \code{$id}, \code{$map.id}, \code{$map.chr},
#'     \code{$map.pos}}{Carried through for \code{\link{primePanel}}.}
#'   \item{\code{$thresholds}}{Named list of the threshold arguments used.}
#'   \item{\code{$removed}}{Named list with one element per filter step,
#'     each a character vector of the removed line or marker names.}
#'   \item{\code{$n.original}}{Named integer: \code{lines} and
#'     \code{markers} before filtering.}
#'   \item{\code{$n.final}}{Named integer: \code{lines} and \code{markers}
#'     after filtering.}
#' }
#'
#' @seealso \code{\link{checkPanel}}, \code{\link{primePanel}}
#' @export
filterPanel <- function(x, ...) UseMethod("filterPanel")

#' @rdname filterPanel
#' @export
filterPanel.default <- function(x, map,
                                 id           = "id",
                                 map.id       = "marker",
                                 map.chr      = "chr",
                                 map.pos      = "pos",
                                 encoding     = "012",
                                 miss.line    = 0.20,
                                 miss.marker  = 0.20,
                                 maf          = 0.05,
                                 dup.lines    = TRUE,
                                 dup.markers  = FALSE,
                                 het          = NULL,
                                 ...) {
    .filterPanel_core(geno = x, map = map, id = id, map.id = map.id,
                      map.chr = map.chr, map.pos = map.pos,
                      encoding = encoding, miss.line = miss.line,
                      miss.marker = miss.marker, maf = maf,
                      dup.lines = dup.lines, dup.markers = dup.markers,
                      het = het)
}

#' @rdname filterPanel
#' @export
filterPanel.checkPanel <- function(x,
                                    miss.line    = 0.20,
                                    miss.marker  = 0.20,
                                    maf          = 0.05,
                                    dup.lines    = TRUE,
                                    dup.markers  = FALSE,
                                    het          = NULL,
                                    ...) {
    .filterPanel_core(geno = x$geno, map = x$map, id = x$id,
                      map.id = x$map.id, map.chr = x$map.chr,
                      map.pos = x$map.pos, encoding = x$encoding,
                      miss.line = miss.line, miss.marker = miss.marker,
                      maf = maf, dup.lines = dup.lines,
                      dup.markers = dup.markers, het = het)
}

# =============================================================================
# Internal workhorse
# =============================================================================
.filterPanel_core <- function(geno, map, id, map.id, map.chr, map.pos,
                               encoding, miss.line, miss.marker, maf,
                               dup.lines, dup.markers, het) {

    encoding <- match.arg(encoding, c("012", "pm1"))

    # Extract matrix and IDs
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

    n.orig.lines   <- nrow(geno.mat)
    n.orig.markers <- ncol(geno.mat)
    removed        <- list()

    # -------------------------------------------------------------------------
    # Step 1: Map consistency
    # -------------------------------------------------------------------------
    map.markers  <- as.character(map[[map.id]])
    geno.markers <- colnames(geno.mat)
    not_in_map   <- setdiff(geno.markers, map.markers)
    if (length(not_in_map)) {
        geno.mat <- geno.mat[, !colnames(geno.mat) %in% not_in_map, drop = FALSE]
        map      <- map[map[[map.id]] %in% colnames(geno.mat), , drop = FALSE]
    }
    removed$map_consistency <- not_in_map

    # -------------------------------------------------------------------------
    # Step 2: Duplicate lines
    # -------------------------------------------------------------------------
    if (dup.lines) {
        line_sigs <- apply(geno.mat, 1L, function(x)
            paste(ifelse(is.na(x), "NA", as.character(x)), collapse = "|"))
        dup_idx  <- duplicated(line_sigs)
        dup_nms  <- rownames(geno.mat)[dup_idx]
        if (length(dup_nms)) {
            geno.mat <- geno.mat[!dup_idx, , drop = FALSE]
        }
        removed$dup_lines <- dup_nms
    } else {
        removed$dup_lines <- character(0L)
    }

    # -------------------------------------------------------------------------
    # Step 3: Duplicate markers
    # -------------------------------------------------------------------------
    if (dup.markers) {
        mark_sigs <- apply(geno.mat, 2L, function(x)
            paste(ifelse(is.na(x), "NA", as.character(x)), collapse = "|"))
        dup_midx  <- duplicated(mark_sigs)
        dup_mnms  <- colnames(geno.mat)[dup_midx]
        if (length(dup_mnms)) {
            geno.mat <- geno.mat[, !dup_midx, drop = FALSE]
            map      <- map[map[[map.id]] %in% colnames(geno.mat), , drop = FALSE]
        }
        removed$dup_markers <- dup_mnms
    } else {
        removed$dup_markers <- character(0L)
    }

    # -------------------------------------------------------------------------
    # Step 4: Marker missingness
    # -------------------------------------------------------------------------
    if (!is.null(miss.marker)) {
        mm       <- colMeans(is.na(geno.mat))
        drop_m   <- names(mm)[mm > miss.marker]
        if (length(drop_m)) {
            geno.mat <- geno.mat[, !colnames(geno.mat) %in% drop_m, drop = FALSE]
            map      <- map[map[[map.id]] %in% colnames(geno.mat), , drop = FALSE]
        }
        removed$miss_marker <- drop_m
    } else {
        removed$miss_marker <- character(0L)
    }

    # -------------------------------------------------------------------------
    # Step 5: Line missingness
    # -------------------------------------------------------------------------
    if (!is.null(miss.line)) {
        ml     <- rowMeans(is.na(geno.mat))
        drop_l <- names(ml)[ml > miss.line]
        if (length(drop_l))
            geno.mat <- geno.mat[!rownames(geno.mat) %in% drop_l, , drop = FALSE]
        removed$miss_line <- drop_l
    } else {
        removed$miss_line <- character(0L)
    }

    # -------------------------------------------------------------------------
    # Step 6: MAF  (on cleaned line set)
    # -------------------------------------------------------------------------
    if (!is.null(maf) && maf > 0) {
        col_means <- colMeans(geno.mat, na.rm = TRUE)
        p         <- if (encoding == "012") col_means / 2
                     else (col_means + 1) / 2
        maf_vec   <- pmin(p, 1 - p)
        drop_maf  <- names(maf_vec)[!is.na(maf_vec) & maf_vec < maf]
        if (length(drop_maf)) {
            geno.mat <- geno.mat[, !colnames(geno.mat) %in% drop_maf, drop = FALSE]
            map      <- map[map[[map.id]] %in% colnames(geno.mat), , drop = FALSE]
        }
        removed$maf <- drop_maf
    } else {
        removed$maf <- character(0L)
    }

    # -------------------------------------------------------------------------
    # Step 7: Heterozygosity  (on MAF-filtered data)
    # -------------------------------------------------------------------------
    if (!is.null(het)) {
        het_line <- if (encoding == "012")
            rowMeans(geno.mat == 1L, na.rm = TRUE)
        else
            rowMeans(geno.mat == 0L, na.rm = TRUE)
        drop_het <- names(het_line)[!is.na(het_line) & het_line > het]
        if (length(drop_het))
            geno.mat <- geno.mat[!rownames(geno.mat) %in% drop_het, , drop = FALSE]
        removed$het <- drop_het
    } else {
        removed$het <- character(0L)
    }

    if (nrow(geno.mat) == 0L)
        warning("All lines were removed by the filters. ",
                "Consider relaxing one or more thresholds.")
    if (ncol(geno.mat) == 0L)
        warning("All markers were removed by the filters. ",
                "Consider relaxing one or more thresholds.")

    structure(
        list(
            geno       = geno.mat,
            map        = map,
            encoding   = encoding,
            id         = id,
            map.id     = map.id,
            map.chr    = map.chr,
            map.pos    = map.pos,
            thresholds = list(miss.line = miss.line, miss.marker = miss.marker,
                              maf = maf, dup.lines = dup.lines,
                              dup.markers = dup.markers, het = het),
            removed    = removed,
            n.original = c(lines = n.orig.lines, markers = n.orig.markers),
            n.final    = c(lines   = nrow(geno.mat),
                           markers = ncol(geno.mat))
        ),
        class = "filteredPanel"
    )
}

# =============================================================================
# print.filteredPanel
# =============================================================================

#' Print a filteredPanel summary
#'
#' @param x A \code{"filteredPanel"} object from \code{\link{filterPanel}}.
#' @param \dots Currently ignored.
#' @return \code{x} invisibly.
#' @exportS3Method
print.filteredPanel <- function(x, ...) {
    cat("\nFilter Summary\n")
    cat(strrep("=", 52), "\n", sep = "")
    cat(sprintf("  Original : %d lines,  %d markers\n",
                x$n.original["lines"], x$n.original["markers"]))
    cat("\n")

    steps <- list(
        list(key = "map_consistency", label = "Map consistency",
             what = "marker"),
        list(key = "dup_lines",       label = "Duplicate lines",
             what = "line"),
        list(key = "dup_markers",     label = "Duplicate markers",
             what = "marker"),
        list(key = "miss_marker",
             label = sprintf("Marker missingness > %.0f%%",
                             100 * (x$thresholds$miss.marker %||% 0)),
             what = "marker"),
        list(key = "miss_line",
             label = sprintf("Line missingness > %.0f%%",
                             100 * (x$thresholds$miss.line %||% 0)),
             what = "line"),
        list(key = "maf",
             label = sprintf("MAF < %.2f", x$thresholds$maf %||% 0),
             what = "marker"),
        list(key = "het",
             label = sprintf("Heterozygosity > %.0f%%",
                             100 * (x$thresholds$het %||% 0)),
             what = "line")
    )

    for (i in seq_along(steps)) {
        s   <- steps[[i]]
        rem <- x$removed[[s$key]]
        if (is.null(rem)) {
            cat(sprintf("  Step %d [%-35s]: skipped\n", i, s$label))
        } else {
            cat(sprintf("  Step %d [%-35s]: removed %d %s(s)\n",
                        i, s$label, length(rem), s$what))
        }
    }

    cat("\n")
    cat(sprintf("  Remaining: %d lines,  %d markers\n",
                x$n.final["lines"], x$n.final["markers"]))
    cat("\n")
    invisible(x)
}

# Null-coalescing helper (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b
