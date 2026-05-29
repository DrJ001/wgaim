# =============================================================================
# filterPanel.R
# Filter a raw marker panel based on data quality thresholds.
#
# filterPanel() is an S3 generic that accepts:
#   - A raw geno matrix / data.frame  (filterPanel.default)
#   - A "checkPanel"   object         (filterPanel.checkPanel)
#   - A "filteredPanel" object        (filterPanel.filteredPanel)
#
# Filtering is controlled by a single `steps` named list.  The list element
# order determines execution order.  The supported step names and their
# sentinel "skip" values are:
#
#   map         = FALSE  (logical; TRUE runs the step)
#   miss.marker = NULL   (numeric threshold or NULL to skip)
#   miss.line   = NULL   (numeric threshold or NULL to skip)
#   het.line    = NULL   (numeric threshold or NULL to skip)
#   het.marker  = NULL   (numeric threshold or NULL to skip)
#   dup.lines   = FALSE  (logical)
#   dup.markers = FALSE  (logical)
#   maf         = NULL   (numeric threshold or NULL to skip)
#
# Map consistency is always inserted as the very first step, even if the
# user omits it from a custom steps list.
#
# When the input is a "filteredPanel" object, a new pass is appended to
# $history rather than replacing it, so the full filtering audit trail is
# preserved across multiple calls.
# =============================================================================

# Valid step names (used for validation and label lookup)
.valid_step_names <- c("map", "miss.marker", "miss.line", "het.line",
                       "het.marker", "dup.lines", "dup.markers", "maf")

# Default steps list (canonical reference; matches old function-argument defaults)
.default_filter_steps <- list(
    map         = TRUE,
    miss.marker = 0.20,
    miss.line   = 0.20,
    het.line    = NULL,
    het.marker  = NULL,
    dup.lines   = TRUE,
    dup.markers = FALSE,
    maf         = 0.05
)

# =============================================================================
# S3 generic
# =============================================================================

#' Filter a Marker Panel by Data Quality Thresholds
#'
#' @description
#' Applies a user-defined sequence of data quality filters to a raw marker
#' genotype matrix and its associated genetic map.  Filtering is controlled
#' by a single \code{steps} named list whose \emph{element order} determines
#' the execution order.  The default \code{steps} list applies a statistically
#' principled workflow:
#'
#' \enumerate{
#'   \item \strong{Map consistency} (\code{map}) — markers absent from
#'     \code{map} are dropped; always executed first even if omitted from a
#'     custom \code{steps} list.
#'   \item \strong{Marker missingness} (\code{miss.marker}) — markers with a
#'     missing rate above the threshold are removed first, because a failed
#'     genotyping assay inflates apparent line-level missingness.
#'   \item \strong{Line missingness} (\code{miss.line}) — lines with a missing
#'     rate above the threshold are removed on the cleaned marker set.
#'   \item \strong{Line heterozygosity} (\code{het.line}) — lines with a
#'     heterozygosity rate above the threshold are removed.  Excess
#'     heterozygosity typically indicates a mislabelled or contaminated
#'     sample.  Skipped by default (\code{NULL}).
#'   \item \strong{Marker heterozygosity} (\code{het.marker}) — markers with
#'     a heterozygosity rate above the threshold are removed.  High per-marker
#'     heterozygosity suggests a paralogous locus or genotyping artefact.
#'     Skipped by default (\code{NULL}).
#'   \item \strong{Duplicate lines} (\code{dup.lines}) — lines with identical
#'     genotype profiles are dropped (second and subsequent copies), performed
#'     after quality filters so that clean copies are retained.
#'   \item \strong{Duplicate markers} (\code{dup.markers}) — markers with
#'     identical genotype profiles are dropped.  Skipped by default
#'     (\code{FALSE}).
#'   \item \strong{MAF} (\code{maf}) — markers with minor allele frequency
#'     below the threshold are removed last, so that allele frequencies are
#'     computed on the fully cleaned dataset.
#' }
#'
#' @section Custom step lists:
#' Supply a partial \code{steps} list to run only those steps, in the order
#' given.  Use \code{modifyList(formals(filterPanel.default)$steps, list(...))}
#' to adjust individual thresholds while keeping all other defaults.
#'
#' @section Multi-pass filtering:
#' Passing a \code{"filteredPanel"} object as \code{geno} runs an additional
#' filtering pass and \emph{appends} the new pass to \code{$history}.
#' \code{print()} displays the full audit trail across all passes.
#'
#' @param geno A numeric matrix (\code{lines x markers}) with row names
#'   identifying lines, a \code{data.frame} with a line identifier column, a
#'   \code{"checkPanel"} object from \code{\link{checkPanel}}, or a
#'   \code{"filteredPanel"} object from a previous \code{filterPanel()} call.
#' @param map A \code{data.frame} containing the genetic map.  Not required
#'   when \code{geno} is a \code{"checkPanel"} or \code{"filteredPanel"} object
#'   (the stored map is used automatically).
#' @param id Character string naming the line identifier column when
#'   \code{geno} is a \code{data.frame}.  Default \code{"id"}.
#' @param map.id Character string naming the marker column in \code{map}.
#'   Default \code{"marker"}.
#' @param map.chr Character string naming the chromosome column in \code{map}.
#'   Default \code{"chr"}.
#' @param map.pos Character string naming the marker position column in
#'   \code{map}.  Positions are used as supplied — no unit conversion is
#'   applied.  Default \code{"pos"}.
#' @param encoding Character string specifying the genotype encoding:
#'   \code{"012"} (default) or \code{"pm1"}.
#' @param steps A named list controlling which filters are applied and in what
#'   order.  Each element name must be one of \code{"map"},
#'   \code{"miss.marker"}, \code{"miss.line"}, \code{"het.line"},
#'   \code{"het.marker"}, \code{"dup.lines"}, \code{"dup.markers"},
#'   \code{"maf"}.  Numeric steps are skipped when their value is
#'   \code{NULL}; logical steps (\code{map}, \code{dup.lines},
#'   \code{dup.markers}) are skipped when \code{FALSE}.  Map consistency is
#'   always prepended as step 1 if not already first in the list.  The
#'   default list reproduces the original fixed workflow.
#' @param \dots Currently ignored; reserved for future use.
#'
#' @return An object of class \code{"filteredPanel"} — a list containing:
#' \describe{
#'   \item{\code{$geno}}{Filtered genotype matrix (current state).}
#'   \item{\code{$map}}{Filtered map data frame (current state).}
#'   \item{\code{$encoding}, \code{$id}, \code{$map.id}, \code{$map.chr},
#'     \code{$map.pos}}{Carried through for \code{\link{primePanel}}.}
#'   \item{\code{$history}}{A list of pass records, one per \code{filterPanel}
#'     call.  Each pass contains \code{$pass} (integer), \code{$steps} (the
#'     steps list used), \code{$removed} (named list of removed items per
#'     step), \code{$n.before}, and \code{$n.after} (named integer vectors).}
#'   \item{\code{$n.original}}{Named integer vector (\code{lines},
#'     \code{markers}) recording the dimensions \emph{before the very first
#'     filtering pass}.}
#'   \item{\code{$n.final}}{Named integer vector (\code{lines},
#'     \code{markers}) recording the current dimensions after all passes.}
#' }
#'
#' @seealso \code{\link{checkPanel}}, \code{\link{primePanel}}
#' @export
filterPanel <- function(geno, ...) UseMethod("filterPanel")


# =============================================================================
# filterPanel.default  —  raw matrix or data.frame
# =============================================================================

#' @rdname filterPanel
#' @export
filterPanel.default <- function(geno,
                                map,
                                id       = "id",
                                map.id   = "marker",
                                map.chr  = "chr",
                                map.pos  = "pos",
                                encoding = "012",
                                steps    = list(
                                    map         = TRUE,
                                    miss.marker = 0.20,
                                    miss.line   = 0.20,
                                    het.line    = NULL,
                                    het.marker  = NULL,
                                    dup.lines   = TRUE,
                                    dup.markers = FALSE,
                                    maf         = 0.05
                                ), ...) {

    geno.mat <- .extract_geno_matrix(geno, id)
    n.orig   <- c(lines = nrow(geno.mat), markers = ncol(geno.mat))

    result   <- .run_steps(geno.mat, map, encoding, id, map.id, map.chr,
                           map.pos, steps)

    .build_filteredPanel(
        result     = result,
        encoding   = encoding,
        id         = id,
        map.id     = map.id,
        map.chr    = map.chr,
        map.pos    = map.pos,
        steps      = steps,
        prev_hist  = list(),
        n.original = n.orig,
        n.before   = n.orig
    )
}


# =============================================================================
# filterPanel.checkPanel  —  accepts a checkPanel object
# =============================================================================

#' @rdname filterPanel
#' @export
filterPanel.checkPanel <- function(geno,
                                   steps = list(
                                       map         = TRUE,
                                       miss.marker = 0.20,
                                       miss.line   = 0.20,
                                       het.line    = NULL,
                                       het.marker  = NULL,
                                       dup.lines   = TRUE,
                                       dup.markers = FALSE,
                                       maf         = 0.05
                                   ), ...) {
    chk <- geno
    filterPanel.default(
        geno     = chk$geno,
        map      = chk$map,
        id       = chk$id,
        map.id   = chk$map.id,
        map.chr  = chk$map.chr,
        map.pos  = chk$map.pos,
        encoding = chk$encoding,
        steps    = steps
    )
}


# =============================================================================
# filterPanel.filteredPanel  —  additional pass; appends to $history
# =============================================================================

#' @rdname filterPanel
#' @export
filterPanel.filteredPanel <- function(geno,
                                      steps = list(
                                          map         = TRUE,
                                          miss.marker = 0.20,
                                          miss.line   = 0.20,
                                          het.line    = NULL,
                                          het.marker  = NULL,
                                          dup.lines   = TRUE,
                                          dup.markers = FALSE,
                                          maf         = 0.05
                                      ), ...) {
    prev     <- geno
    n.before <- prev$n.final

    result   <- .run_steps(prev$geno, prev$map, prev$encoding,
                           prev$id, prev$map.id, prev$map.chr,
                           prev$map.pos, steps)

    .build_filteredPanel(
        result     = result,
        encoding   = prev$encoding,
        id         = prev$id,
        map.id     = prev$map.id,
        map.chr    = prev$map.chr,
        map.pos    = prev$map.pos,
        steps      = steps,
        prev_hist  = prev$history,
        n.original = prev$n.original,
        n.before   = n.before
    )
}


# =============================================================================
# .run_steps()  —  workhorse: iterate over the steps list
# =============================================================================

.run_steps <- function(geno.mat, map, encoding, id, map.id, map.chr,
                       map.pos, steps) {

    encoding <- match.arg(encoding, c("012", "pm1"))

    # --- Validate step names --------------------------------------------------
    bad <- setdiff(names(steps), .valid_step_names)
    if (length(bad))
        stop("Unknown step name(s) in 'steps': ",
             paste(bad, collapse = ", "),
             ".  Valid names: ",
             paste(.valid_step_names, collapse = ", "), ".")

    # --- Ensure map consistency is always the first step ---------------------
    if (!identical(names(steps)[1L], "map")) {
        steps <- c(list(map = TRUE), steps[names(steps) != "map"])
    }

    for (col in c(map.id, map.chr, map.pos))
        if (!col %in% names(map))
            stop("Column '", col, "' not found in map.")

    # Helper: is a genotype call heterozygous under the current encoding?
    .is_het <- function(mat)
        if (encoding == "012") mat == 1L else mat == 0L

    removed     <- list()
    steps_run   <- list()   # records what actually ran (including forced map)

    for (nm in names(steps)) {
        val <- steps[[nm]]

        # ---- map consistency ------------------------------------------------
        if (nm == "map") {
            steps_run[[nm]] <- val
            if (isFALSE(val)) { removed[[nm]] <- character(0L); next }
            map.markers  <- as.character(map[[map.id]])
            geno.markers <- colnames(geno.mat)
            drop_m       <- setdiff(geno.markers, map.markers)
            if (length(drop_m)) {
                geno.mat <- geno.mat[, !colnames(geno.mat) %in% drop_m,
                                     drop = FALSE]
                map      <- map[map[[map.id]] %in% colnames(geno.mat), ,
                                drop = FALSE]
            }
            removed[["map_consistency"]] <- drop_m
            next
        }

        # ---- marker missingness ---------------------------------------------
        if (nm == "miss.marker") {
            steps_run[[nm]] <- val
            if (is.null(val)) { removed[[nm]] <- character(0L); next }
            mm     <- colMeans(is.na(geno.mat))
            drop_m <- names(mm)[mm > val]
            if (length(drop_m)) {
                geno.mat <- geno.mat[, !colnames(geno.mat) %in% drop_m,
                                     drop = FALSE]
                map      <- map[map[[map.id]] %in% colnames(geno.mat), ,
                                drop = FALSE]
            }
            removed[[nm]] <- drop_m
            next
        }

        # ---- line missingness -----------------------------------------------
        if (nm == "miss.line") {
            steps_run[[nm]] <- val
            if (is.null(val)) { removed[[nm]] <- character(0L); next }
            ml     <- rowMeans(is.na(geno.mat))
            drop_l <- names(ml)[ml > val]
            if (length(drop_l))
                geno.mat <- geno.mat[!rownames(geno.mat) %in% drop_l, ,
                                     drop = FALSE]
            removed[[nm]] <- drop_l
            next
        }

        # ---- line heterozygosity --------------------------------------------
        if (nm == "het.line") {
            steps_run[[nm]] <- val
            if (is.null(val)) { removed[[nm]] <- character(0L); next }
            hl     <- rowMeans(.is_het(geno.mat), na.rm = TRUE)
            drop_l <- names(hl)[!is.na(hl) & hl > val]
            if (length(drop_l))
                geno.mat <- geno.mat[!rownames(geno.mat) %in% drop_l, ,
                                     drop = FALSE]
            removed[[nm]] <- drop_l
            next
        }

        # ---- marker heterozygosity ------------------------------------------
        if (nm == "het.marker") {
            steps_run[[nm]] <- val
            if (is.null(val)) { removed[[nm]] <- character(0L); next }
            hm     <- colMeans(.is_het(geno.mat), na.rm = TRUE)
            drop_m <- names(hm)[!is.na(hm) & hm > val]
            if (length(drop_m)) {
                geno.mat <- geno.mat[, !colnames(geno.mat) %in% drop_m,
                                     drop = FALSE]
                map      <- map[map[[map.id]] %in% colnames(geno.mat), ,
                                drop = FALSE]
            }
            removed[[nm]] <- drop_m
            next
        }

        # ---- duplicate lines ------------------------------------------------
        if (nm == "dup.lines") {
            steps_run[[nm]] <- val
            if (isFALSE(val)) { removed[[nm]] <- character(0L); next }
            line_sigs <- apply(geno.mat, 1L, function(x)
                paste(ifelse(is.na(x), "NA", as.character(x)), collapse = "|"))
            dup_idx <- duplicated(line_sigs)
            dup_nms <- rownames(geno.mat)[dup_idx]
            if (length(dup_nms))
                geno.mat <- geno.mat[!dup_idx, , drop = FALSE]
            removed[[nm]] <- dup_nms
            next
        }

        # ---- duplicate markers ----------------------------------------------
        if (nm == "dup.markers") {
            steps_run[[nm]] <- val
            if (isFALSE(val)) { removed[[nm]] <- character(0L); next }
            mark_sigs <- apply(geno.mat, 2L, function(x)
                paste(ifelse(is.na(x), "NA", as.character(x)), collapse = "|"))
            dup_midx <- duplicated(mark_sigs)
            dup_mnms <- colnames(geno.mat)[dup_midx]
            if (length(dup_mnms)) {
                geno.mat <- geno.mat[, !dup_midx, drop = FALSE]
                map      <- map[map[[map.id]] %in% colnames(geno.mat), ,
                                drop = FALSE]
            }
            removed[[nm]] <- dup_mnms
            next
        }

        # ---- MAF ------------------------------------------------------------
        if (nm == "maf") {
            steps_run[[nm]] <- val
            if (is.null(val) || val <= 0) { removed[[nm]] <- character(0L); next }
            col_means <- colMeans(geno.mat, na.rm = TRUE)
            p_vec     <- if (encoding == "012") col_means / 2
                         else (col_means + 1) / 2
            maf_vec   <- pmin(p_vec, 1 - p_vec)
            drop_maf  <- names(maf_vec)[!is.na(maf_vec) & maf_vec < val]
            if (length(drop_maf)) {
                geno.mat <- geno.mat[, !colnames(geno.mat) %in% drop_maf,
                                     drop = FALSE]
                map      <- map[map[[map.id]] %in% colnames(geno.mat), ,
                                drop = FALSE]
            }
            removed[[nm]] <- drop_maf
            next
        }
    }

    if (nrow(geno.mat) == 0L)
        warning("All lines were removed by the filters. ",
                "Consider relaxing one or more thresholds.")
    if (ncol(geno.mat) == 0L)
        warning("All markers were removed by the filters. ",
                "Consider relaxing one or more thresholds.")

    list(geno.mat  = geno.mat,
         map       = map,
         removed   = removed,
         steps_run = steps_run)
}


# =============================================================================
# .build_filteredPanel()  —  assemble the return object
# =============================================================================

.build_filteredPanel <- function(result, encoding, id, map.id, map.chr,
                                 map.pos, steps, prev_hist,
                                 n.original, n.before) {

    n.after <- c(lines   = nrow(result$geno.mat),
                 markers = ncol(result$geno.mat))

    new_pass <- list(
        pass     = length(prev_hist) + 1L,
        steps    = result$steps_run,
        removed  = result$removed,
        n.before = n.before,
        n.after  = n.after
    )

    structure(
        list(
            geno       = result$geno.mat,
            map        = result$map,
            encoding   = encoding,
            id         = id,
            map.id     = map.id,
            map.chr    = map.chr,
            map.pos    = map.pos,
            history    = c(prev_hist, list(new_pass)),
            n.original = n.original,
            n.final    = n.after
        ),
        class = "filteredPanel"
    )
}


# =============================================================================
# .extract_geno_matrix()  —  shared matrix extraction
# =============================================================================

.extract_geno_matrix <- function(geno, id) {
    if (is.data.frame(geno)) {
        if (!id %in% names(geno))
            stop("Column '", id, "' not found in geno data frame.")
        ids      <- as.character(geno[[id]])
        geno.mat <- as.matrix(geno[, !names(geno) %in% id, drop = FALSE])
        rownames(geno.mat) <- ids
    } else {
        geno.mat <- as.matrix(geno)
        if (is.null(rownames(geno.mat)))
            stop("geno matrix must have row names identifying lines.")
    }
    storage.mode(geno.mat) <- "numeric"
    geno.mat
}


# =============================================================================
# .step_label()  —  human-readable label for a step
# =============================================================================

.step_label <- function(nm, val) {
    switch(nm,
        map         = "Map consistency",
        miss.marker = sprintf("Marker missingness > %.0f%%", 100 * val),
        miss.line   = sprintf("Line missingness > %.0f%%",   100 * val),
        het.line    = sprintf("Line heterozygosity > %.0f%%",   100 * val),
        het.marker  = sprintf("Marker heterozygosity > %.0f%%", 100 * val),
        dup.lines   = "Duplicate lines",
        dup.markers = "Duplicate markers",
        maf         = sprintf("MAF < %.2f", val),
        nm   # fallback: use name as-is
    )
}

.step_what <- function(nm) {
    switch(nm,
        map         = "marker",
        miss.marker = "marker",
        miss.line   = "line",
        het.line    = "line",
        het.marker  = "marker",
        dup.lines   = "line",
        dup.markers = "marker",
        maf         = "marker",
        "item"
    )
}


# =============================================================================
# print.filteredPanel
# =============================================================================

#' Print a filteredPanel summary
#'
#' Displays the full filtering history across all passes, one section per
#' \code{filterPanel()} call, followed by an overall summary line.
#'
#' @param x A \code{"filteredPanel"} object from \code{\link{filterPanel}}.
#' @param \dots Currently ignored.
#' @return \code{x} invisibly.
#' @exportS3Method
print.filteredPanel <- function(x, ...) {

    n_passes <- length(x$history)
    cat("\nFilter History",
        if (n_passes > 1L) sprintf("  [%d passes]", n_passes) else "",
        "\n", sep = "")
    cat(strrep("=", 52), "\n", sep = "")
    cat(sprintf("  Original : %d lines,  %d markers\n",
                x$n.original["lines"], x$n.original["markers"]))

    for (h in x$history) {
        cat("\n")
        if (n_passes > 1L)
            cat(sprintf("  -- Pass %d ", h$pass),
                strrep("-", max(0L, 40L - nchar(h$pass))), "\n", sep = "")

        step_num <- 0L
        # map_consistency is stored under "map_consistency" key; all others
        # use their step name as the removed key.
        removed_keys <- names(h$removed)

        for (nm in names(h$steps)) {
            step_num  <- step_num + 1L
            val       <- h$steps[[nm]]

            # Determine the removed-list key for this step name
            rem_key <- if (nm == "map") "map_consistency" else nm
            rem     <- h$removed[[rem_key]]

            # Is this step active?
            active <- if (nm %in% c("map", "dup.lines", "dup.markers"))
                          isTRUE(val)
                      else
                          !is.null(val)

            if (!active) {
                cat(sprintf("  Step %d [%-35s]: skipped\n",
                            step_num, .step_label(nm, val)))
            } else {
                cat(sprintf("  Step %d [%-35s]: removed %d %s(s)\n",
                            step_num,
                            .step_label(nm, val),
                            length(rem),
                            .step_what(nm)))
            }
        }

        cat(sprintf("  After    : %d lines,  %d markers\n",
                    h$n.after["lines"], h$n.after["markers"]))
    }

    cat("\n")
    if (n_passes > 1L)
        cat(sprintf("  Final    : %d lines,  %d markers  (from %d lines, %d markers)\n",
                    x$n.final["lines"],    x$n.final["markers"],
                    x$n.original["lines"], x$n.original["markers"]))
    else
        cat(sprintf("  Remaining: %d lines,  %d markers\n",
                    x$n.final["lines"], x$n.final["markers"]))
    cat("\n")
    invisible(x)
}

# Null-coalescing helper (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b
