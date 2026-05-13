# =============================================================================
# aimTable.R
# Stack QTL/GWAS summary tables from multiple fitted models into one
# combined 'super-table' data.frame ready for export.
# =============================================================================

#' Stack QTL summary tables from multiple fitted models
#'
#' Combines the \code{summary()} output from two or more \code{qtlAim} or
#' \code{gwasAim} objects into a single stacked \code{data.frame} with a
#' leading \code{Trait} column, ready for export to a spreadsheet or a
#' LaTeX table via \pkg{xtable}.
#'
#' All objects supplied via \code{...} must satisfy two rules: (1) they must
#' all be of the same class (all \code{qtlAim} \emph{or} all \code{gwasAim};
#' mixing is not permitted); and (2) for \code{qtlAim} objects, all models
#' must have been run with the same \code{gen.type} (all \code{"interval"}
#' \emph{or} all \code{"marker"}), ensuring a consistent column layout in the
#' stacked table.  Objects for which no QTL / significant markers were
#' detected are silently dropped from the result with a console message.
#'
#' @param \dots One or more fitted objects of class \code{"qtlAim"} or
#'   \code{"gwasAim"}.
#' @param genObj The genomic data object required by \code{summary()}: a
#'   \code{"wgCross"} produced by \code{\link{primeCross}} (for
#'   \code{qtlAim} objects) or a \code{"wgPanel"} produced by
#'   \code{\link{primePanel}} (for \code{gwasAim} objects).  Supply a
#'   single shared object when all models used the same genomic data, or a
#'   \code{list} of objects with one element per model in \code{...}.
#' @param labels Optional character vector of trait labels, one per model
#'   in \code{...}.  If \code{NULL} (default), the left-hand side of the
#'   \code{fixed} component of each model's call is used.
#' @param columns Either \code{"all"} (default) to include every column
#'   returned by \code{summary()}, or a numeric vector of column indices to
#'   retain (applied \emph{after} the leading \code{Trait} column is
#'   inserted, so column 1 in \code{summary()} output corresponds to
#'   \code{columns = 1}).
#' @param LOD Logical; passed through to \code{summary.qtlAim()} /
#'   \code{summary.gwasAim()}.  If \code{TRUE} (default), a LOD score
#'   column is appended.
#'
#' @return A \code{data.frame} of class \code{c("aimTable", "data.frame")}
#'   with a leading \code{Trait} column followed by the columns from the
#'   individual \code{summary()} tables.  Rows with no detected QTL are
#'   excluded.  The attribute \code{"obj.class"} records the common object
#'   class (\code{"qtlAim"} or \code{"gwasAim"}).
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{gwasAim}},
#'   \code{\link{summary.qtlAim}}, \code{\link{summary.gwasAim}}
#'
#' @examples
#' \dontrun{
#' ## Two QTL models sharing the same cross object
#' gw.qtl  <- qtlAim(gw.asf,  genObj = genoObj, merge.by = "Genotype")
#' yld.qtl <- qtlAim(yld.asf, genObj = genoObj, merge.by = "Genotype")
#'
#' at <- aimTable(gw.qtl, yld.qtl,
#'                genObj = genoObj,
#'                labels = c("Grain weight", "Yield"))
#' print(at)
#'
#' ## Export to LaTeX (requires xtable)
#' # print(xtable::xtable(at), include.rownames = FALSE)
#'
#' ## Two GWAS models
#' gw.gwas  <- gwasAim(gw.asf,  genObj = panelObj, merge.by = "Genotype")
#' yld.gwas <- gwasAim(yld.asf, genObj = panelObj, merge.by = "Genotype")
#'
#' at2 <- aimTable(gw.gwas, yld.gwas, genObj = panelObj)
#' print(at2)
#' }
#'
#' @export
aimTable <- function(..., genObj, labels = NULL, columns = "all", LOD = TRUE) {

    dots <- list(...)

    # ---- basic argument checks -----------------------------------------------

    if (length(dots) == 0L)
        stop("Supply at least one 'qtlAim' or 'gwasAim' object via '...'.")

    if (missing(genObj))
        stop("'genObj' is a required argument.")

    # Determine and enforce a single common object class
    obj_classes <- vapply(dots, function(el) {
        if (inherits(el, "qtlAim"))  return("qtlAim")
        if (inherits(el, "gwasAim")) return("gwasAim")
        stop("All objects in '...' must be of class 'qtlAim' or 'gwasAim'.")
    }, character(1L))

    if (length(unique(obj_classes)) > 1L)
        stop("All objects in '...' must share the same class ",
             "(all 'qtlAim' or all 'gwasAim'; mixing is not permitted).")

    obj_class <- obj_classes[1L]

    # For qtlAim objects, enforce a common gen.type (interval vs marker) so
    # that all summary() calls return the same column layout.
    if (obj_class == "qtlAim") {
        gen_types <- vapply(dots, function(el) el$QTL$type, character(1L))
        if (length(unique(gen_types)) > 1L)
            stop("All 'qtlAim' objects must have been analysed with the same ",
                 "gen.type. Found: ",
                 paste(unique(gen_types), collapse = " and "), ".")
    }

    # ---- genObj — allow a single shared object or a per-model list -----------

    if (!is.list(genObj) ||
        inherits(genObj, "wgCross") ||
        inherits(genObj, "wgPanel")) {
        # Single shared genomic object — replicate for every model
        genObj_list <- rep(list(genObj), length(dots))
    } else {
        genObj_list <- genObj
        if (length(genObj_list) != length(dots))
            stop("When 'genObj' is a list it must have one element per model ",
                 "supplied in '...'.")
    }

    # Validate individual genObj types against the object class
    expected_geno <- if (obj_class == "qtlAim") "wgCross" else "wgPanel"
    invisible(lapply(genObj_list, function(g) {
        if (!inherits(g, expected_geno))
            stop("'genObj' must be of class '", expected_geno, "' for '",
                 obj_class, "' objects.")
    }))

    # ---- labels --------------------------------------------------------------

    if (!is.null(labels)) {
        if (length(labels) != length(dots))
            stop("'labels' must have the same length as the number of models ",
                 "supplied in '...'.")
    } else {
        labels <- vapply(dots, function(el) {
            tryCatch(
                deparse(el$call$fixed[[2L]]),
                error = function(e) "unknown"
            )
        }, character(1L))
    }

    # ---- build individual summary tables -------------------------------------

    has_qtl <- vapply(dots, function(el) !is.null(el$QTL$effects), logical(1L))

    if (!any(has_qtl)) {
        message("No models have detected QTL / significant markers. ",
                "Returning an empty table.")
        empty <- data.frame(Trait = character(0L))
        class(empty) <- c("aimTable", "data.frame")
        attr(empty, "obj.class") <- obj_class
        return(invisible(empty))
    }

    if (any(!has_qtl)) {
        dropped <- labels[!has_qtl]
        message("The following trait(s) had no detected QTL and were excluded: ",
                paste(dropped, collapse = ", "))
    }

    active_dots   <- dots[has_qtl]
    active_genObj <- genObj_list[has_qtl]
    active_labels <- labels[has_qtl]

    summ_list <- mapply(
        function(obj, gobj) summary(obj, gobj, LOD = LOD),
        active_dots, active_genObj,
        SIMPLIFY = FALSE
    )

    # ---- column subsetting ---------------------------------------------------

    if (!identical(columns, "all")) {
        if (!is.numeric(columns))
            stop("'columns' must be either \"all\" or a numeric vector of ",
                 "column indices.")
        summ_list <- lapply(summ_list, function(s) {
            idx <- columns[columns >= 1L & columns <= ncol(s)]
            if (length(idx) == 0L)
                stop("'columns' contains no valid column indices for the ",
                     "summary table (which has ", ncol(s), " columns).")
            s[, idx, drop = FALSE]
        })
    }

    # ---- stack into super-table ----------------------------------------------

    nrows  <- vapply(summ_list, nrow, integer(1L))
    qtab   <- do.call(rbind.data.frame, summ_list)
    qtab   <- cbind.data.frame(Trait = rep(active_labels, times = nrows),
                                qtab,
                                stringsAsFactors = FALSE)
    rownames(qtab) <- NULL

    class(qtab)          <- c("aimTable", "data.frame")
    attr(qtab, "obj.class") <- obj_class
    qtab
}
