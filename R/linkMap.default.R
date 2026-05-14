# =============================================================================
# linkMap.default.R
# Backward-compatible dispatcher for list-based multi-model calls.
#
# Handles:
#   linkMap(list(qtl_yld, qtl_tgw), genObj = x)  # old qtlAim API
#   linkMap(list(gwas1, gwas2),      genObj = x)  # old gwasAim API
#
# Unwraps the list and delegates to the appropriate typed method, which
# detects the extra models via its own ... inspection path.
#
# The preferred new API passes models directly:
#   linkMap(qtl_yld, qtl_tgw, genObj = x)
#   linkMap(gwas1,   gwas2,   genObj = x)
# =============================================================================

#' @rdname linkMap
#' @export
linkMap.default <- function(object, ...) {
    if (!is.list(object))
        stop("No linkMap method for class '",
             paste(class(object), collapse = ", "), "'.")
    if (length(object) == 0L)
        stop("'object' is an empty list.")
    types <- unique(vapply(object, function(x) class(x)[1L], character(1L)))
    if (length(types) != 1L)
        stop("All elements of 'object' must have the same class.")
    fn <- switch(types,
        qtlAim  = linkMap.qtlAim,
        gwasAim = linkMap.gwasAim,
        stop("No multi-model linkMap method for class '", types, "'.")
    )
    ## Because ... precedes genObj in the typed method signatures,
    ## those params must be named in every call.  Positional extras therefore
    ## land safely in ..., never accidentally matching a named parameter.
    extras <- unname(object[-1L])
    do.call(fn, c(list(object[[1L]]), extras, list(...)))
}
