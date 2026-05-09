# =============================================================================
# engine_validate.R
# Internal shared validation for wgAim analyses.
# =============================================================================

#' @keywords internal
.validateModel <- function(baseModel, merge.by, method, selection, breakout) {
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
        if (!baseModel$converge)
            stop("Base model not converged: Check base model before proceeding.")
    }
    asremlEnv <- lapply(baseModel$formulae, function(el) attr(el, ".Environment"))
    phenoData <- eval(baseModel$call$data)
    if (is.null(merge.by))
        stop("merge.by: name of the column linking phenotypic and genotypic data is required.")
    if (!(method %in% c("fixed", "random")))
        stop("method must be either \"fixed\" or \"random\".")
    if (!(selection %in% c("interval", "chromosome")))
        stop("selection must be either \"interval\" or \"chromosome\".")
    if (!is.numeric(breakout) | breakout < -1 | breakout == 0)
        stop("breakout must be -1 or a positive integer.")
    list(baseModel = baseModel, asremlEnv = asremlEnv, phenoData = phenoData)
}
