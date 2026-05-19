# =============================================================================
# engine_effects.R
# Internal functions for merging selected effects into phenoData and
# adding them as fixed or random terms to the ASReml models.
# =============================================================================

#' @keywords internal
.mergeEffect <- function(phenoData, genoData, qtl, merge.by) {
    tmp <- cbind.data.frame(rownames(genoData), genoData[, qtl])
    qtl.x <- gsub("Chr\\.", "X.", qtl)
    names(tmp) <- c(merge.by, qtl.x)
    phenoData <- cbind.data.frame(ord = 1:nrow(phenoData), phenoData)
    phenoData <- merge(phenoData, tmp, by = merge.by, all.x = TRUE, all.y = FALSE)
    phenoData <- phenoData[order(phenoData$ord), ]
    phenoData <- phenoData[, -2]
    list(phenoData = phenoData, qtl.x = qtl.x)
}

#' @keywords internal
#' NOTE: phenoData must be a named parameter here so that ASReml's update()
#' can resolve quote(phenoData) in the model call from this function's frame.
#'
#' Trait : character name of the Trait factor (NULL = univariate).
#'   NULL  -> adds qtl.x only (current behaviour, unchanged).
#'   non-NULL -> adds qtl.x AND Trait:qtl.x to capture both main and
#'               interaction effects for each detected QTL.
.addEffect <- function(baseModel, qtlModel, phenoData, merge.by, qtl.x,
                       method, iter, Trait = NULL, ...) {
    # Build the fixed/random term string.
    # Univariate: just qtl.x
    # Multivariate: qtl.x + Trait:qtl.x (main + interaction)
    fix.terms <- if (is.null(Trait))
        qtl.x
    else
        paste(qtl.x, paste(Trait, qtl.x, sep = ":"), sep = " + ")

    if (method == "random") {
        ran.form <- formula(paste("~ . +", qtl.x))   # random path: main only (no interaction)
        cat("\nRandom Effects QTL Model Iteration (", iter, "):\n")
        cat("========================================\n")
        baseModel <- .vModify(baseModel, merge.by)
        baseModel <- update(baseModel, random. = ran.form, ...)
        cat("\nRandom Effects QTL plus Interval/Marker Model Iteration (", iter, "):\n")
        cat("=============================================================\n")
        qtlModel <- .vModify(qtlModel, merge.by)
        qtlModel <- update(qtlModel, random. = ran.form, ...)
        list.coefs <- qtlModel$coefficients$random
        zind   <- grep("X\\.", rownames(list.coefs))
        coefs  <- list.coefs[zind, 1]
        names(coefs) <- rownames(list.coefs)[zind]
        vcoefs <- qtlModel$vcoeff$random[zind]
    } else {
        fix.form <- as.formula(paste(". ~ . +", fix.terms))
        cat("\nFixed Effects QTL Model Iteration (", iter, "):\n")
        cat("========================================\n")
        baseModel <- update(baseModel, fixed. = fix.form, ...)
        cat("\nFixed Effects QTL plus Interval/Marker Model Iteration (", iter, "):\n")
        cat("============================================================\n")
        qtlModel <- update(qtlModel, fixed. = fix.form, ...)
        list.coefs <- qtlModel$coefficients$fixed
        zind   <- grep("X\\.", rownames(list.coefs))
        coefs  <- rev(list.coefs[zind, 1])
        names(coefs) <- rev(rownames(list.coefs)[zind])
        vcoefs <- rev(qtlModel$vcoeff$fixed[zind])
    }
    list(baseModel = baseModel, qtlModel = qtlModel, coefs = coefs, vcoefs = vcoefs)
}
