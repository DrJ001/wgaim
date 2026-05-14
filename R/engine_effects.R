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
.addEffect <- function(baseModel, qtlModel, phenoData, merge.by, qtl.x, method, iter, ...) {
    if (method == "random") {
        ran.form <- formula(paste("~ . +", qtl.x))
        cat("\nRandom Effects QTL Model Iteration (", iter, "):\n")
        cat("========================================\n")
        baseModel <- .vModify(baseModel, merge.by)
        baseModel <- update(baseModel, random. = ran.form, ...)
        cat("\nRandom Effects QTL plus Interval/Marker Model Iteration (", iter, "):\n")
        cat("=============================================================\n")
        qtlModel <- .vModify(qtlModel, merge.by)
        qtlModel <- update(qtlModel, random. = ran.form, ...)
        list.coefs <- qtlModel$coefficients$random
        zind <- grep("X\\.", rownames(list.coefs))
        coefs <- list.coefs[zind, 1]
        names(coefs) <- rownames(list.coefs)[zind]
        vcoefs <- qtlModel$vcoeff$random[zind]
    } else {
        fix.form <- as.formula(paste(". ~ . +", qtl.x))
        cat("\nFixed Effects QTL Model Iteration (", iter, "):\n")
        cat("========================================\n")
        baseModel <- update(baseModel, fixed. = fix.form, ...)
        cat("\nFixed Effects QTL plus Interval/Marker Model Iteration (", iter, "):\n")
        cat("============================================================\n")
        qtlModel <- update(qtlModel, fixed. = fix.form, ...)
        list.coefs <- qtlModel$coefficients$fixed
        zind <- grep("X\\.", rownames(list.coefs))
        coefs <- rev(list.coefs[zind, 1])
        names(coefs) <- rev(rownames(list.coefs)[zind])
        vcoefs <- rev(qtlModel$vcoeff$fixed[zind])
    }
    list(baseModel = baseModel, qtlModel = qtlModel, coefs = coefs, vcoefs = vcoefs)
}
