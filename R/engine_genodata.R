# =============================================================================
# engine_genodata.R
# Internal functions for genotype data construction and line-fixing.
# =============================================================================

#' @keywords internal
.buildGenoData <- function(intervalObj, gen.type, glines, plines) {
    if (gen.type == "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$interval.data)
    else
        gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    nint <- lapply(gdat, function(el) 1:ncol(el))
    lint <- unlist(lapply(nint, length))
    mnams <- paste("Chr", rep(names(intervalObj$geno), times = lint), unlist(nint), sep = ".")
    dimnames(genoData) <- list(as.character(glines), mnams)
    genoData <- genoData[rownames(genoData) %in% as.character(plines), , drop = FALSE]
    state <- rep(1, ncol(genoData))
    names(state) <- mnams
    list(genoData = genoData, mnams = mnams, state = state)
}

#' @keywords internal
.fixLines <- function(baseModel, phenoData, genoData, merge.by, plines, fix.lines, ...) {
    rterms <- unlist(strsplit(deparse(baseModel$call$random[[2]]), " \\+ "))
    rterms <- rterms[-grep(merge.by, rterms)]
    whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)
    genetic.term <- merge.by
    if (!all(whg) & fix.lines) {
        phenoData$Gomit <- phenoData$Gsave <- plines
        levels(phenoData$Gsave)[!whg] <- NA
        levels(phenoData$Gomit)[whg] <- "GEN"
        fix.form <- as.formula(". ~ Gomit + .")
        ran.base <- formula(paste(c("~ Gsave", rterms), collapse = " + "))
        baseModel$call$data <- quote(phenoData)
        cat("\nFixing lines and updating initial base model:\n")
        cat("============================================\n")
        baseModel <- update(baseModel, fixed. = fix.form, random. = ran.base, ...)
        merge.by <- "Gsave"
    }
    list(baseModel = baseModel, phenoData = phenoData, merge.by = merge.by,
         rterms = rterms, genetic.term = genetic.term)
}
