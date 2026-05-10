# =============================================================================
# summary.qtlAim.R
# S3 summary method for qtlAim objects.
# =============================================================================

#' @describeIn qtlAim Produce a detailed summary table of detected QTL, sorted
#'   by chromosome and position. Returns a \code{data.frame} with columns for
#'   chromosome, flanking markers and their cM positions, inferred interval
#'   midpoint, QTL effect size, p-value or posterior probability, percentage of
#'   phenotypic variance explained, and (optionally) LOD score. Returns
#'   \code{NULL} invisibly if no QTL were detected.
#' @param object A \code{qtlAim} object.
#' @param intervalObj The \code{"interval"} or \code{"cross"} object used in
#'   the analysis, required to resolve QTL positions.
#' @param LOD Logical. If \code{TRUE} (default), a LOD score column is
#'   appended to the summary table. LOD is computed from the Wald z-statistic
#'   as \eqn{0.5 \log_{10}(\exp(z^2))}.
#' @export
summary.qtlAim <- function(object, intervalObj, LOD = TRUE, ...) {
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj must be of class \"cross\".")
    if (is.null(qtle <- object$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
        return()
    }
    sigma2 <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1
    if (object$QTL$type == "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$interval.data)
    else
        gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    gterm <- object$QTL$diag$genetic.term
    scale <- object$QTL$diag$rel.scale
    dimnames(genoData) <- list(as.character(intervalObj$pheno[[gterm]]),
                               names(object$QTL$diag$state))
    genoSub <- genoData[, as.logical(object$QTL$diag$state), drop = FALSE]
    if ("Gsave" %in% names(object$mf))
        gterm <- "Gsave"
    genoSub <- genoSub[rownames(genoSub) %in% levels(object$mf[[gterm]]), , drop = FALSE]
    coef.mark <- mean(apply(genoSub, 1, function(el) sum(el * el)), na.rm = TRUE)
    mark.terms <- paste("mbf*.*ints*", paste0("vm\\(", gterm, "*"), sep = "|")
    oth.terms  <- object$vparameters[-grep(mark.terms, names(object$vparameters))]
    var.mark   <- sigma2 * object$vparameters[grep(mark.terms, names(object$vparameters))] / scale
    var.res    <- sigma2 * oth.terms[grep(gterm, names(oth.terms))]
    if (object$QTL$method == "random") {
        var.est  <- sigma2 * object$vparameters[grep("X\\.", names(object$vparameters))]
        coef.est <- apply(genoData[, object$QTL$qtl, drop = FALSE]^2, 2, mean, na.rm = TRUE)
    } else {
        var.est  <- qtle^2
        coef.est <- rep(1, length(qtle))
    }
    var.all  <- sum(c(coef.est, coef.mark, 1) * c(var.est, var.mark, var.res))
    perc.var <- round(100 * (coef.est * var.est) / var.all, 1)
    zrat <- qtle / sqrt(object$QTL$veffects * sigma2)
    if (object$QTL$method == "random") {
        pvalue <- round((1 - pchisq(zrat^2, 1)) / 2, 4)
        pvalue[pvalue < 0.0001] <- "<0.0001"
    } else {
        pvalue <- round(2 * (1 - pnorm(abs(zrat))), 4)
        pvalue[as.numeric(pvalue) < 0.0001] <- "<0.0001"
    }
    lod <- round(0.5 * log10(exp(zrat^2)), 2)
    qtlm <- getQTL(object, intervalObj)
    if (object$QTL$type == "interval") {
        qtab <- data.frame(
            Chromosome    = qtlm[, 1],
            "Left Marker" = qtlm[, 3],
            "dist(cM)"    = as.numeric(qtlm[, 4]),
            "Infer. Marker" = qtlm[, 5],
            "dist(cM)"    = as.numeric(qtlm[, 6]),
            "Right Marker" = qtlm[, 7],
            "dist(cM)"    = as.numeric(qtlm[, 8]),
            Size          = round(qtle, 4),
            Prob          = pvalue,
            "Perc.Var"    = perc.var,
            check.names   = FALSE,
            stringsAsFactors = FALSE
        )
    } else {
        qtab <- data.frame(
            Chromosome  = qtlm[, 1],
            Marker      = qtlm[, 3],
            "dist(cM)"  = as.numeric(qtlm[, 4]),
            Size        = round(qtle, 4),
            Prob        = pvalue,
            "Perc.Var"  = perc.var,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    }
    if (LOD)
        qtab$LOD <- lod
    # Sort by chromosome then by inferred/marker position (3rd column for interval, 2nd for marker)
    pos.col <- if (object$QTL$type == "interval") 3 else 2
    qtab <- qtab[order(as.numeric(qtab$Chromosome),
                       as.numeric(qtab[, pos.col]),
                       na.last = TRUE), ]
    rownames(qtab) <- NULL
    qtab
}
