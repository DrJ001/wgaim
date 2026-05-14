# =============================================================================
# summary.gwasAim.R
# S3 summary method for gwasAim objects.
# Returns a data frame: Chromosome, Marker, dist(cM), Size, Pvalue, %Var, LOD
# =============================================================================

#' @describeIn gwasAim Produce a summary table of significant markers, sorted
#'   by chromosome and cM position. Returns a \code{data.frame} with columns
#'   for chromosome, marker name, cM position, effect size, p-value, percentage
#'   of phenotypic variance explained, and (optionally) LOD score. The
#'   significance threshold and marker count are stored as an attribute.
#' @param object A \code{gwasAim} object.
#' @param genObj The \code{"wgPanel"} object passed to \code{gwasAim},
#'   produced by \code{\link{primePanel}}.
#' @param LOD Logical. If \code{TRUE} (default), a LOD score column is appended.
#' @export
summary.gwasAim <- function(object, genObj, LOD = TRUE, ...) {
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgPanel"))
        stop("genObj must be of class \"wgPanel\" produced by primePanel().")
    if (is.null(qtle <- object$QTL$effects)) {
        cat("No significant markers detected.\n")
        return(invisible(NULL))
    }
    sigma2 <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1

    # Always marker mode for GWAS
    gdat     <- lapply(genObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    gterm    <- object$QTL$diag$genetic.term
    scale    <- object$QTL$diag$rel.scale
    dimnames(genoData) <- list(as.character(genObj$pheno[[gterm]]),
                               names(object$QTL$diag$state))
    genoSub <- genoData[, as.logical(object$QTL$diag$state), drop = FALSE]
    if ("Gsave" %in% names(object$mf)) gterm <- "Gsave"
    genoSub <- genoSub[rownames(genoSub) %in% levels(object$mf[[gterm]]), , drop = FALSE]

    coef.mark  <- mean(apply(genoSub, 1, function(el) sum(el * el)), na.rm = TRUE)
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
        pvalue <- round((1 - pchisq(zrat^2, 1)) / 2, 6)
        pvalue[pvalue < 1e-6] <- "<1e-06"
    } else {
        pvalue <- round(2 * (1 - pnorm(abs(zrat))), 6)
        pvalue[as.numeric(pvalue) < 1e-6] <- "<1e-06"
    }
    lod  <- round(0.5 * log10(exp(zrat^2)), 2)
    qtlm <- getQTL(object, genObj)   # reuses existing getQTL (marker mode)

    qtab <- data.frame(
        Chromosome = qtlm[, 1],
        Marker     = qtlm[, 3],
        "dist(cM)" = as.numeric(qtlm[, 4]),
        Size       = round(qtle, 4),
        Pvalue     = pvalue,
        "Perc.Var" = perc.var,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    if (LOD) qtab$LOD <- lod

    # Sort by chromosome then position -- extract leading integer from chromosome
    # name so alphanumeric names (e.g. "Chr1", "1A") sort numerically.
    chr_lead <- as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1",
                               as.character(qtab$Chromosome)))
    qtab <- qtab[order(chr_lead,
                       as.character(qtab$Chromosome),
                       qtab[["dist(cM)"]],
                       na.last = TRUE), ]
    rownames(qtab) <- NULL

    qtab
}
