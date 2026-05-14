# =============================================================================
# summary.qtlAim.R
# S3 summary method for qtlAim objects.
# =============================================================================

#' Summary of a Fitted \code{qtlAim} Model
#'
#' @description
#' Produces a detailed summary table of the significant QTL detected by
#' \code{\link{qtlAim}}, sorted by chromosome and position.
#'
#' For interval-type analyses the table contains chromosome, left flanking
#' marker and its cM position, inferred interval midpoint and its cM position,
#' right flanking marker and its cM position, estimated additive effect,
#' p-value, and percentage of phenotypic variance explained. For marker-type
#' analyses the flanking-marker columns are replaced by a single marker name
#' and cM position. An optional LOD score column can be appended.
#'
#' Returns \code{NULL} invisibly (with a message) when no significant QTL
#' were detected.
#'
#' @param object A fitted object of class \code{"qtlAim"}, as returned by
#'   \code{\link{qtlAim}}.
#' @param genObj The \code{"wgCross"} object used in the original
#'   \code{\link{qtlAim}} call, produced by \code{\link{primeCross}}.
#'   Required to resolve flanking marker names and cM positions.
#' @param LOD Logical. If \code{TRUE} (default), a LOD score column is
#'   appended to the summary table. LOD is computed from the Wald z-statistic
#'   as \eqn{0.5 \log_{10}(\exp(z^2))}.
#' @param \dots Currently unused.
#'
#' @return A \code{data.frame} with one row per detected QTL and columns
#'   depending on \code{gen.type}:
#' \describe{
#'   \item{Interval type (11 or 12 columns)}{Chromosome, Left Marker,
#'     dist(cM), Infer. Marker, dist(cM), Right Marker, dist(cM), Size,
#'     Prob, Perc.Var, and optionally LOD.}
#'   \item{Marker type (6 or 7 columns)}{Chromosome, Marker, dist(cM),
#'     Size, Prob, Perc.Var, and optionally LOD.}
#' }
#' Returns \code{NULL} invisibly if no QTL were detected.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{print.qtlAim}},
#'   \code{\link{aimTrace}}, \code{\link{getQTL}}
#'
#' @examples
#' \dontrun{
#' # After running qtlAim():
#' summary(qtl.fit, genObj = genoRxK)
#' summary(qtl.fit, genObj = genoRxK, LOD = FALSE)
#' }
#'
#' @export
summary.qtlAim <- function(object, genObj, LOD = TRUE, ...) {
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("genObj must be of class \"wgCross\" produced by primeCross().")
    if (is.null(qtle <- object$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
        return()
    }
    sigma2 <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1
    if (object$QTL$type == "interval")
        gdat <- lapply(genObj$geno, function(el) el$interval.data)
    else
        gdat <- lapply(genObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    gterm <- object$QTL$diag$genetic.term
    scale <- object$QTL$diag$rel.scale
    dimnames(genoData) <- list(as.character(genObj$pheno[[gterm]]),
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
    qtlm <- getQTL(object, genObj)
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
    # Sort by chromosome then by inferred/marker position.
    # Extract the leading integer from the chromosome name as the primary key
    # so that alphanumeric names like "1A", "2B", "3D" sort numerically
    # (1 before 2 before 3) without coercion warnings.  Full chromosome name
    # breaks ties within the same number (e.g. "3A" before "3B" before "3D").
    pos.col  <- if (object$QTL$type == "interval") 3 else 2
    chr_lead <- as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1",
                               as.character(qtab$Chromosome)))
    pos_vals <- suppressWarnings(as.numeric(qtab[, pos.col]))
    qtab <- qtab[order(chr_lead,
                       as.character(qtab$Chromosome),
                       pos_vals,
                       na.last = TRUE), ]
    rownames(qtab) <- NULL
    qtab
}
