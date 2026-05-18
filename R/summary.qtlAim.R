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

    # ------------------------------------------------------------------
    # Multivariate: use the final PRUNED model's fixed coefficients.
    #
    # object$QTL$effects comes from coef.list[[iter-1]] — the loop model
    # which had both X.chr.idx AND Trial:X.chr.idx simultaneously, causing
    # ASReml to alias Trial1's interaction coefficient to zero.
    # After waldTest pruning each QTL is either:
    #   MAIN      : X.chr.idx only         (one coefficient, not aliased)
    #   INTERACTION: Trial:X.chr.idx only  (ntrait coefficients, none aliased)
    # so object$coefficients$fixed contains the correct non-aliased values.
    #
    # Univariate: no Trait slot -> use object$QTL$effects as before.
    # ------------------------------------------------------------------
    is.mv <- !is.null(object$QTL$Trait)
    if (is.mv) {
        all.fc <- object$coefficients$fixed
        zind   <- grep("X\\.", rownames(all.fc))
        qtle   <- setNames(rev(all.fc[zind, 1L]), rev(rownames(all.fc)[zind]))
        veff   <- rev(object$vcoeff$fixed[zind])
    } else {
        veff <- object$QTL$veffects
    }

    enams  <- names(qtle)
    qtl.x  <- sub("^.*:(X\\..*)$", "\\1", enams)   # always extract X.chr.idx part
    qtl.x[!grepl(":", enams)] <- enams[!grepl(":", enams)]  # main effects unchanged
    trait.lab <- rep("MAIN", length(enams))
    if (is.mv) {
        int.rows <- grepl(":", enams)
        trait.lab[int.rows] <- sub("^(.*):X\\..*$", "\\1", enams[int.rows])
        prefix <- paste0(object$QTL$Trait, "_")
        trait.lab <- gsub(prefix, "", trait.lab, fixed = TRUE)
    }

    if (object$QTL$method == "random") {
        var.est  <- sigma2 * object$vparameters[grep("X\\.", names(object$vparameters))]
        coef.est <- apply(genoData[, object$QTL$qtl, drop = FALSE]^2, 2, mean, na.rm = TRUE)
    } else {
        var.est  <- qtle^2
        coef.est <- rep(1, length(qtle))
    }
    # Perc.var: one value per unique QTL, broadcast to all its rows.
    # For MAIN QTL:        var contribution = effect^2
    # For INTERACTION QTL: var contribution = mean(per-trial effects^2)
    # (no bare X.chr.idx key exists for interaction QTL in the final model)
    # var.mark is summed across traits for multivariate to avoid length mismatch.
    unique.qtl.x <- unique(qtl.x)
    if (object$QTL$method == "random") {
        var.est.qtl  <- var.est[unique.qtl.x]
        coef.est.qtl <- coef.est[unique.qtl.x]
    } else {
        var.est.qtl  <- vapply(unique.qtl.x, function(uq)
            mean(qtle[qtl.x == uq]^2), numeric(1L))
        coef.est.qtl <- rep(1, length(unique.qtl.x))
    }
    # Scale all variance components to per-trial averages before computing
    # Perc.Var.  var.mark and var.res are vectors of length ntrait for
    # multivariate models (one value per trial); summing them would inflate
    # the denominator relative to the per-trial-average var.est.qtl numerator.
    # Using mean() keeps all three components on the same per-trial scale so
    # the percentages are interpretable and sum to <= 100%.
    var.mark.scalar  <- mean(var.mark)
    var.res.scalar   <- mean(var.res)
    var.all  <- sum(var.est.qtl) + coef.mark * var.mark.scalar + var.res.scalar
    perc.var.per.qtl <- setNames(
        round(100 * var.est.qtl / var.all, 1), unique.qtl.x)
    perc.var <- perc.var.per.qtl[match(qtl.x, unique.qtl.x)]

    zrat <- qtle / sqrt(veff * sigma2)
    if (object$QTL$method == "random") {
        pvalue <- round((1 - pchisq(zrat^2, 1)) / 2, 4)
        pvalue[pvalue < 0.0001] <- "<0.0001"
    } else {
        pvalue <- round(2 * (1 - pnorm(abs(zrat))), 4)
        pvalue[as.numeric(pvalue) < 0.0001] <- "<0.0001"
    }
    lod <- round(0.5 * log10(exp(zrat^2)), 2)

    # getQTL() uses object$QTL$qtl for row order and count.
    # Temporarily override it with qtl.x (converted to "Chr." prefix format) so
    # getQTL() returns exactly as many rows as qtle, in the same order --
    # including repeated rows for multivariate interaction QTL.
    orig.qtl       <- object$QTL$qtl
    object$QTL$qtl <- sub("^X\\.", "Chr.", qtl.x)
    qtlm           <- getQTL(object, genObj)
    object$QTL$qtl <- orig.qtl

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

    # For multivariate analyses, prepend the Trait column
    if (is.mv)
        qtab <- cbind(Trait = trait.lab, qtab)

    # Sort by chromosome then by inferred/marker position (then Trait level).
    # Extract leading integer from chromosome name for numeric-safe sorting.
    pos.col  <- if (is.mv) {
        if (object$QTL$type == "interval") 4L else 3L
    } else {
        if (object$QTL$type == "interval") 3L else 2L
    }
    chr_col  <- if (is.mv) "Chromosome" else "Chromosome"
    chr_lead <- as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1",
                               as.character(qtab$Chromosome)))
    pos_vals <- suppressWarnings(as.numeric(qtab[, pos.col]))
    sort.keys <- if (is.mv)
        list(chr_lead, as.character(qtab$Chromosome), pos_vals, qtab$Trait)
    else
        list(chr_lead, as.character(qtab$Chromosome), pos_vals)
    qtab <- qtab[do.call(order, c(sort.keys, list(na.last = TRUE))), ]
    rownames(qtab) <- NULL
    qtab
}
