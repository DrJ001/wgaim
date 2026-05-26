# =============================================================================
# summary.gwasAim.R
# S3 summary method for gwasAim objects.
# Returns a data frame: Chromosome, Marker, dist(cM), Size, Pvalue, %Var, LOD
# =============================================================================

#' Summary of a Fitted \code{gwasAim} Model
#'
#' @description
#' Produces a summary table of significant markers detected by
#' \code{\link{gwasAim}}, sorted by chromosome and cM position.
#'
#' The table contains chromosome, marker name, cM position, estimated
#' additive effect, p-value, percentage of phenotypic variance explained, and
#' (optionally) a LOD score column. Returns \code{NULL} invisibly (with a
#' message) when no significant markers were detected.
#'
#' @param object A fitted object of class \code{"gwasAim"}, as returned by
#'   \code{\link{gwasAim}}.
#' @param genObj The \code{"wgPanel"} object used in the original
#'   \code{\link{gwasAim}} call, produced by \code{\link{primePanel}}.
#'   Required to resolve marker names and cM positions.
#' @param LOD Logical. If \code{TRUE} (default), a LOD score column is
#'   appended to the summary table. LOD is computed as
#'   \eqn{0.5 \log_{10}(\exp(z^2))}.
#' @param \dots Currently unused.
#'
#' @return A \code{data.frame} with one row per significant marker and columns
#'   Chromosome, Marker, dist(cM), Size, Pvalue, Perc.Var, and optionally LOD
#'   (6 or 7 columns univariate).  For multivariate models a leading
#'   \code{Trait} column is prepended (7 or 8 columns), with rows repeated
#'   once per trial level for each detected marker.
#'   Returns \code{NULL} invisibly if no markers were detected.
#'
#' @seealso \code{\link{gwasAim}}, \code{\link{print.gwasAim}},
#'   \code{\link{aimTrace}}, \code{\link{getQTL}}
#'
#' @examples
#' \dontrun{
#' # After running gwasAim():
#' summary(gwas.fit, genObj = panel)
#' summary(gwas.fit, genObj = panel, LOD = FALSE)
#' }
#'
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

    # Multivariate: use the final PRUNED model's fixed coefficients.
    # object$QTL$effects comes from the loop model which had both
    # X.chr.idx AND Trial:X.chr.idx, causing Trial1 to be aliased to 0.
    # The pruned model has either X.chr.idx (MAIN) or Trial:X.chr.idx
    # (INTERACTION with both trials non-aliased). Univariate: use effects as-is.
    is.mv <- !is.null(object$QTL$Trait)
    if (is.mv) {
        all.fc <- object$coefficients$fixed
        zind   <- grep("X\\.", rownames(all.fc))
        qtle   <- setNames(rev(all.fc[zind, 1L]), rev(rownames(all.fc)[zind]))
        veff   <- rev(object$vcoeff$fixed[zind])
    } else {
        veff <- object$QTL$veffects
    }

    enams <- names(qtle)
    qtl.x <- sub("^.*:(X\\..*)$", "\\1", enams)
    qtl.x[!grepl(":", enams)] <- enams[!grepl(":", enams)]
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
    # Interaction QTL have no bare X.chr.idx key, so use mean(effects^2).
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
        pvalue <- round((1 - pchisq(zrat^2, 1)) / 2, 6)
        pvalue[pvalue < 1e-6] <- "<1e-06"
    } else {
        pvalue <- round(2 * (1 - pnorm(abs(zrat))), 6)
        pvalue[as.numeric(pvalue) < 1e-6] <- "<1e-06"
    }
    lod  <- round(0.5 * log10(exp(zrat^2)), 2)

    # getQTL() uses object$QTL$qtl for row order and count.
    # Temporarily override it with qtl.x (converted to "Chr." prefix format) so
    # getQTL() returns exactly as many rows as qtle, in the same order --
    # including repeated rows for multivariate interaction QTL.
    orig.qtl       <- object$QTL$qtl
    object$QTL$qtl <- sub("^X\\.", "Chr.", qtl.x)
    qtlm           <- getQTL(object, genObj)
    object$QTL$qtl <- orig.qtl

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

    if (is.mv)
        qtab <- cbind(Trait = trait.lab, qtab)

    # Sort by chromosome, position, then Trait level
    pos.col  <- if (is.mv) 3L else 2L
    chr_lead <- as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1",
                               as.character(qtab$Chromosome)))
    pos_vals <- suppressWarnings(as.numeric(qtab[, pos.col + 1L]))
    sort.keys <- if (is.mv)
        list(chr_lead, as.character(qtab$Chromosome), pos_vals, qtab$Trait)
    else
        list(chr_lead, as.character(qtab$Chromosome), pos_vals)
    qtab <- qtab[do.call(order, c(sort.keys, list(na.last = TRUE))), ]
    rownames(qtab) <- NULL

    qtab
}
