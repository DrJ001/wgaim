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
    # Multivariate: decompose effects vector into QTL x Trait rows.
    # This must happen before perc.var so qtl.x is available.
    # effect names are either "X.chr.idx" (main) or "Trait_level:X.chr.idx"
    # (interaction). We label each row with its Trait level, using "MAIN"
    # for main-effect-only QTL.
    # Univariate: no Trait slot -> all rows labelled "MAIN" (column suppressed).
    # ------------------------------------------------------------------
    is.mv  <- !is.null(object$QTL$Trait)
    enams  <- names(qtle)
    qtl.x  <- sub("^.*:(X\\..*)$", "\\1", enams)   # always extract X.chr.idx part
    qtl.x[!grepl(":", enams)] <- enams[!grepl(":", enams)]  # main effects unchanged
    trait.lab <- rep("MAIN", length(enams))
    if (is.mv) {
        int.rows <- grepl(":", enams)
        trait.lab[int.rows] <- sub(paste0("^(.*):X\\..*$"), "\\1", enams[int.rows])
        # Strip the Trait column name prefix (e.g. "Site_A" -> "A")
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
    # For multivariate analyses the effects vector contains both main and
    # interaction rows sharing the same underlying QTL. Compute perc.var once
    # per unique QTL (using the main-effect row) then broadcast to all rows.
    unique.qtl.x <- unique(qtl.x)   # bare X.chr.idx keys, one per detected QTL
    var.est.qtl  <- if (object$QTL$method == "random")
        var.est[unique.qtl.x]
    else
        qtle[unique.qtl.x]^2
    coef.est.qtl <- if (object$QTL$method == "random")
        coef.est[unique.qtl.x]
    else
        rep(1, length(unique.qtl.x))
    var.all  <- sum(c(coef.est.qtl, coef.mark, 1) * c(var.est.qtl, var.mark, var.res))
    perc.var.per.qtl <- setNames(
        round(100 * (coef.est.qtl * var.est.qtl) / var.all, 1), unique.qtl.x)
    perc.var <- perc.var.per.qtl[match(qtl.x, unique.qtl.x)]

    zrat <- qtle / sqrt(object$QTL$veffects * sigma2)
    if (object$QTL$method == "random") {
        pvalue <- round((1 - pchisq(zrat^2, 1)) / 2, 4)
        pvalue[pvalue < 0.0001] <- "<0.0001"
    } else {
        pvalue <- round(2 * (1 - pnorm(abs(zrat))), 4)
        pvalue[as.numeric(pvalue) < 0.0001] <- "<0.0001"
    }
    lod <- round(0.5 * log10(exp(zrat^2)), 2)

    # Temporarily remap QTL$effects names to bare X.chr.idx for getQTL()
    orig.effects <- object$QTL$effects
    object$QTL$effects <- setNames(qtle, qtl.x)
    qtlm <- getQTL(object, genObj)
    object$QTL$effects <- orig.effects   # restore

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
