# =============================================================================
# primeCross.R
# Prepare a qtl cross object for whole-genome analysis with qtlAim / gpAim.
#
# primeCross() enriches a qtl cross object with imputed genotypes and
# optionally interval midpoint data, returning a "wgCross" object ready
# for the wgAim analysis functions.
#
# Helper functions (not exported):
#   imputeGen() -- Martinez-Curnow genotype imputation
#   lambdaf()   -- interval weight matrix for arbitrary cM grids
# =============================================================================

#' Prepare a Cross Object for Whole-Genome Analysis
#'
#' @description
#' Converts a \code{qtl} cross object into a \code{"wgCross"} object suitable
#' for use with \code{\link{qtlAim}} and \code{\link{gpAim}}.
#'
#' Two analysis types are supported, selected via the \code{type} argument:
#'
#' \describe{
#'   \item{\code{"interval"} (default)}{Missing marker genotypes are imputed,
#'     co-located markers are collapsed via \code{\link{fixMap}}, and interval
#'     midpoint scores are computed for each inter-marker interval. This is
#'     the recommended type for QTL detection in biparental populations.}
#'   \item{\code{"marker"}}{Missing marker genotypes are imputed and
#'     co-located markers are collapsed, but no interval data are generated.
#'     Use this type when working directly with marker scores rather than
#'     interval midpoints.}
#' }
#'
#' The resulting object retains the full \code{qtl} cross structure and can
#' still be used with all \code{qtl} package functions.
#'
#' @param object A \code{qtl} cross object. Supported population types are
#'   \code{"bc"}, \code{"dh"}, \code{"f2"}, and \code{"riself"}.
#' @param type Character string, either \code{"interval"} (default) or
#'   \code{"marker"}. Controls whether interval markers are generated
#'   according to the rules of the argument \code{infer} below in addition to imputed marker scores.
#' @param impute Character string selecting the imputation method for missing
#'   marker genotypes. One of \code{"MartinezCurnow"} (default) or
#'   \code{"Broman"}. Partial matching is allowed.
#' @param consensus.mark Logical. If \code{TRUE} (default), co-located markers
#'   are collapsed to consensus genotypes via \code{\link{fixMap}} before
#'   imputation.
#' @param id Character string naming the column in \code{object$pheno} that
#'   uniquely identifies genotypic rows. Default \code{"id"}.
#' @param subset Character vector of chromosome names to retain. If
#'   \code{NULL} (default) all chromosomes are used.
#' @param infer Either \code{"mid"} (default) to place inferred interval
#'   markers at interval midpoints, or a numeric cM step size passed to
#'   \code{lambdaf()} for a finer-resolution grid. Only used when
#'   \code{type = "interval"}.
#'
#' @return An object of class \code{c("wgCross", cls, "cross")} where
#'   \code{cls} is the original population class (\code{"bc"}, \code{"dh"},
#'   etc.). The \code{$geno} list elements are augmented with:
#'   \describe{
#'     \item{\code{$imputed.data}}{Imputed marker genotype matrix (lines x
#'       markers) in additive \eqn{\pm 1} coding.}
#'     \item{\code{$dist}}{Inter-marker distances in Morgans.}
#'     \item{\code{$theta}}{Recombination fractions between adjacent markers.}
#'     \item{\code{$interval.data}}{Interval midpoint score matrix (lines x
#'       intervals). Only present when \code{type = "interval"}.}
#'     \item{\code{$inferred.map}}{Named numeric vector of inferred interval
#'       midpoint positions in cM. Only present when \code{type = "interval"}.}
#'   }
#'   The \code{$type} field records the value of \code{type} for use by
#'   downstream analysis functions.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{gpAim}}, \code{\link{fixMap}}
#'
#' @export
primeCross <- function(object, type = "interval", impute = "MartinezCurnow",
    consensus.mark = TRUE, id = "id", subset = NULL, infer = "mid")
{
    type <- match.arg(type, c("interval", "marker"))
    cls  <- class(object)[1]
    if (!(cls %in% c("bc", "dh", "f2", "riself")))
        stop("This function is restricted to populations inheriting from ",
             "classes \"bc\", \"dh\", \"f2\", \"riself\".")
    object <- drop.nullmarkers(object)
    if (!(id %in% names(object$pheno)))
        stop("The unique identifier for the genotypic rows, ",
             deparse(substitute(id)),
             ", cannot be found in the genotypic data.")
    if (!is.null(subset))
        object <- subset(object, chr = subset)
    if (consensus.mark) {
        tpheno <- object$pheno
        object <- fixMap(object, rd = 3)
        object$pheno <- tpheno
    }
    lid   <- as.character(object$pheno[[id]])
    mtype <- c("Broman", "MartinezCurnow")
    if (is.na(tp <- pmatch(impute, mtype)))
        stop("impute must be one of \"Broman\" or \"MartinezCurnow\". ",
             "Partial matching is allowed.")
    impute <- mtype[tp]
    if (impute == "Broman")
        object <- argmax.geno(object)

    object$geno <- lapply(object$geno, function(el, impute, lid, cls, type,
                                                 infer) {
        row.names(el$data) <- as.character(lid)
        if (!is.null(el$argmax))
            el$imputed.data <- el$argmax
        else
            el$imputed.data <- el$data
        if (cls %in% "f2") {
            el$imputed.data[el$imputed.data == 3] <- -1
            el$imputed.data[el$imputed.data == 2] <-  0
        } else {
            el$imputed.data[el$imputed.data == 2] <- -1
        }
        if (length(el$map) == 1) {
            el$dist  <- 0
            el$theta <- 0
            names(el$dist) <- names(el$map)
            el$imputed.data[is.na(el$imputed.data)] <- 0
            if (type == "interval") {
                el$interval.data <- as.matrix(el$imputed.data / 2, ncol = 1)
                dimnames(el$interval.data)[[2]] <- names(el$map)
                el$inferred.map  <- setNames(el$map, names(el$map))
            }
        } else {
            el$dist  <- diff(el$map) / 100
            el$theta <- 0.5 * (1 - exp(-2 * el$dist))
            if (cls %in% "riself")
                el$theta <- (el$theta / 2) / (1 - el$theta)
            if (impute == "MartinezCurnow")
                el$imputed.data <- imputeGen(el$theta, el$imputed.data,
                                             dom = FALSE)$add
            dimnames(el$imputed.data)[[1]] <- dimnames(el$data)[[1]]
            if (type == "interval") {
                if (infer == "mid") {
                    lend    <- length(el$dist)
                    elambda <- el$theta / (2 * el$dist * (1 - el$theta))
                    lambda  <- matrix(0, nrow = length(el$map), ncol = lend)
                    inams   <- paste("midm", seq_len(ncol(lambda)), sep = "")
                    mind    <- cbind(c(seq_len(lend), seq_len(lend) + 1),
                                    c(seq_len(lend), seq_len(lend)))
                    lambda[mind]    <- rep(elambda, 2)
                    el$inferred.map <- cumsum(diff(el$map)) - diff(el$map) / 2
                    names(el$inferred.map) <- inams
                } else {
                    lam    <- lambdaf(el, infer)
                    lambda <- lam$lambda
                    inams  <- paste("ifm", seq_len(ncol(lambda)), sep = "")
                    el$inferred.map <- lam$imark
                    names(el$inferred.map) <- inams
                }
                el$interval.data <- el$imputed.data %*% lambda
                dimnames(el$interval.data)[[2]] <- inams
            }
        }
        el
    }, impute, lid, cls, type, infer)

    if (type == "interval") {
        for (i in seq_along(names(object$geno))) {
            wchr <- names(object$geno)[i]
            names(object$geno[[wchr]]$inferred.map) <-
                paste(wchr, names(object$geno[[wchr]]$inferred.map), sep = ".")
            colnames(object$geno[[wchr]]$interval.data) <-
                names(object$geno[[wchr]]$inferred.map)
        }
    }
    if (length(grep("\\.", names(object$geno)))) {
        warning("Removing \".\" from linkage group names.")
        names(object$geno) <- gsub("\\.", "", names(object$geno))
    }
    object$type  <- type
    class(object) <- c("wgCross", class(object))
    object
}

#' Collapse Co-located Markers to Consensus Genotypes
#'
#' @description
#' Identifies markers at identical map positions (within \code{rd} decimal
#' places) and collapses each co-located group to a single consensus marker.
#' Markers where lines disagree within a bin are set to \code{NA}. Redundant
#' markers are dropped via \code{qtl::drop.markers()}.
#'
#' @param full.data A \code{qtl} cross object.
#' @param rd Integer. Number of decimal places for rounding map positions
#'   when detecting co-location. Default \code{3}.
#'
#' @return The cross object with co-located markers collapsed and a
#'   \code{$colocated.markers} data frame appended recording which markers
#'   were merged and into which bin.
#'
#' @seealso \code{\link{primeCross}}
#'
#' @export
fixMap <- function(full.data, rd = 3)
{
    drop.mark <- lapply(full.data$geno, function(el) {
        emap  <- round(el$map, rd)
        um    <- unique(emap)
        dmark <- clist <- NA
        if (length(um) != length(emap)) {
            pm  <- pmatch(emap, um, duplicates.ok = TRUE)
            pmt <- table(pm)
            nums <- as.numeric(names(table(pm))[pmt > 1])
            pm[!(pm %in% nums)] <- nums[length(nums)] + 1
            clist <- split(emap, pm)
            if (any(pmt == 1)) len <- length(clist) - 1
            else               len <- length(clist)
            combl <- lapply(clist[seq_len(len)], function(cl) {
                names(cl)[1] <- paste0(names(cl)[1], "(C)")
                cbind.data.frame(marker = names(cl), dist = cl)
            })
            clist       <- do.call("rbind", combl)
            clist$bin   <- rep(seq_len(length(combl)),
                               times = sapply(combl, nrow))
            dlist <- split.data.frame(t(el$data), pm)
            dmark <- lapply(dlist[seq_len(len)], function(dl) {
                con <- apply(dl, 2, function(ell) {
                    ell <- ell[!is.na(ell)]
                    if (length(ellu <- unique(ell)) > 1 || !length(ellu)) NA
                    else ellu
                })
                dn  <- dimnames(dl)[[1]]
                con <- matrix(con, nrow = 1, dimnames = list(dn[1]))
                dm  <- dn[-1]
                list(con = con, dm = dm)
            })
            cond <- t(do.call("rbind", lapply(dmark, function(el) el$con)))
            dind <- pmatch(dimnames(cond)[[2]], dimnames(el$data)[[2]])
            cnam <- paste0(dimnames(cond)[[2]], "(C)")
            el$data[, dind]               <- cond
            dimnames(el$data)[[2]][dind]  <- names(el$map)[dind] <- cnam
            dmark <- unlist(lapply(dmark, function(el) el$dm))
        }
        list(chr = el, dmark = dmark, clist = clist)
    })
    full.data$geno <- lapply(drop.mark, function(dm) dm$chr)
    dmarkl  <- unlist(lapply(drop.mark, function(dm) dm$dmark))
    newmap  <- drop.markers(full.data, dmarkl[!is.na(dmarkl)])
    chre    <- unlist(lapply(drop.mark, function(cm) any(is.na(cm))))
    cor.mark <- as.data.frame(do.call("rbind",
        lapply(drop.mark[!chre], function(cm) cm$clist)))
    chrn    <- unlist(lapply(drop.mark[!chre], function(cm) nrow(cm$clist)))
    cor.mark$chr <- rep(names(nmar(full.data))[!chre], times = chrn)
    cor.mark$bin <- paste(cor.mark$chr, cor.mark$bin, sep = ".")
    rownames(cor.mark) <- NULL
    newmap$colocated.markers <- cor.mark
    newmap
}

# --- Internal helpers --------------------------------------------------------

#' @keywords internal
imputeGen <- function(theta, chr, dom = TRUE)
{
    th.f <- function(the, ind) {
        th <- 0
        for (i in ind)
            th <- th + the[i - 1] - 2 * th * the[i - 1]
        th
    }
    dom.gen <- function(thL, thR, thLR, wh) {
        switch(wh,
            a = (2*thL*(1-thL)*thR*(1-thR)) / (1-thLR)^2,
            b = (2*thL*(1-thL)*thR*(1-thR)) / thLR^2,
            c = (thL*(1-thL)*(1-2*thR*(1-thR))) / (thLR*(1-thLR)),
            d = (thR*(1-thR)*(1-2*thL*(1-thL))) / (thLR*(1-thLR)),
            e = ((thL^2+(1-thL)^2)*(thR^2+(1-thR)^2)) / (thLR^2+(1-thLR)^2)
        )
    }
    chrd <- NULL
    if (dom) {
        chrd <- chr + 1
        chrd[chrd %in% 2] <- 0
    }
    wh <- which(is.na(chr), arr.ind = TRUE)
    if (dim(wh)[1] != 0) {
        wh <- wh[order(wh[, "row"], wh[, "col"]), , drop = FALSE]
        sp <- split(wh[, "col"], wh[, "row"])
        lr <- lapply(sp, function(el, n) {
            left <- el - 1; right <- el + 1
            while (any(c(left, right) %in% el)) {
                left[left   %in% el] <- left[left   %in% el] - 1
                right[right %in% el] <- right[right %in% el] + 1
            }
            left[left == 0] <- right[right == n + 1] <- NA
            list(left, right)
        }, n = ncol(chr))
        left  <- unlist(sapply(lr, "[", 1))
        right <- unlist(sapply(lr, "[", 2))
        xL    <- chr[cbind(wh[, "row"], left)]
        xR    <- chr[cbind(wh[, "row"], right)]
        whc   <- wh[, "col"]
        filld <- filla <- c()
        for (i in seq_len(nrow(wh))) {
            if (is.na(left[i]) && is.na(right[i])) {
                filld[i] <- filla[i] <- 0
            } else if (!is.na(left[i]) && is.na(right[i])) {
                tl       <- th.f(theta, whc[i]:(left[i] + 1))
                filla[i] <- xL[i] * (1 - 2*tl)
                if (dom) {
                    if (abs(xL[i]) == 1) filld[i] <- 2*tl*(1-tl)
                    else                 filld[i] <- 1 - 2*tl*(1-tl)
                }
            } else if (is.na(left[i]) && !is.na(right[i])) {
                tr       <- th.f(theta, (whc[i]+1):right[i])
                filla[i] <- xR[i] * (1 - 2*tr)
                if (dom) {
                    if (abs(xR[i]) == 1) filld[i] <- 2*tr*(1-tr)
                    else                 filld[i] <- 1 - 2*tr*(1-tr)
                }
            } else {
                tl  <- th.f(theta, whc[i]:(left[i]+1))
                tr  <- th.f(theta, (whc[i]+1):right[i])
                tlr <- tl + tr - 2*tl*tr
                if (tlr == 0) {
                    filla[i] <- xL[i]
                } else {
                    lambda   <- (tr*(1-tr)*(1-2*tl)) / (tlr*(1-tlr))
                    rho      <- (tl*(1-tl)*(1-2*tr)) / (tlr*(1-tlr))
                    filla[i] <- xL[i]*lambda + xR[i]*rho
                }
                if (dom) {
                    if (abs(xL[i]*xR[i]) > 0) {
                        if (xL[i] == xR[i]) filld[i] <- dom.gen(tl, tr, tlr, "a")
                        else                filld[i] <- dom.gen(tl, tr, tlr, "b")
                    } else {
                        if (abs(xL[i]) == 1)       filld[i] <- dom.gen(tl, tr, tlr, "c")
                        else if (abs(xR[i]) == 1)  filld[i] <- dom.gen(tl, tr, tlr, "d")
                        else                        filld[i] <- dom.gen(tl, tr, tlr, "e")
                    }
                }
            }
        }
        chr[wh] <- filla
        if (dom) chrd[wh] <- filld
    }
    list(add = chr, dom = chrd)
}

#' @keywords internal
lambdaf <- function(el, infer)
{
    ints   <- seq(0, el$map[length(el$map)], by = infer)
    cs     <- split(ints, cut(ints, el$map, include.lowest = TRUE))
    ri     <- findInterval(ints, el$map)
    ci     <- seq_along(ri)
    lm     <- el$map[seq_len(length(el$map) - 1)]
    rm     <- el$map[2:length(el$map)]
    lambda <- matrix(0, nrow = length(el$map), ncol = length(ints))
    laml   <- lamr <- list()
    k      <- 1
    for (i in seq_along(cs)) {
        if (length(cs[[i]])) {
            ltq      <- 0.5 * (1 - exp(-2 * (cs[[i]] - lm[i]) / 100))
            it       <- 0.5 * (1 - exp(-2 * (rm[i]  - lm[i]) / 100))
            den      <- it * (1 - it) * (1 - 2*ltq)
            laml[[k]] <- (1 - it - ltq) * (it - ltq) / den
            lamr[[k]] <- ltq * (1 - ltq) * (1 - 2*it) / den
            k        <- k + 1
        }
    }
    index <- rbind(cbind(ri, ci), cbind(ri + 1, ci))
    lambda[index] <- c(unlist(laml), unlist(lamr))
    list(lambda = lambda, imark = ints)
}
