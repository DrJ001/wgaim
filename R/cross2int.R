#' Convert a cross object to an interval object for use with qtlAim
#'
#' @description
#' Converts a \code{qtl} cross object into an \code{"interval"} object
#' suitable for use with \code{\link{qtlAim}}. Marker data are imputed using
#' either the Martinez-Curnow or Broman method, interval midpoint positions
#' are calculated, and optional consensus-marker collapsing is applied.
#'
#' @param object A \code{qtl} cross object. Supported population types are
#'   \code{"bc"}, \code{"dh"}, \code{"f2"}, and \code{"riself"}.
#' @param impute Character string selecting the imputation method for missing
#'   marker data. One of \code{"MartinezCurnow"} (default) or \code{"Broman"}.
#'   Partial matching is allowed.
#' @param consensus.mark Logical. If \code{TRUE} (default), co-located markers
#'   are collapsed to consensus genotypes via \code{fixMap()} before
#'   processing.
#' @param id Character string naming the column in \code{object$pheno} that
#'   uniquely identifies genotypic rows. Default \code{"id"}.
#' @param subset Character vector of chromosome names to retain. If
#'   \code{NULL} (default) all chromosomes are used.
#' @param infer Either \code{"mid"} (default) to place inferred interval
#'   markers at interval midpoints, or a numeric cM step passed to
#'   \code{lambdaf()} for finer-resolution imputation.
#'
#' @return The input cross object with an additional \code{"interval"} class
#'   and augmented \code{$geno} list elements containing \code{interval.data},
#'   \code{inferred.map}, \code{imputed.data}, \code{dist}, and \code{theta}.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{fixMap}}
#'
#' @export
cross2int <- function(object, impute = "MartinezCurnow", consensus.mark = TRUE,
    id = "id", subset = NULL, infer = "mid")
{
    cls <- class(object)[1]
    if (!(cls %in% c("bc", "dh", "f2", "riself")))
        stop("This function is restricted to populations inheriting from classes \"bc\",\"dh\",\"f2\",\"riself\".")
    object <- drop.nullmarkers(object)
    if (!(id %in% names(object$pheno)))
        stop("The unique identifier for the genotypic rows, ",
            deparse(substitute(id)), ",cannot be found in genotypic data")
    if (!is.null(subset))
        object <- subset(object, chr = subset)
    if (consensus.mark) {
        tpheno <- object$pheno
        object <- fixMap(object, rd = 3)
        object$pheno <- tpheno
    }
    lid <- as.character(object$pheno[[id]])
    mtype <- c("Broman", "MartinezCurnow")
    if (is.na(type <- pmatch(impute, mtype)))
        stop("Missing marker type must be one of \"Broman\" or \"MartinezCurnow\". Partial matching is allowed.")
    impute <- mtype[type]
    if (impute == "Broman")
        object <- argmax.geno(object)
    object$geno <- lapply(object$geno, function(el, impute, lid, cls) {
        row.names(el$data) <- as.character(lid)
        if (!is.null(el$argmax))
            el$imputed.data <- el$argmax
        else
            el$imputed.data <- el$data
        if (cls %in% "f2") {
            el$imputed.data[el$imputed.data == 3] <- -1
            el$imputed.data[el$imputed.data == 2] <- 0
        } else {
            el$imputed.data[el$imputed.data == 2] <- -1
        }
        if (length(el$map) == 1) {
            el$dist <- 0
            el$theta <- 0
            elambda <- 1/2
            names(el$dist) <- names(el$map)
            el$imputed.data[is.na(el$imputed.data)] <- 0
            el$interval.data <- as.matrix(el$imputed.data/2, ncol = 1)
            dimnames(el$interval.data)[[2]] <- names(el$map)
        } else {
            el$dist <- diff(el$map)/100
            el$theta <- 0.5 * (1 - exp(-2 * el$dist))
            if (cls %in% "riself")
                el$theta <- (el$theta/2)/(1 - el$theta)
            if (impute == "MartinezCurnow")
                el$imputed.data <- imputeGen(el$theta, el$imputed.data,
                    dom = FALSE)$add
            dimnames(el$imputed.data)[[1]] <- dimnames(el$data)[[1]]
            if (infer == "mid") {
                lend <- length(el$dist)
                elambda <- el$theta/(2 * el$dist * (1 - el$theta))
                lambda <- matrix(0, nrow = length(el$map), ncol = lend)
                inams <- paste("midm", 1:ncol(lambda), sep = "")
                mind <- cbind(c(1:lend, 1:lend + 1), c(1:lend, 1:lend))
                lambda[mind] <- rep(elambda, 2)
                el$inferred.map <- (cumsum(diff(el$map)) - diff(el$map)/2)
                names(el$inferred.map) <- inams
            } else {
                lam <- lambdaf(el, infer)
                lambda <- lam$lambda
                inams <- paste("ifm", 1:ncol(lambda), sep = "")
                el$inferred.map <- lam$imark
                names(el$inferred.map) <- inams
            }
            el$interval.data <- el$imputed.data %*% lambda
            dimnames(el$interval.data)[[2]] <- inams
        }
        el
    }, impute, lid, cls)
    for (i in seq_along(names(object$geno))) {
        wchr <- names(object$geno)[i]
        names(object$geno[[wchr]]$inferred.map) <-
            paste(wchr, names(object$geno[[wchr]]$inferred.map), sep = ".")
        colnames(object$geno[[wchr]]$interval.data) <-
            names(object$geno[[wchr]]$inferred.map)
    }
    if (length(grep("\\.", names(object$geno)))) {
        warning("Removing \".\" from linkage group names.")
        names(object$geno) <- gsub("\\.", "", names(object$geno))
    }
    class(object) <- c(class(object), "interval")
    object
}

#' Collapse co-located markers to consensus genotypes
#'
#' @description
#' Identifies markers at identical map positions (within \code{rd} decimal
#' places) and collapses each co-located group to a single consensus marker.
#' Markers where lines disagree within a bin are set to \code{NA}. The
#' redundant markers are dropped via \code{qtl::drop.markers()}.
#'
#' @param full.data A \code{qtl} cross object.
#' @param rd Integer. Number of decimal places to use when rounding map
#'   positions for co-location detection. Default \code{3}.
#'
#' @return The cross object with co-located markers collapsed and a
#'   \code{$colocated.markers} data frame appended recording which markers
#'   were merged and into which bin.
#'
#' @seealso \code{\link{cross2int}}
#'
#' @export
fixMap <- function(full.data, rd = 3)
{
    drop.mark <- lapply(full.data$geno, function(el) {
        emap <- round(el$map, rd)
        um <- unique(emap)
        dmark <- clist <- NA
        if (length(um) != length(emap)) {
            pm <- pmatch(emap, um, duplicates.ok = TRUE)
            pmt <- table(pm)
            nums <- as.numeric(names(table(pm))[pmt > 1])
            pm[!(pm %in% nums)] <- nums[length(nums)] + 1
            clist <- split(emap, pm)
            if (any(pmt == 1)) len <- length(clist) - 1
            else len <- length(clist)
            combl <- lapply(clist[1:len], function(cl) {
                names(cl)[1] <- paste(names(cl)[1], "(C)", sep = "")
                cbind.data.frame(marker = names(cl), dist = cl)
            })
            clist <- do.call("rbind", combl)
            clist$bin <- rep(1:length(combl), times = sapply(combl, nrow))
            dlist <- split.data.frame(t(el$data), pm)
            dmark <- lapply(dlist[1:len], function(dl) {
                con <- apply(dl, 2, function(ell) {
                    ell <- ell[!is.na(ell)]
                    if (length(ellu <- unique(ell)) > 1 | !length(ellu))
                        NA
                    else
                        ellu
                })
                dn <- dimnames(dl)[[1]]
                con <- matrix(con, nrow = 1, dimnames = list(dn[1]))
                dm <- dn[-1]
                list(con = con, dm = dm)
            })
            cond <- t(do.call("rbind", lapply(dmark, function(el) el$con)))
            dind <- pmatch(dimnames(cond)[[2]], dimnames(el$data)[[2]])
            cnam <- paste(dimnames(cond)[[2]], "(C)", sep = "")
            el$data[, dind] <- cond
            dimnames(el$data)[[2]][dind] <- names(el$map)[dind] <- cnam
            dmark <- unlist(lapply(dmark, function(el) el$dm))
        }
        list(chr = el, dmark = dmark, clist = clist)
    })
    full.data$geno <- lapply(drop.mark, function(dm) dm$chr)
    dmarkl <- unlist(lapply(drop.mark, function(dm) dm$dmark))
    newmap <- drop.markers(full.data, dmarkl[!is.na(dmarkl)])
    chre <- unlist(lapply(drop.mark, function(cm) any(is.na(cm))))
    cor.mark <- as.data.frame(do.call("rbind",
        lapply(drop.mark[!chre], function(cm) cm$clist)))
    chrn <- unlist(lapply(drop.mark[!chre], function(cm) nrow(cm$clist)))
    cor.mark$chr <- rep(names(nmar(full.data))[!chre], times = chrn)
    cor.mark$bin <- paste(cor.mark$chr, cor.mark$bin, sep = ".")
    rownames(cor.mark) <- NULL
    newmap$colocated.markers <- cor.mark
    newmap
}

# --- Internal helpers for cross2int ---

#' Impute missing genotype data using the Martinez-Curnow method
#'
#' @description
#' Fills missing genotype values in a marker matrix using conditional
#' expectations derived from flanking observed markers and recombination
#' fractions. Optionally also imputes dominance deviations.
#'
#' @param theta Numeric vector of recombination fractions between adjacent
#'   markers (length = number of markers - 1).
#' @param chr Numeric matrix of additive genotype scores (lines x markers),
#'   with \code{NA} for missing values.
#' @param dom Logical. If \code{TRUE} (default), also impute dominance
#'   deviations.
#'
#' @return A list with elements \code{add} (imputed additive matrix) and
#'   \code{dom} (imputed dominance matrix, or \code{NULL} if
#'   \code{dom = FALSE}).
#'
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
            a = (2*thL*(1 - thL)*thR*(1 - thR))/(1 - thLR)^2,
            b = (2*thL*(1 - thL)*thR*(1 - thR))/(thLR)^2,
            c = (thL*(1 - thL)*(1 - 2*thR*(1 - thR)))/(thLR*(1 - thLR)),
            d = (thR*(1 - thR)*(1 - 2*thL*(1 - thL)))/(thLR*(1 - thLR)),
            e = ((thL^2 + (1 - thL)^2)*(thR^2 + (1 - thR)^2))/(thLR^2 + (1 - thLR)^2)
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
                left[left %in% el]   <- left[left %in% el] - 1
                right[right %in% el] <- right[right %in% el] + 1
            }
            left[left == 0] <- right[right == n + 1] <- NA
            list(left, right)
        }, n = ncol(chr))
        left  <- unlist(sapply(lr, "[", 1))
        right <- unlist(sapply(lr, "[", 2))
        xL  <- chr[cbind(wh[, "row"], left)]
        xR  <- chr[cbind(wh[, "row"], right)]
        whc <- wh[, "col"]
        filld <- filla <- c()
        for (i in 1:nrow(wh)) {
            if (is.na(left[i]) & is.na(right[i])) {
                filld[i] <- filla[i] <- 0
            } else if (!is.na(left[i]) & is.na(right[i])) {
                tl <- th.f(theta, whc[i]:(left[i] + 1))
                filla[i] <- xL[i] * (1 - 2*tl)
                if (dom) {
                    if (abs(xL[i]) == 1) filld[i] <- 2*tl*(1 - tl)
                    else                 filld[i] <- 1 - 2*tl*(1 - tl)
                }
            } else if (is.na(left[i]) & !is.na(right[i])) {
                tr <- th.f(theta, (whc[i] + 1):(right[i]))
                filla[i] <- xR[i] * (1 - 2*tr)
                if (dom) {
                    if (abs(xR[i]) == 1) filld[i] <- 2*tr*(1 - tr)
                    else                 filld[i] <- 1 - 2*tr*(1 - tr)
                }
            } else {
                tl  <- th.f(theta, whc[i]:(left[i] + 1))
                tr  <- th.f(theta, (whc[i] + 1):right[i])
                tlr <- tl + tr - 2*tl*tr
                if (tlr == 0) {
                    filla[i] <- xL[i]
                } else {
                    lambda   <- (tr*(1 - tr)*(1 - 2*tl))/(tlr*(1 - tlr))
                    rho      <- (tl*(1 - tl)*(1 - 2*tr))/(tlr*(1 - tlr))
                    filla[i] <- xL[i]*lambda + xR[i]*rho
                }
                if (dom) {
                    if (abs(xL[i]*xR[i]) > 0) {
                        if (xL[i] == xR[i]) filld[i] <- dom.gen(tl, tr, tlr, "a")
                        else                 filld[i] <- dom.gen(tl, tr, tlr, "b")
                    } else {
                        if (abs(xL[i]) == 1)      filld[i] <- dom.gen(tl, tr, tlr, "c")
                        else if (abs(xR[i]) == 1) filld[i] <- dom.gen(tl, tr, tlr, "d")
                        else                       filld[i] <- dom.gen(tl, tr, tlr, "e")
                    }
                }
            }
        }
        chr[wh] <- filla
        if (dom) chrd[wh] <- filld
    }
    list(add = chr, dom = chrd)
}

#' Compute interval imputation weights at arbitrary cM positions
#'
#' @description
#' Given a chromosome element and a cM step size, computes the lambda weight
#' matrix used to project imputed marker genotypes onto a grid of inferred
#' positions. Used internally by \code{cross2int()} when \code{infer} is
#' numeric rather than \code{"mid"}.
#'
#' @param el A single chromosome element from \code{object$geno}, containing
#'   at least \code{el$map} and \code{el$dist}.
#' @param infer Numeric. The cM step size for the inferred position grid.
#'
#' @return A list with elements \code{lambda} (weight matrix,
#'   markers x inferred positions) and \code{imark} (inferred cM positions).
#'
#' @keywords internal
lambdaf <- function(el, infer)
{
    ints <- seq(0, el$map[length(el$map)], by = infer)
    cs   <- split(ints, cut(ints, el$map, include.lowest = TRUE))
    ri   <- findInterval(ints, el$map)
    ci   <- 1:length(ri)
    lm   <- el$map[1:(length(el$map) - 1)]
    rm   <- el$map[2:(length(el$map))]
    lambda <- matrix(0, nrow = length(el$map), ncol = length(ints))
    laml <- lamr <- list()
    k <- 1
    for (i in 1:length(cs)) {
        if (length(cs[[i]])) {
            ltq <- 0.5 * (1 - exp(-2 * (cs[[i]] - lm[i])/100))
            it  <- 0.5 * (1 - exp(-2 * (rm[i]  - lm[i])/100))
            den <- it * (1 - it) * (1 - 2*ltq)
            laml[[k]] <- (1 - it - ltq) * (it - ltq) / den
            lamr[[k]] <- ltq * (1 - ltq) * (1 - 2*it) / den
            k <- k + 1
        }
    }
    index <- rbind(cbind(ri, ci), cbind(ri + 1, ci))
    lambda[index] <- c(unlist(laml), unlist(lamr))
    list(lambda = lambda, imark = ints)
}
