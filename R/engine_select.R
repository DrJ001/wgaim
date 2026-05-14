# =============================================================================
# engine_select.R
# Internal functions for interval/marker selection and significance testing.
# =============================================================================

#' @keywords internal
# Trait    : character name of the trait factor column (NULL = univariate)
# ntrait   : integer, length(levels(phenoData[[Trait]])); 1 when Trait = NULL
.qtlSelect <- function(asm, phenoData, intervalObj, gen.type, selection,
                       exclusion.window, state, verboseLev,
                       Trait = NULL, ntrait = 1L) {
    sigma2 <- asm$sigma2
    if (asm$vparameters.con[length(asm$vparameters.con)] == 4)
        sigma2 <- 1

    # ------------------------------------------------------------------
    # Compute Ga (ntrait x ntrait genetic covariance matrix) and Ginv.
    # For ntrait == 1: Ga = avar (scalar wrapped in 1x1 matrix).
    # For ntrait > 1: Ga extracted from vparameters based on variance
    #   structure (diag / corh / fa(k)).
    # ------------------------------------------------------------------
    sterms <- "vm.*covObj|mbf.*ints"
    vp     <- asm$vparameters

    if (ntrait == 1L) {
        avar <- vp[grep(sterms, names(vp))] * sigma2
        Ga   <- matrix(avar, 1L, 1L)
    } else {
        rterms.all <- attr(terms.formula(asm$call$random), "term.labels")
        mterm      <- rterms.all[grep(sterms, rterms.all)]
        vpars      <- vp[grep(sterms, names(vp))] * sigma2
        if (grepl("diag", mterm)) {
            Ga <- diag(vpars, ntrait)
        } else if (grepl("corh", mterm)) {
            # corh(Trait): vparameters = c(var1, var2, cor12)
            sds <- sqrt(vpars[1:2])
            Ga  <- diag(vpars[1:2])
            Ga[1L, 2L] <- Ga[2L, 1L] <- vpars[3L] * prod(sds)
        } else if (grepl("fa", mterm)) {
            n.fa <- as.integer(sub(".*,\\s*(\\d+)\\).*", "\\1", mterm))
            psi  <- vpars[1:ntrait]
            Lam  <- matrix(vpars[(ntrait + 1L):((n.fa + 1L) * ntrait)], ncol = n.fa)
            Ga   <- Lam %*% t(Lam) + diag(psi)
        } else {
            # Fallback: treat as unstructured
            Ga <- diag(vpars[1:ntrait])
        }
        avar <- NA_real_   # not used directly when ntrait > 1
    }
    Ginv <- if (ntrait == 1L) matrix(1 / Ga[1L, 1L], 1L, 1L) else MASS::ginv(Ga)

    # ------------------------------------------------------------------
    # Predict conditional BLUPs (atilde) and compute posterior variance
    # ------------------------------------------------------------------
    cat(" Predict step for outlier statistics \n")
    cat("=====================================\n")

    if (!is.null(cov.env <- attr(intervalObj, "env"))) {
        # vm path
        rterms.all <- attr(terms.formula(asm$call$random), "term.labels")
        vmterm.raw <- rterms.all[grep("vm.*covObj", rterms.all)]
        # vmterm.raw is now "diag(Trait):vm(merge.by, covObj)" for multivariate,
        # or "vm(merge.by, covObj)" for univariate.
        # For predict():
        #   classify -- the factor cross to predict over: "merge.by" (univariate)
        #               or "merge.by:Trait" (multivariate)
        #   only     -- the full random term label as it appears in the formula

        if (ntrait == 1L) {
            pv      <- predict(asm, classify = vmterm.raw,
                               only = vmterm.raw, vcov = TRUE, data = phenoData)
            atilde  <- pv$pvals[, "predicted.value"]
            qtilde  <- as.vector(cov.env$trans %*% atilde)
            vatilde <- avar * cov.env$relm - as.matrix(pv$vcov)
            qhalf   <- cov.env$trans %*% vatilde
            vqtilde <- colSums(t(qhalf) * t(cov.env$trans))
        } else {
            # classify: the two factors involved -- merge.by and Trait
            # only: the full term label (including diag/corh/fa prefix)
            classify.term <- paste(merge.by, Trait, sep = ":")
            pv   <- predict(asm, classify = classify.term,
                            only = vmterm.raw, vcov = TRUE, data = phenoData)
            # Sort: merge.by-major order (all Trait levels for line 1, then line 2 ...)
            ord    <- order(pv$pvals[[merge.by]], pv$pvals[[Trait]])
            atilde <- matrix(pv$pvals[ord, "predicted.value"],
                             ncol = ntrait, byrow = FALSE)
            atilde[is.na(atilde)] <- 0
            pev    <- as.matrix(pv$vcov)[ord, ord]
            pev[is.na(pev)] <- 0
            # Back-transform via Cholesky inverse: qtilde.mat is nmarkers x ntrait
            qtilde.mat <- cov.env$trans %*% atilde
            vatilde    <- kronecker(Ga, cov.env$relm) - pev
            vqtilde    <- .compute_vqtilde(cov.env$trans, Ginv, vatilde, ntrait)
            # Scalar outlier stat: t(qtilde_i) %*% Ginv %*% qtilde_i
            qtilde <- apply(qtilde.mat, 1L, function(q) sum(q * (Ginv %*% q)))
        }
    } else {
        # mbf path
        if (ntrait == 1L) {
            avar    <- vp[grep("mbf.*ints", names(vp))] * sigma2
            mbf     <- grep("mbf", rownames(asm$coefficients$random))
            qtilde  <- asm$coefficients$random[mbf, 1L]
            pevar   <- sigma2 * asm$vcoeff$random[mbf]
            vqtilde <- avar - pevar
        } else {
            # mbf multivariate: coefficients are stacked [trait1 lines, trait2 lines, ...]
            mbf     <- grep("mbf", rownames(asm$coefficients$random))
            nlines  <- length(mbf) %/% ntrait
            atilde  <- matrix(asm$coefficients$random[mbf, 1L],
                               nrow = nlines, ncol = ntrait, byrow = FALSE)
            pevar.v <- sigma2 * asm$vcoeff$random[mbf]
            # Approximate vatilde as kronecker diagonal
            vatilde <- kronecker(Ga, diag(nlines)) -
                       diag(pevar.v)
            qtilde  <- apply(atilde, 1L, function(q) sum(q * (Ginv %*% q)))
            vqtilde <- vapply(seq_len(nlines), function(i) {
                Ti   <- kronecker(diag(ntrait), matrix(diag(nlines)[i, ], nrow = 1L))
                tmp2 <- Ti %*% vatilde %*% t(Ti)
                sum(diag(Ginv %*% tmp2))
            }, numeric(1L))
        }
    }

    # ------------------------------------------------------------------
    # Outlier statistics (scalar per interval/marker regardless of ntrait)
    # ------------------------------------------------------------------
    gnams <- names(state)[as.logical(state)]
    if (ntrait == 1L) {
        names(qtilde) <- names(vqtilde) <- gnams
        oint <- ifelse(!is.na(qtilde^2 / vqtilde), qtilde^2 / vqtilde, 0)
    } else {
        names(qtilde) <- names(vqtilde) <- gnams
        oint <- ifelse(!is.na(qtilde / vqtilde), qtilde / vqtilde, 0)
    }
    names(oint) <- gnams
    ochr <- NULL

    if (selection == "chromosome") {
        chr.names <- names(intervalObj$geno)
        nochr     <- length(chr.names)
        allc      <- sapply(strsplit(gnams, "\\."), "[", 2)
        ochr      <- numeric(nochr)
        for (c in seq_len(nochr)) {
            whc     <- allc %in% chr.names[c]
            if (ntrait == 1L) {
                nums <- qtilde[whc]^2
            } else {
                nums <- qtilde[whc]   # already quadratic forms
            }
            dens    <- vqtilde[whc]
            ochr[c] <- ifelse(!is.na(sum(nums) / sum(dens)), sum(nums) / sum(dens), 0)
        }
        names(ochr) <- chr.names
        mchr <- chr.names[ochr == max(ochr)]
        cint <- allc %in% mchr
        chri <- oint[cint]
        mint <- (1:length(chri))[chri == max(chri)]
        qtl  <- names(chri)[mint]
        if (verboseLev > 0) {
            cgen <- ifelse(gen.type == "marker", "Marker", "Interval")
            cat("\n Selection of chromosome using the AOM statistic\n")
            cat("=============================================== \n")
            for (i in seq_len(nochr))
                cat(" Chromosome ", chr.names[i], "Outlier Statistic ", ochr[i], "\n")
            cat("============================================= \n\n")
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for (i in seq_along(chri))
                cat(cgen, names(chri)[i], "Outlier Statistic ", chri[i], "\n")
            cat("=============================================== \n\n")
        }
    } else {
        qtl <- names(oint)[oint == max(oint)]
        qsp <- unlist(strsplit(qtl, split = "\\."))
        mint <- as.numeric(qsp[3])
        mchr <- qsp[2]
        if (verboseLev > 0) {
            cgen <- ifelse(gen.type == "marker", "Marker", "Interval")
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for (i in seq_along(oint))
                cat(cgen, names(oint)[i], "Outlier Statistic ", oint[i], "\n")
            cat("=============================================== \n\n")
        }
    }

    # Take first if tied; fill out full state vectors; apply exclusion window
    qtl   <- qtl[1]
    blups <- tint <- state
    tint[as.logical(state)]  <- oint
    # blups: for ntrait==1 standardised BLUP; for ntrait>1 sqrt of quad form
    blups[as.logical(state)] <- sqrt(abs(qtilde)) / sqrt(pmax(abs(vqtilde), .Machine$double.eps))
    oint  <- tint
    schr  <- sapply(strsplit(names(state), "\\."), "[", 2)
    wnams <- names(state)[schr %in% mchr]
    inums <- as.numeric(sapply(strsplit(wnams, "\\."), "[", 3))
    if (gen.type == "marker")
        dists <- intervalObj$geno[[mchr]]$map
    else
        dists <- intervalObj$geno[[mchr]]$inferred.map
    dists <- dists[inums]
    exc   <- wnams[abs(dists - dists[mint]) <= exclusion.window]
    state[exc] <- 0
    list(state = state, qtl = qtl, ochr = ochr, oint = oint, blups = blups)
}

#' @keywords internal
# ntrait == 1 : standard one-sided boundary LRT on a scalar variance parameter.
# ntrait >  1 : boundary LRT on a covariance matrix; uses pchisq.mixture().
.lrtTest <- function(qtlModel, baseModel, TypeI, ntrait = 1L) {
    baseLogL <- baseModel$loglik
    stat     <- 2 * (qtlModel$loglik - baseLogL)
    if (ntrait == 1L) {
        pvalue <- (1 - pchisq(stat, 1)) / 2
        pass   <- stat >= qchisq(1 - 2 * TypeI, 1)
    } else {
        pvalue <- 1 - pchisq.mixture(stat, ntrait = ntrait)
        pass   <- pvalue <= TypeI
    }
    list(baseLogL = baseLogL, stat = stat, pvalue = pvalue, pass = pass)
}
