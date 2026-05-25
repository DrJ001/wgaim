# =============================================================================
# engine_select.R
# Internal functions for interval/marker selection and significance testing.
# =============================================================================

#' @keywords internal
# Trait    : character name of the trait factor column (NULL = univariate)
# ntrait   : integer, length(levels(phenoData[[Trait]])); 1 when Trait = NULL
.qtlSelect <- function(asm, phenoData, intervalObj, gen.type, selection,
                       exclusion.window, state, verboseLev, merge.by = NULL,
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
        } else if (grepl("corh|corgh", mterm)) {
            # corh/corgh vparameters in ASReml are stored:
            #   correlations first (names end in ".cor"), then variances.
            # Variances are ratios to sigma2 and need * sigma2.
            # Correlations are in [-1,1] and must NOT be multiplied by sigma2.
            vpar.nms  <- names(vp)[grep(sterms, names(vp))]
            vpars.raw <- vp[grep(sterms, names(vp))]
            is.cor    <- grepl("\\.cor$", vpar.nms)
            var.pars  <- pmax(vpars.raw[!is.cor] * sigma2, 0)
            cor.pars  <- vpars.raw[is.cor]
            sds  <- sqrt(var.pars)
            Ga   <- diag(var.pars)
            idx  <- 0L
            for (col in seq_len(ntrait - 1L)) {
                for (row in (col + 1L):ntrait) {
                    idx <- idx + 1L
                    Ga[row, col] <- Ga[col, row] <- cor.pars[idx] * sds[row] * sds[col]
                }
            }
        } else if (grepl("fa", mterm)) {
            # fa(Trait, k) vparameters are ordered:
            #   [1:ntrait]       specific variances Psi_j  (ratios to sigma2)
            #   [ntrait+1 : end] factor loadings Lambda_jk (raw scale, not /sigma2)
            # ASReml includes the identifiability-constrained loading (con=4,
            # value=0) in vparameters, so the total is always ntrait*(n.fa+1).
            # Separate by name: !var = Psi; !fa[k] = Lambda.
            n.fa      <- as.integer(sub(".*,\\s*(\\d+)\\).*", "\\1", mterm))
            vpar.nms  <- names(vp)[grep(sterms, names(vp))]
            vpars.raw <- vp[grep(sterms, names(vp))]
            is.loading <- grepl("!fa[0-9]+$", vpar.nms)
            psi <- pmax(vpars.raw[!is.loading] * sigma2, 0)
            Lam <- matrix(vpars.raw[is.loading], ncol = n.fa)
            Ga  <- Lam %*% t(Lam) + diag(as.numeric(psi))
            # When specific variances hit the boundary (zero), Ga is rank-deficient.
            # MASS::ginv() computes the Moore-Penrose pseudoinverse, which is the
            # correct generalised inverse Ga^- without ad-hoc regularisation.
        } else {
            # Fallback: treat as unstructured
            Ga <- diag(vpars[1:ntrait])
        }
        avar <- NA_real_   # not used directly when ntrait > 1
    }
    Ginv <- if (ntrait == 1L) matrix(1 / Ga[1L, 1L], 1L, 1L) else MASS::ginv(Ga)

    # Trait levels used for column names of scaled.mat (ntrait > 1 only)
    trait.levels <- if (!is.null(Trait)) levels(phenoData[[Trait]]) else NULL

    # scaled.mat: nmarkers x ntrait signed z-scores; populated inside vm/mbf branches
    scaled.mat <- NULL

    # ------------------------------------------------------------------
    # Predict conditional BLUPs (atilde) and compute posterior variance
    # ------------------------------------------------------------------
    cat(" Predict step for outlier statistics \n")
    cat("=====================================\n")

    if (!is.null(cov.env <- attr(intervalObj, "env"))) {
        # vm path
        rterms.all <- attr(terms.formula(asm$call$random), "term.labels")
        vmterm.raw <- rterms.all[grep("vm.*covObj", rterms.all)]
        # vmterm.only: G.param key used for predict(only=) in the UV case.
        # MV predict(only=) requires the raw formula string (vmterm.raw) because
        # ASReml matches on the formula term, not the G.param key.
        vmterm.only <- names(asm$G.param)[grep("vm.*covObj", names(asm$G.param))]

        if (ntrait == 1L) {
            pv      <- predict(asm, classify = vmterm.only,
                               only = vmterm.only, vcov = TRUE, data = phenoData, maxit = 1)
            atilde  <- pv$pvals[, "predicted.value"]
            qtilde  <- as.vector(cov.env$trans %*% atilde)
            vatilde <- avar * cov.env$relm - as.matrix(pv$vcov)
            qhalf   <- cov.env$trans %*% vatilde
            vqtilde <- colSums(t(qhalf) * t(cov.env$trans))
        } else {
            # MV: use the raw formula term string for `only`, not the G.param key.
            classify.term <- paste(Trait, merge.by, sep = ":")
            pv   <- predict(asm, classify = classify.term,
                            only = vmterm.raw, vcov = TRUE, data = phenoData, maxit = 1)
            # Sort Trait-major: all lines for Trait 1 first, then Trait 2, ...
            # This matches kronecker(Ga, relm) which is block-structured as
            # Ga[i,j] * relm in Trait-major order.
            # byrow=FALSE then correctly fills atilde: column k = all lines for trait k.
            nlines <- nrow(pv$pvals) %/% ntrait
            ord    <- order(pv$pvals[[Trait]], pv$pvals[[merge.by]])
            atilde <- matrix(pv$pvals[ord, "predicted.value"],
                             nrow = nlines, ncol = ntrait, byrow = FALSE)
            atilde[is.na(atilde)] <- 0
            pev    <- as.matrix(pv$vcov)[ord, ord]
            pev[is.na(pev)] <- 0
            # Back-transform: qtilde.mat is nmarkers x ntrait
            qtilde.mat <- cov.env$trans %*% atilde
            vatilde    <- kronecker(Ga, cov.env$relm) - pev

            # Use the Rcpp implementation for speed and correctness.
            vqtilde    <- compute_vqtilde(cov.env$trans, Ginv,
                                          as.matrix(vatilde), ntrait)
            # Scalar outlier stat: t(qtilde_i) %*% Ginv %*% qtilde_i
            qtilde <- apply(qtilde.mat, 1L, function(q) sum(q * (Ginv %*% q)))
            # Signed scaled BLUPs (nmarkers x ntrait): per-trial z-scores.
            # For trial j, posterior variance = diag(trans %*% vatilde_jj %*% t(trans))
            # where vatilde_jj is the j-th nlines x nlines diagonal block of vatilde.
            scaled.mat <- matrix(0, nrow = nrow(qtilde.mat), ncol = ntrait,
                                 dimnames = list(rownames(qtilde.mat), trait.levels))
            for (j in seq_len(ntrait)) {
                rows_j  <- ((j - 1L) * nlines + 1L):(j * nlines)
                v_j     <- rowSums((cov.env$trans %*% vatilde[rows_j, rows_j]) *
                                    cov.env$trans)
                v_j     <- pmax(v_j, 0)
                scaled.mat[, j] <- ifelse(v_j > 0, qtilde.mat[, j] / sqrt(v_j), 0)
            }
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
            # Signed scaled BLUPs (nmarkers x ntrait): per-trial z-scores.
            # vatilde diagonal blocks: block j = rows/cols ((j-1)*nlines+1):(j*nlines).
            scaled.mat <- matrix(0, nrow = nlines, ncol = ntrait,
                                 dimnames = list(NULL, trait.levels))
            for (j in seq_len(ntrait)) {
                rows_j  <- ((j - 1L) * nlines + 1L):(j * nlines)
                v_j     <- pmax(diag(vatilde[rows_j, rows_j]), 0)
                scaled.mat[, j] <- ifelse(v_j > 0, atilde[, j] / sqrt(v_j), 0)
            }
        }
    }

    # ------------------------------------------------------------------
    # Outlier statistics (scalar per interval/marker regardless of ntrait)
    # ------------------------------------------------------------------
    # vqtilde is the posterior variance diagonal: prior - PEV.  At late
    # iterations, after most markers are excluded, the remaining covObj can
    # be nearly rank-deficient and PEV can numerically exceed the prior,
    # yielding negative vqtilde values.  Floor at zero so that oint is
    # always non-negative (a negative outlier stat is meaningless).
    gnams <- names(state)[as.logical(state)]
    vqtilde <- pmax(vqtilde, 0)
    names(qtilde) <- names(vqtilde) <- gnams
    # Outlier statistic: two forms depending on ntrait.
    #
    # ntrait == 1: qtilde is a signed scalar BLUP; squaring gives a
    #   chi-squared-scale Wald statistic: oint = qtilde^2 / vqtilde.
    #
    # ntrait > 1: qtilde = t(q_i) %*% Ginv %*% q_i is already a quadratic
    #   form (>= 0, chi-squared-scale); vqtilde = tr(Ginv Ti vatilde Ti') is
    #   its variance.  No further squaring: oint = qtilde / vqtilde.
    oint <- if (ntrait == 1L)
        ifelse(vqtilde > 0, qtilde^2 / vqtilde, 0)
    else
        ifelse(vqtilde > 0, qtilde  / vqtilde, 0)
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
    qtl  <- qtl[1]
    tint <- state
    tint[as.logical(state)] <- oint
    # blups: signed z-scores for BLUP/contrast plots.
    #   ntrait == 1: named vector (length = all markers); zero where vqtilde == 0.
    #   ntrait >  1: named matrix (all markers x ntrait); zero rows for excluded markers.
    if (ntrait == 1L) {
        blups <- state
        blups[as.logical(state)] <- ifelse(
            vqtilde > 0,
            qtilde / sqrt(vqtilde),
            0
        )
    } else {
        rownames(scaled.mat) <- gnams
        blups <- matrix(0, nrow = length(state), ncol = ntrait,
                        dimnames = list(names(state), trait.levels))
        blups[as.logical(state), ] <- scaled.mat
    }
    oint <- tint
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
# Boundary LRT comparing qtlModel (with composite vm/mbf term) against
# baseModel (same fixed QTL effects and residual structure, no vm/mbf term).
#
# ntrait == 1 : one-sided test of a scalar variance parameter at its lower
#   boundary.  Null: 1/2*chi^2(0) + 1/2*chi^2(1).
#
# ntrait > 1  : boundary test of the ntrait x ntrait genetic covariance matrix
#   Ga of the vm term being zero.  Null: pchisq.mixture(stat, ntrait).
#
#   Justification: under H0 (Ga = 0), only the ntrait diagonal variance
#   parameters are at their boundary (>= 0).  Correlation and factor-loading
#   parameters are non-identifiable at Ga = 0 and do not contribute degrees
#   of freedom to the null distribution (Self & Liang 1987, non-regular case).
#   Consequently pchisq.mixture(stat, ntrait) is the exact null distribution
#   for diag(Trait):vm() and an asymptotically valid approximation for
#   corgh(Trait):vm() and fa(Trait,k):vm().
#
#   For ntrait > 1 the caller is responsible for passing a qtlModel whose
#   vm/mbf term uses diag(Trait): (built via .buildDiagVmModel()).  This
#   removes correlation/loading parameters that are non-identifiable under H0
#   asymptotically but weakly identifiable in finite samples, which would
#   otherwise inflate the LRT statistic and cause anti-conservative tests.
#   baseModel (no vm term, same residual structure and QTL fixed effects)
#   is the correct null and requires no modification.
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
