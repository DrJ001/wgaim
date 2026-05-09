# =============================================================================
# engine_select.R
# Internal functions for interval/marker selection and significance testing.
# =============================================================================

#' @keywords internal
.qtlSelect <- function(asm, phenoData, intervalObj, gen.type, selection,
                       exclusion.window, state, verboseLev) {
    sigma2 <- asm$sigma2
    if (asm$vparameters.con[length(asm$vparameters.con)] == 4)
        sigma2 <- 1
    if (!is.null(cov.env <- attr(intervalObj, "env"))) {
        cat(" Predict step for outlier statistics \n")
        cat("=====================================\n")
        rterms <- attr(terms.formula(asm$call$random), "term.labels")
        vmterm <- rterms[grep("vm.*covObj", rterms)]
        pv <- predict(asm, classify = vmterm, only = vmterm, vcov = TRUE, data = phenoData)
        avar <- asm$vparameters[grep("vm.*covObj", names(asm$vparameters))] * sigma2
        atilde <- pv$pvals[, "predicted.value"]
        qtilde <- as.vector(cov.env$trans %*% atilde)
        vatilde <- avar * cov.env$relm - as.matrix(pv$vcov)
        qhalf <- cov.env$trans %*% vatilde
        vqtilde <- colSums(t(qhalf) * t(cov.env$trans))
    } else {
        avar <- asm$vparameters[grep("mbf.*ints", names(asm$vparameters))] * sigma2
        mbf <- grep("mbf", rownames(asm$coefficients$random))
        qtilde <- asm$coefficients$random[mbf, 1]
        pevar <- sigma2 * asm$vcoeff$random[mbf]
        vqtilde <- avar - pevar
    }
    gnams <- names(state)[as.logical(state)]
    names(qtilde) <- names(vqtilde) <- gnams
    oint <- ifelse(!is.na(qtilde^2 / vqtilde), qtilde^2 / vqtilde, 0)
    names(oint) <- gnams
    ochr <- NULL
    if (selection == "chromosome") {
        chr.names <- names(intervalObj$geno)
        nochr <- length(chr.names)
        allc <- sapply(strsplit(gnams, "\\."), "[", 2)
        ochr <- c()
        for (c in 1:nochr) {
            whc <- allc %in% chr.names[c]
            nums <- qtilde[whc]^2
            dens <- vqtilde[whc]
            ochr[c] <- ifelse(!is.na(sum(nums) / sum(dens)), sum(nums) / sum(dens), 0)
        }
        names(ochr) <- chr.names
        mchr <- chr.names[ochr == max(ochr)]
        cint <- allc %in% mchr
        chri <- oint[cint]
        mint <- (1:length(chri))[chri == max(chri)]
        qtl <- names(chri)[mint]
        if (verboseLev > 0) {
            cgen <- ifelse(gen.type == "marker", "Marker", "Interval")
            cat("\n Selection of chromosome using the AOM statistic\n")
            cat("=============================================== \n")
            for (i in 1:nochr)
                cat(" Chromosome ", chr.names[i], "Outlier Statistic ", ochr[i], "\n")
            cat("============================================= \n\n")
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for (i in 1:length(chri))
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
            for (i in 1:length(oint))
                cat(cgen, names(oint)[i], "Outlier Statistic ", oint[i], "\n")
            cat("=============================================== \n\n")
        }
    }
    # Take first if tied, fill out full state vectors, apply exclusion window
    qtl <- qtl[1]
    blups <- tint <- state
    tint[as.logical(state)] <- oint
    blups[as.logical(state)] <- qtilde / sqrt(abs(vqtilde))
    oint <- tint
    schr <- sapply(strsplit(names(state), "\\."), "[", 2)
    wnams <- names(state)[schr %in% mchr]
    inums <- as.numeric(sapply(strsplit(wnams, "\\."), "[", 3))
    if (gen.type == "marker")
        dists <- intervalObj$geno[[mchr]]$map
    else
        dists <- intervalObj$geno[[mchr]]$inferred.map
    dists <- dists[inums]
    exc <- wnams[abs(dists - dists[mint]) <= exclusion.window]
    state[exc] <- 0
    list(state = state, qtl = qtl, ochr = ochr, oint = oint, blups = blups)
}

#' @keywords internal
.lrtTest <- function(qtlModel, baseModel, TypeI) {
    baseLogL <- baseModel$loglik
    stat <- 2 * (qtlModel$loglik - baseLogL)
    pvalue <- (1 - pchisq(stat, 1)) / 2
    pass <- stat >= qchisq(1 - 2 * TypeI, 1)
    list(baseLogL = baseLogL, stat = stat, pvalue = pvalue, pass = pass)
}
