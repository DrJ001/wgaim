wgaim <- function (baseModel, ...)
UseMethod("wgaim")

wgaim.default <- function(baseModel, ...)
  stop("Currently the only supported method is \"asreml\"")

wgaim.asreml <- function (baseModel, intervalObj, merge.by = NULL, fix.lines = TRUE, gen.type = "interval", method = "fixed", selection = "interval", force = FALSE, exclusion.window = 20, breakout = -1, TypeI = 0.05, trace = TRUE, verboseLev = 0, ...)
{
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
        if(!baseModel$converge)
            stop("Base model not converged: Check base model before proceeding with QTL analysis.")
    }
    asremlEnv <- lapply(baseModel$formulae, function(el) attr(el, ".Environment"))
    phenoData <- eval(baseModel$call$data)
    if (missing(phenoData))
        stop("phenoData is a required argument.")
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj is not of class \"interval\"")
    if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
    if (is.null(glines <- intervalObj$pheno[, merge.by]))
        stop("Genotypic data does not contain column \"", merge.by,
             "\".")
    if (is.null(plines <- phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
             "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in genotypic \"", merge.by, "\" column do not match any names in phenotypic \"",
             merge.by, "\" column.")
    if (!(method %in% c("fixed","random")))
        stop("Method has to be either \"fixed\" or \"random\" (see ?wgaim.asreml).")
    if (!(selection %in% c("interval","chromosome")))
        stop("Selection method has to be either \"interval\" or \"chromosome\" (see ?wgaim.asreml).")
    if(!is.numeric(breakout) | breakout < -1 | breakout == 0)
        stop("breakout argument must be -1 or a positive integer.")
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }
    if(gen.type %in% "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$interval.data)
    else gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    nint <- lapply(gdat, function(el) 1:ncol(el))
    lint <- unlist(lapply(nint, length))
    mnams <- paste("Chr", rep(names(intervalObj$geno), times = lint), unlist(nint), sep = ".")
    dimnames(genoData) <- list(as.character(glines), mnams)
    genoData <- genoData[rownames(genoData) %in% as.character(plines),]
    rterms <- unlist(strsplit(deparse(baseModel$call$random[[2]]), " \\+ "))
    rterms <- rterms[-grep(merge.by, rterms)]
    whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)
    genetic.term <- merge.by
    vm <- FALSE
    if(!all(whg) & fix.lines){
#        !all(whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)))
        phenoData$Gomit <- phenoData$Gsave <- plines
        levels(phenoData$Gsave)[!whg] <- NA
        levels(phenoData$Gomit)[whg] <- "GEN"
        fix.form <- as.formula(paste(". ~ Gomit + .", sep = ""))
        ran.base <- formula(paste(c("~ Gsave", rterms), collapse = " + "))
        baseModel$call$data <- quote(phenoData)
        cat("\nFixing lines and updating initial base model:\n")
        cat("============================================\n")
        baseModel <- update(baseModel, fixed. = fix.form, random. = ran.base, ...)
        merge.by <- "Gsave"
    }
    qtlModel <- baseModel
    if((ncol(genoData) > nrow(genoData)) & !force){
        cov.env <- constructCM(genoData)
        covObj <- cov.env$relm
        vmterms <- c(paste("vm","(",merge.by,", covObj)", sep = ""), merge.by)
        ran.form <- as.formula(paste(c("~", vmterms, rterms), collapse = " + "))
        attr(intervalObj, "env") <- cov.env
        vm <- TRUE
    } else {
        covObj <- cbind.data.frame(rownames(genoData), genoData)
        names(covObj)[1] <- merge.by
        qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
        qtlModel$call$mbf$ints$cov <- "covObj"
        ran.form <- as.formula(paste(c("~ mbf('ints')", merge.by, rterms), collapse = " + "))
    }
    assign("covObj", covObj, envir = parent.frame())
    cat("\nRandom Effects Interval/Marker Model Iteration (1):\n")
    cat("============================================\n")
    qtlModel$call$data <- quote(phenoData)
    qtlModel <- update(qtlModel, random. = ran.form, ...)
    ldiag <- coef.list <- vcoef.list <- list()
    qtl <- c(); iter <- 1
    state <- rep(1, ncol(genoData))
    names(state) <- mnams
    repeat {
        selq <- qtlSelect(qtlModel, phenoData, intervalObj, gen.type, selection, exclusion.window, state, verboseLev)
        state <- selq$state
        ldiag$oint[[iter]] <- selq$oint
        ldiag$ochr[[iter]] <- selq$ochr
        ldiag$blups[[iter]] <- selq$blups
        baseLogL <- baseModel$loglik
        stat <- 2 * (qtlModel$loglik - baseLogL)
        ldiag$lik[[iter]] <- c(baseLogL, qtlModel$loglik, stat, (1 - pchisq(stat, 1))/2)
        if ((stat < qchisq(1 - 2 * TypeI, 1)) | (breakout == iter))
            break
        qtl[iter] <- selq$qtl
        cqtl <- strsplit(qtl[iter], "\\.")
        wchr <- sapply(cqtl, "[", 2)
        wint <- sapply(cqtl, "[", 3)
        message("Found QTL on chromosome ", wchr, " ", gen.type, " ", wint)
        tmp <- cbind.data.frame(rownames(genoData), genoData[, qtl[iter]])
        qtl.x <- gsub("Chr\\.", "X.", qtl[iter])
        names(tmp) <- c(merge.by, qtl.x)
        phenoData <- cbind.data.frame(ord = 1:nrow(phenoData), phenoData)
        phenoData <- merge(phenoData, tmp, by = merge.by, all.x = TRUE, all.y = FALSE)
        phenoData <- phenoData[order(phenoData$ord), ]
        phenoData <- phenoData[, -2]
        mout <- (1:ncol(genoData))[!as.logical(state)]
        genoSub <- genoData[, -mout]
        if((ncol(genoSub) > nrow(genoSub)) & !force){
            cov.env <- constructCM(genoSub)
            covObj <- cov.env$relm
            attr(intervalObj, "env") <- cov.env
        }
        else {
            covObj <- cbind.data.frame(rownames(genoSub), genoSub)
            names(covObj)[1] <- merge.by
            if(is.null(qtlModel$call$mbf$ints) & vm){
                attr(intervalObj, "env") <- NULL
                rterms <- unlist(strsplit(deparse(qtlModel$call$random[[2]]), " \\+ "))
                rterms <- rterms[!(rterms %in% vmterms)]
                qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
                qtlModel$call$mbf$ints$cov <- "covObj"
                ran.form <- as.formula(paste(c("~ mbf('ints')", merge.by, rterms), collapse = " + "))
                qtlModel$call$random <- ran.form
            }
        }
        assign("covObj", covObj, envir = parent.frame())
        qtlModel$call$data <- baseModel$call$data <- quote(phenoData)
        if(method == "random"){
            ran.form <- formula(paste("~ . +", qtl.x, sep = ""))
            cat("\nRandom Effects QTL Model Iteration (", iter, "):\n")
            cat("========================================\n")
            baseModel <- vModify(baseModel, merge.by)
            baseModel <- update(baseModel, random. = ran.form, ...)
            cat("\nRandom Effects QTL plus Interval/Marker Model Iteration (", iter,"):\n")
            cat("=============================================================\n")
            qtlModel <- vModify(qtlModel, merge.by)
            qtlModel <- update(qtlModel, random. = ran.form, ...)
            list.coefs <- qtlModel$coefficients$random
            zind <- grep("X\\.", rownames(list.coefs))
            sub.list <- list.coefs[zind, 1]
            names(sub.list) <- rownames(list.coefs)[zind]
            coef.list[[iter]] <- sub.list
            vcoef.list[[iter]] <- qtlModel$vcoeff$random[zind]
        }
        else {
            fix.form <- as.formula(paste(". ~ . +", qtl.x, sep = ""))
            cat("\nFixed Effects QTL Model Iteration (", iter, "):\n")
            cat("========================================\n")
            baseModel <- update(baseModel, fixed. = fix.form, ...)
            cat("\nFixed Effects QTL plus Interval/Marker Model Iteration (", iter, "):\n")
            cat("============================================================\n")
            qtlModel <- update(qtlModel, fixed. = fix.form, ...)
            list.coefs <- qtlModel$coefficients$fixed
            zind <- grep("X\\.", rownames(list.coefs))
            sub.list <- rev(list.coefs[zind, 1])
            names(sub.list) <- rev(rownames(list.coefs)[zind])
            coef.list[[iter]] <- sub.list
            vcoef.list[[iter]] <- rev(qtlModel$vcoeff$fixed[zind])
        }
        iter <- iter + 1
    }
    qtl.list <- list()
    qtl.list$selection <- selection
    qtl.list$method <- method
    qtl.list$type <- gen.type
    qtl.list$diag <- ldiag
    qtl.list$iterations <- iter
    if (length(qtl)) {
        qtl.list$diag$coef.list <- coef.list
        qtl.list$diag$vcoef.list <- vcoef.list
        qtl.list$diag$lik.mat <- matrix(unlist(ldiag$lik), ncol = 4, byrow = TRUE)
        dimnames(qtl.list$diag$lik.mat)[[2]] <- c("L0","L1","Statistic","Pvalue")
        qtl.list$diag$state <- state
        qtl.list$diag$genetic.term <- genetic.term
        qtl.list$diag$rel.scale <- 1
        if(exists("cov.env")) qtl.list$diag$rel.scale <- cov.env$scale
        qtl.list$breakout <- ifelse(breakout != -1, TRUE, FALSE)
        qtl.list$qtl <- qtl
        qtl.list$effects <- coef.list[[iter - 1]]
        qtl.list$veffects <- vcoef.list[[iter - 1]]
    }
    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = parent.frame())
    qtlModel <- envFix(qtlModel, asremlEnv)
    ## baseModel$call$data <- as.name(data.name)
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("wgaim", "asreml")
    qtlModel
}

vModify <- function (model, merge.by){
    namg <- names(model$G.param)
    terms <- paste("mbf*.ints*", paste("vm\\(", merge.by, "*", sep = ""), "X\\.", sep = "|")
    nterm <- grep(terms, namg)
    if(length(nterm)){
        for (i in nterm){
            con.term <- model$G.param[[i]][[1]]$con == "B"
            if (any(con.term)) {
                model$G.param[[i]][[1]]$con[con.term] <- "P"
                model$G.param[[i]][[1]]$initial[con.term] <- 0.1
            }
        }
    }
    model
}

envFix <- function(model, asremlEnv){
    for(i in names(asremlEnv)){
        attr(model$formulae[[i]], ".Environment") <- asremlEnv[[i]]
        if(i %in% names(model$call))
            environment(model$call[[i]]) <- asremlEnv[[i]]
    }
    for(i in names(attributes(model$mf)$model.terms))
        attr(attributes(model$mf)$model.terms[[i]]$Terms.obj, ".Environment") <- NULL
    attributes(model$mf)$mbf.env <- attributes(model$mf)$points.env <- NULL
    model
}

constructCM <- function(genoData, scale.method = "diag") {
        tg <- t(genoData)
        relm <- crossprod(tg)
        scale <- mean(diag(relm))
        relm <- relm/scale
        attr(relm, "rowNames") <- dimnames(relm)[[2]] <- rownames(genoData)
        ch <- chol(relm)
        chol.inv <- chol2inv(ch)
        rm.env <- new.env()
        rm.env$trans <- (tg %*% chol.inv)/scale
        rm.env$relm <- relm
        rm.env$scale <- scale
        rm.env
}

getQTL <- function (object, intervalObj)
{
  spe <- strsplit(names(object$QTL$effects),"\\.")
  wchr <- sapply(spe, "[", 2)
  wint <- as.numeric(sapply(spe, "[", 3))
  qtlm <- matrix(ncol = 6, nrow = length(wchr))
  for (i in 1:length(wchr)) {
    lhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
    qtlm[i, 1:4] <-  c(wchr[i], wint[i], names(lhmark), round(lhmark, 2))
    if(object$QTL$type == "interval"){
        if(length(intervalObj$geno[[wchr[i]]]$map) > 1)
            rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i] + 1]
        else rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
        qtlm[i, 5:6] <- c(names(rhmark), round(rhmark, 2))
   }
   else qtlm <- qtlm[,-c(5:6)]
 }
  qtlm
}

summary.wgaim <- function (object, intervalObj, LOD = TRUE, ...)
{
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (is.null(qtle <- object$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
        return()
    }
    sigma2 <- object$sigma2
    if(object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1
    if (object$QTL$type == "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$interval.data)
    else gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    gterm <- object$QTL$diag$genetic.term
    scale <- object$QTL$diag$rel.scale
    dimnames(genoData) <- list(as.character(intervalObj$pheno[[gterm]]), names(object$QTL$diag$state))
    genoSub <- genoData[,as.logical(object$QTL$diag$state)]
    if("Gsave" %in% names(object$mf))
        gterm <- "Gsave"
    genoSub <- genoSub[rownames(genoSub) %in% levels(object$mf[[gterm]]),]
    coef.mark <- c(mean(apply(genoSub,1,function(el) sum(el*el)), na.rm=TRUE))
    mark.terms <- paste("mbf*.*ints*", paste("vm\\(", gterm, "*", sep = ""), sep = "|")
    oth.terms <- object$vparameters[-grep(mark.terms, names(object$vparameters))]
    var.mark <- sigma2*object$vparameters[grep(mark.terms, names(object$vparameters))]/scale
    var.res <- sigma2*oth.terms[grep(gterm, names(oth.terms))]
    if(object$QTL$method == "random"){
        var.est <- sigma2*object$vparameters[grep("X\\.", names(object$vparameters))]
        coef.est <- apply(genoData[,object$QTL$qtl, drop = FALSE]^2, 2, mean, na.rm = TRUE)
    } else {
        var.est <- qtle^2
        coef.est <- rep(1, length(qtle))
    }
    var.all <- sum(c(coef.est, coef.mark, 1)*c(var.est, var.mark, var.res))
    perc.var <- round(100*(coef.est*var.est)/var.all, 1)
    zrat <- qtle/sqrt(object$QTL$veffects * sigma2)
    if (object$QTL$method == "random") {
        pvalue <- round((1 - pchisq(zrat^2, df=1))/2, 4)
        pname <- "Prob"
    }
    else {
        pvalue <- round(2 * (1 - pnorm(abs(zrat))), 4)
        pname <- "Pvalue"
    }
    qtlmat <- as.data.frame(matrix(getQTL(object, intervalObj), nrow = length(qtle)))
    qtlmat <- cbind.data.frame(qtlmat[, c(1, 3:ncol(qtlmat))], round(qtle, 4), pvalue, perc.var)
    add.lab <- c("Size", pname, "% Var")
    if (object$QTL$type == "interval")
        collab <- c("Chromosome", "Left Marker", "dist(cM)","Right Marker", "dist(cM)")
    else collab <- c("Chromosome", "Marker", "dist(cM)")
    names(qtlmat) <- c(collab, add.lab)
    if (LOD){
        lod <- round(0.5 * log(exp(zrat^2), base = 10), 4)
        qtlmat <- cbind.data.frame(qtlmat, LOD = lod)
    }
    qtlmat <- qtlmat[order(qtlmat[, 1], as.numeric(as.character(qtlmat[,3]))), ]
    rownames(qtlmat) <- as.character(1:length(qtle))
    qtlmat
}


qtlTable <- function (..., intervalObj = NULL, labels = NULL, columns = "all")
{
    dots <- list(...)
    if (is.null(intervalObj))
        stop("Argument intervalObj cannot be NULL")
    nams <- unlist(lapply(dots, function(el) el$QTL$type))
    if (length(unique(nams)) > 1)
        stop("Models must have been analysed with the same genetic type (see ?wgaim.asreml).")
    mnams <- unlist(lapply(dots, function(el) el$QTL$method))
    if(length(unique(mnams)) > 1)
        stop("Models must have been analysed using the same method (see ?wgaim.asreml).")
    if (!is.null(labels)) {
        if (length(labels) != length(dots))
            stop("Length of labels is not equal to the number of models.")
    }
    else labels <- unlist(lapply(dots, function(el) deparse(el$call$fixed[[2]])))
    olist <- unlist(lapply(dots, function(el) if(!is.null(el$QTL$effects)) TRUE else FALSE))
    dots <- dots[olist]
    labels <- labels[olist]
    qlist <- lapply(dots, function(el, intervalObj, columns) {
        summ <- summary(el, intervalObj)
        if (is.numeric(columns))
            summ <- summ[, columns]
        summ
    }, intervalObj, columns)
    qr <- sapply(qlist, nrow)
    qt <- do.call("rbind.data.frame", qlist)
    qt <- cbind.data.frame(Trait = rep(labels, times = qr), qt)
    qt
}

print.wgaim <- function (x, intervalObj, ...)
{
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj is not of class \"interval\"")
    if (is.null(x$QTL$effects))
        cat("There are no significant putative QTL's\n")
    else {
        qtlm <- getQTL(x, intervalObj)
        for (z in 1:nrow(qtlm)) {
            int <- paste(qtlm[z, 1], qtlm[z, 2], sep = ".")
            if(x$QTL$type == "interval")
                cat("\nPutative QTL found on the interval", int,
                    "\nLeft-hand marker is", qtlm[z, 3], "\nRight-hand marker is",
                    qtlm[z, 5], "\n")
            else cat("\nPutative QTL found close to marker", int,
                    "\nMarker is", qtlm[z, 3], "\n")
        }
    }
}

tr <- function(object, ...)
    UseMethod("tr")

tr.wgaim <- function (object, iter = 1:length(object$QTL$effects), lik.out = TRUE, ...)
{
    dots <- list(...)
    if (!is.na(pmatch("digits", names(dots))))
        dig <- dots$digits
    else dig <- options()$digits
    cl <- object$QTL$diag$coef.list
    vl <- object$QTL$diag$vcoef.list
    sigma2 <- object$sigma2
    if(object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1
    zrl <- lapply(1:length(cl), function(i, cl, vl, sigma2)
                  cl[[i]]/(sqrt(vl[[i]]*sigma2)), cl = cl, vl = vl, sigma2 = sigma2)
    if (any(ret <- is.na(pmatch(iter, 1:length(zrl))))) {
        warning("\"iter\" values outside expected range .. using ones that are in iteration range.")
        iter <- iter[!ret]
    }
    if(object$QTL$method == "random")
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- (1 - pchisq(el^2, df=1))/2
            pv <- round(pv, dig)
            c(pv, rep(NA, len - length(pv)))
    }, len = length(zrl), dig = dig)
    else
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- round(2 * (1 - pnorm(abs(el))), dig)
            c(pv, rep(NA, len - length(pv)))
        }, len = length(zrl), dig = dig)
    qtlmat <- do.call("rbind", pvals)
    qnams <- gsub("X\\.", "", names(cl))
    dimnames(qtlmat) <- list(paste("Iter", 1:length(zrl), sep = "."),
        qnams)
    cat("\nIncremental QTL P-value Matrix.\n")
    cat("===============================\n")
    qtlmat[qtlmat < 0.001] <- "<0.001"
    qtlmat[is.na(qtlmat)] <- ""
    print.default(qtlmat[iter, 1:iter[length(iter)]], quote = FALSE,
        right = TRUE, ...)
    if (lik.out) {
        cat("\nLikelihood Ratio Test of Additive Variance Parameter.\n")
        cat("====================================================\n")
        dmat <- round(as.matrix(object$QTL$diag$lik.mat), dig)
        dimnames(dmat)[[1]] <- paste("Iter", 1:(length(zrl) +
            1), sep = ".")
        dmat[, 4][dmat[, 4] < 0.001] <- "<0.001"
        print.default(dmat, quote = FALSE, right = TRUE, ...)
    }
}

qtlSelect <- function(asm, phenoData, intervalObj, gen.type, selection, exclusion.window, state, verboseLev) {
    sigma2 <- asm$sigma2
    if(asm$vparameters.con[length(asm$vparameters.con)] == 4)
        sigma2 <- 1
    if(!is.null(cov.env <- attr(intervalObj, "env"))) {
        cat(" Predict step for outlier statistics \n")
        cat("=====================================\n")
        rterms <- attr(terms.formula(asm$call$random), "term.labels")
        vmterm <- rterms[grep("vm.*covObj", rterms)]
        pv <- predict(asm, classify = vmterm, only = vmterm, vcov=TRUE, data = phenoData)
        avar <- asm$vparameters[grep("vm.*covObj", names(asm$vparameters))]*sigma2
        atilde <- pv$pvals[, 'predicted.value']
        qtilde <- as.vector(cov.env$trans %*% atilde)
        vatilde <- avar * cov.env$relm - as.matrix(pv$vcov)
        qhalf <- cov.env$trans %*% vatilde
        vqtilde <- colSums(t(qhalf)*t(cov.env$trans))
    }
    else {
        avar <- asm$vparameters[grep("mbf.*ints", names(asm$vparameters))]*sigma2
        mbf <- grep("mbf", rownames(asm$coefficients$random))
        qtilde <- asm$coefficients$random[mbf, 1]
        pevar <- sigma2*asm$vcoeff$random[mbf]
        vqtilde <- avar - pevar
    }
    gnams <- names(state)[as.logical(state)]
    names(qtilde) <- names(vqtilde) <- gnams
    oint <- ifelse(!is.na(qtilde^2/vqtilde), qtilde^2/vqtilde, 0)
    names(oint) <- gnams
    ochr <- NULL
    if(selection == "chromosome"){
        chr.names <- names(intervalObj$geno)
        nochr <- length(chr.names)
        allc <- sapply(strsplit(gnams, '\\.'), "[", 2)
        ochr <- c()
        for(c in 1:nochr){
            whc <- allc %in% chr.names[c]
            cqtilde <- qtilde[whc]
            nums <- cqtilde * cqtilde
            dens <- vqtilde[whc]
            ochr[c] <- ifelse(!is.na(sum(nums)/sum(dens)),sum(nums)/sum(dens),0)
        }
        names(ochr) <- chr.names
        mchr <- chr.names[ochr == max(ochr)]
        cint <- allc %in% mchr
        chri <- oint[cint]
        mint <- (1:length(chri))[chri == max(chri)]
        qtl <- names(chri)[mint]
        if(verboseLev > 0) {
            cat("\n Selection of chromosome using the AOM statistic\n")
            cat("=============================================== \n")
            for(i in 1:nochr)
                cat(" Chromosome ", chr.names[i], "Outlier Statistic ", ochr[i], "\n")
            cat("============================================= \n\n")
            cgen <- "Interval"
            if(gen.type == "marker") cgen <- "Marker"
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for(i in 1:length(chri))
                cat(cgen, names(chri)[i], "Outlier Statistic ", chri[i],"\n")
            cat("=============================================== \n\n")
        }
    } else {
        qtl <- names(oint)[oint == max(oint)]
        qsp <- unlist(strsplit(qtl, split="\\."))
        mint <- as.numeric(qsp[3]); mchr <- qsp[2]
        if(verboseLev > 0) {
            cgen <- "Interval"
            if(gen.type == "marker") cgen <- "Marker"
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for(i in 1:length(oint))
                cat(cgen, names(oint)[i], "Outlier Statistic ", oint[i],"\n")
            cat("=============================================== \n\n")
        }
    }

    ## fill out interval stats and update state
    qtl <- qtl[1]
    blups <- tint <- state
    tint[as.logical(state)] <- oint
    blups[as.logical(state)] <- qtilde/sqrt(abs(vqtilde))
    oint <- tint
    ## exclusion window
    schr <- sapply(strsplit(names(state), "\\."), "[", 2)
    wnams <- names(state)[schr %in% mchr]
    inums <- as.numeric(sapply(strsplit(wnams, "\\."),"[", 3))
    dists <- intervalObj$geno[[mchr]]$map
    if((gen.type == "interval") & (length(dists) > 1))
        dists <- dists[2:length(dists)] - diff(dists)/2
    dists <- dists[inums]
    exc <- wnams[abs(dists - dists[mint]) <= exclusion.window]
    state[exc] <- 0
    list(state = state, qtl = qtl, ochr = ochr, oint = oint, blups = blups)
}

cross2int <- function(object, impute = "MartinezCurnow", consensus.mark = TRUE, id = "id",
                      subset = NULL){
  cls <- class(object)[1]
  if(!(cls %in% c("bc","dh","f2","riself")))
    stop("This function is restricted to populations inheriting from classes \"bc\",\"dh\",\"f2\",\"riself\".")
  object <- drop.nullmarkers(object)
  if(!(id %in% names(object$pheno)))
    stop("The unique identifier for the genotypic rows, ", deparse(substitute(id)), ",cannot be found in genotypic data")
  if(!is.null(subset))
    object <- subset(object, chr = subset)
  if(consensus.mark) {
    tpheno <- object$pheno
    object <- fixMap(object, rd = 3)
    object$pheno <- tpheno
  }
  lid <- as.character(object$pheno[[id]])
  mtype <- c("Broman", "MartinezCurnow")
  if(is.na(type <- pmatch(impute, mtype)))
      stop("Missing marker type must be one of \"Broman\" or \"MartinezCurnow\". Partial matching is allowed.")
  impute <- mtype[type]
  if(impute == "Broman")
      object <- argmax.geno(object)
  object$geno <- lapply(object$geno, function(el, impute, lid, cls){
      row.names(el$data) <- as.character(lid)
      if(!is.null(el$argmax))
          el$imputed.data <- el$argmax
      else el$imputed.data <- el$data
      if(cls %in% "f2"){
          el$imputed.data[el$imputed.data == 3] <- -1
          el$imputed.data[el$imputed.data == 2] <- 0
      } else
          el$imputed.data[el$imputed.data == 2] <- - 1
      if(length(el$map) == 1){
          el$dist <- 0; el$theta <- 0; elambda <- 1/2
          names(el$dist) <- names(el$map)
          el$imputed.data[is.na(el$imputed.data)] <- 0
          el$interval.data <- as.matrix(el$imputed.data/2, ncol = 1)
          dimnames(el$interval.data)[[2]] <- names(el$map)
      }
      else {
          el$dist <- diff(el$map)/100
          el$theta <- 0.5*(1-exp(-2*el$dist))
          if(cls %in% "riself")
              el$theta <- (el$theta/2)/(1 - el$theta)
          elambda <- el$theta/(2*el$dist*(1-el$theta))
          if(impute == "MartinezCurnow")
              el$imputed.data <- imputeGen(el$theta, el$imputed.data, dom = FALSE)$add
          dimnames(el$imputed.data)[[1]] <- dimnames(el$data)[[1]]
          lambda <- addiag(elambda,-1) + addiag(c(elambda,0),0)
          lambda <- lambda[,-dim(lambda)[2]]
          el$interval.data <- el$imputed.data %*% lambda
          dimnames(el$interval.data)[[2]] <- names(el$dist)
      }
      el
  }, impute, lid, cls)
  if(length(grep("\\.", names(object$geno)))){
      warning("Removing \".\" from linkage group names.")
      names(object$geno) <- gsub("\\.", "", names(object$geno))
  }
  class(object) <- c(class(object), "interval")
  object
}

fixMap <- function (full.data, rd = 3)
{
    drop.mark <- lapply(full.data$geno, function(el) {
        emap <- round(el$map, rd)
        um <- unique(emap)
        dmark <- clist <- NA
        if(length(um) != length(emap)){
            pm <-  pmatch(emap, um, duplicates.ok = TRUE)
            pmt <- table(pm)
            nums <- as.numeric(names(table(pm))[pmt > 1])
            pm[!(pm %in% nums)] <- nums[length(nums)] + 1
            clist <- split(emap, pm)
            if(any(pmt == 1)) len <- length(clist) - 1
            else len <- length(clist)
            combl <- lapply(clist[1:len], function(cl){
                names(cl)[1] <- paste(names(cl)[1], "(C)", sep = "")
                cbind.data.frame(marker = names(cl), dist = cl)
                })
            clist <- do.call("rbind", combl)
            clist$bin <- rep(1:length(combl), times = sapply(combl, nrow))
            dlist <- split.data.frame(t(el$data), pm)
            dmark <- lapply(dlist[1:len], function(dl) {
                con <- apply(dl, 2, function(ell){
                    ell <- ell[!is.na(ell)]
                    if(length(ellu <- unique(ell)) > 1 | !length(ellu))
                       NA else ellu
                })
                dn <- dimnames(dl)[[1]]
                con <- matrix(con, nrow = 1, dimnames = list(dn[1]))
                dm <- dn[-1]
                list(con = con, dm = dm)
               })
            cond <- t(do.call("rbind", lapply(dmark, function(el) el$con)))
            dind <- pmatch(dimnames(cond)[[2]], dimnames(el$data)[[2]])
            cnam <- paste(dimnames(cond)[[2]], "(C)", sep = "")
            el$data[,dind] <- cond
            dimnames(el$data)[[2]][dind] <- names(el$map)[dind] <- cnam
            dmark <- unlist(lapply(dmark, function(el) el$dm))
        }
        list(chr = el, dmark = dmark, clist = clist)
    })
    full.data$geno <- lapply(drop.mark, function(dm) dm$chr)
    dmarkl <- unlist(lapply(drop.mark, function(dm) dm$dmark))
    newmap <- drop.markers(full.data, dmarkl[!is.na(dmarkl)])
    chre <- unlist(lapply(drop.mark, function(cm) any(is.na(cm))))
    cor.mark <- as.data.frame(do.call("rbind", lapply(drop.mark[!chre], function(cm) cm$clist)))
    chrn <- unlist(lapply(drop.mark[!chre], function(cm) nrow(cm$clist)))
    cor.mark$chr <- rep(names(nmar(full.data))[!chre], times = chrn)
    cor.mark$bin <- paste(cor.mark$chr, cor.mark$bin, sep = ".")
    rownames(cor.mark) <- NULL
    newmap$colocated.markers <- cor.mark
    newmap
}



addiag <- function(x = 1, di = 0, nrow.arg, ncol.arg = n)
{
    if(is.matrix(x)) {
        k <- ifelse(col(x) == (row(x) + di), TRUE, FALSE)
        return(x[k])
    }
    if(missing(x))
        n <- nrow.arg
    else if(length(x) == 1 && missing(di) && missing(nrow.arg) && missing(ncol.arg)) {
        n <- as.integer(x)
        x <- 1
    }
    else n <- length(x)
    if(!missing(nrow.arg))
        n <- nrow.arg
    k <- abs(di)
    p <- ncol.arg + k
    n <- n + k
    m <- matrix(0, n, p)
    k <- ifelse(col(m) == (row(m) + di), TRUE, FALSE)
    m[k] <- x
    m
}

imputeGen <-  function (theta, chr, dom = TRUE){
    th.f <- function(the, ind){
        th <- 0
        for(i in ind)
            th <- th + the[i - 1] - 2*th*the[i - 1]
        th
    }
    dom.gen <- function(thL, thR, thLR, wh){
        switch(wh, a = (2*thL*(1 - thL)*thR*(1 - thR))/(1 - thLR)^2,
               b = (2*thL*(1 - thL)*thR*(1 - thR))/(thLR)^2,
               c = (thL*(1- thL)*(1 - 2*thR*(1 - thR)))/(thLR*(1 - thLR)),
               d = (thR*(1- thR)*(1 - 2*thL*(1 - thL)))/(thLR*(1 - thLR)),
               e = ((thL^2 + (1- thL)^2)*(thR^2 + (1 - thR)^2))/(thLR^2 + (1 - thLR)^2)
               )}
    chrd <- NULL
    if(dom){
        chrd <- chr + 1
        chrd[chrd %in% 2] <- 0
    }
    wh <- which(is.na(chr), arr.ind = TRUE)
    if(dim(wh)[1] != 0){
        wh <- wh[order(wh[,"row"], wh[,"col"]),,drop = FALSE]
        sp <- split(wh[,"col"], wh[,"row"])
        lr <- lapply(sp, function(el, n){
            left <- el - 1; right <- el + 1
            while(any(c(left,right) %in% el)){
                left[left %in% el] <- left[left %in% el] - 1
                right[right %in% el] <- right[right %in% el] + 1
            }
            left[left == 0] <- right[right == n+1] <- NA
            list(left,right)
        }, n = ncol(chr))
        left <- unlist(sapply(lr, "[", 1))
        right <- unlist(sapply(lr, "[", 2))
        xL <- chr[cbind(wh[,"row"], left)]
        xR <- chr[cbind(wh[,"row"], right)]
        whc <- wh[,"col"]
        filld <- filla <- c()
        for(i in 1:nrow(wh)){
            if(is.na(left[i]) & is.na(right[i]))
                filld[i] <- filla[i] <- 0
            else if(!is.na(left[i]) & is.na(right[i])){
                tl <- th.f(theta, whc[i]:(left[i] + 1))
                filla[i] <- xL[i]*(1 - 2*tl)
                if(dom){
                    if(abs(xL[i]) == 1)
                        filld[i] <- 2*tl*(1 - tl)
                    else filld[i] <- 1 - 2*tl*(1 - tl)
                }
            } else if(is.na(left[i]) & !is.na(right[i])){
                tr <- th.f(theta, (whc[i] + 1):(right[i]))
                filla[i] <- xR[i]*(1 - 2*tr)
                if(dom){
                    if(abs(xR[i]) == 1)
                        filld[i] <- 2*tr*(1 - tr)
                    else filld[i] <- 1 - 2*tr*(1 - tr)
                }
            } else {
                tl <- th.f(theta, whc[i]:(left[i] + 1))
                tr <- th.f(theta, (whc[i] + 1):right[i])
                tlr <- tl + tr - 2*tl*tr
                if(tlr == 0)
                    filla[i] <- xL[i]
                else {
                    lambda <- (tr*(1 - tr)*(1 - 2*tl))/(tlr*(1 - tlr))
                    rho <- (tl*(1 - tl)*(1 - 2*tr))/(tlr*(1 - tlr))
                    filla[i] <- xL[i]*lambda + xR[i]*rho
                }
                if(dom){
                    if(abs(xL[i]*xR[i]) > 0){
                        if(xL[i] == xR[i])
                            filld[i] <- dom.gen(tl, tr, tlr, "a")
                        else filld[i] <- dom.gen(tl, tr, tlr, "b")
                    }
                    else {
                        if(abs(xL[i]) == 1)
                            filld[i] <- dom.gen(tl, tr, tlr, "c")
                        else {
                            if(abs(xR[i]) == 1)
                                filld[i] <- dom.gen(tl, tr, tlr, "d")
                            else filld[i] <- dom.gen(tl, tr, tlr, "e")
                        }
                    }
                }
            }
        }
        chr[wh] <- filla
        if(dom)
            chrd[wh] <- filld
    }
    list(add = chr, dom = chrd)
}

linkMap <- function(object, ...)
    UseMethod("linkMap")

linkMap.cross <- function(object, chr, chr.dist, marker.names = "markers", tick = FALSE, squash = TRUE, m.cex = 0.6, ...){
  circ <- function(x, y, shiftx = 0, shifty = 0, ely = 1, elx = 1)
    ((x - shiftx)^2)/elx + ((y - shifty)^2)/ely
  dots <- list(...)
  old.xpd <- par("xpd")
  par(xpd = TRUE)
  on.exit(par(xpd = old.xpd))
  map <- pull.map(object)
  if(!missing(chr)) {
    if(any(is.na(pmatch(chr, names(map)))))
      stop("Some names of chromosome(s) subset do not match names of map.")
    map <- map[chr]
  }
  n.chr <- length(map)
  mt <- list()
  if(!missing(chr.dist)) {
    if(!all(names(chr.dist) %in% c("start","end")))
        stop("names of chr.dist must be \"start\" and/or \"end\"")
    dname <- names(chr.dist)
     if(!is.null(chr.dist$start)) {
        if(length(chr.dist$start) == 1)
            chr.dist$start <- rep(chr.dist$start, length(map))
        if(length(chr.dist$start) != length(map))
            stop("Length of user specified chromosome starting distances need to be 1 or equal to the number of chromosomes")
    } else chr.dist$start <- 0
    if(!is.null(chr.dist$end)) {
        if(length(chr.dist$end) == 1)
            chr.dist$end <- rep(chr.dist$end, length(map))
        if(length(chr.dist$end) != length(map))
            stop("Length of user specified chromosome ending distances need to be 1 or equal to the number of chromosomes")
    } else chr.dist$end <- unlist(lapply(map, max))
    chr.dist <- do.call("cbind.data.frame", chr.dist)
    nm <- names(map)
    map <- lapply(1:nrow(chr.dist), function(el, map, chr.dist){
        tmap <- map[[el]]
        if(max(tmap) < chr.dist$start[el]) NULL
        else tmap[(tmap >= chr.dist$start[el]) & (tmap <= chr.dist$end[el])]
    }, map, chr.dist)
    names(map) <- nm
    map <- map[!sapply(map, is.null)]
    n.chr <- length(map)
  }
  maxlen <- max(unlist(lapply(map, max)))
  minlen <- min(unlist(lapply(map, min)))
  omap <- map
  if(is.null(marker.names)) {
    chrpos <- 1:n.chr
    thelim <- range(chrpos) + c(-0.5, 0.5)
  }
  else {
    if(!is.null(fmark <- attr(object, "flanking"))){
        map <- lapply(map, function(el, fmark){
            el <- el[names(el) %in% fmark]
            if(length(el) == 0) NULL
            else el
        }, fmark)
        omap <- omap[!sapply(map, is.null)]
        map <- map[!sapply(map, is.null)]
        n.chr <- length(map)
    }
    if(all(is.na(pmatch(marker.names, c("markers","dist")))))
      stop("marker.names argument must be either \"dist\", or \"markers\".")
    if(!is.na(pmatch("cex", names(dots))))
      dots$cex <- NULL
#    else cex <- par("cex")
    if(!squash)
      chrpos <- seq(1, n.chr * 3, by = 3)
    else
      chrpos <- seq(1, n.chr * 2, by = 2)
    thelim <- range(chrpos) + c(-1.6, 1.35)
    for(i in 1:n.chr){
        mt[[i]] <- map[[i]]
        if(length(mt[[i]]) > 1){
            conv <- par("pin")[2]/maxlen
            for(j in 1:(length(mt[[i]]) - 1)){
                ch <- mt[[i]][j + 1]*conv - (mt[[i]][j]*conv + 10*par("csi")*m.cex/9)
                if(ch < 0){
                    temp <- mt[[i]][j + 1]*conv + abs(ch)
                    mt[[i]][j + 1] <- temp/conv
                }
            }
        }
    }
    maxlen <- max(c(unlist(lapply(omap, max)),unlist(lapply(mt, max))))
    names(mt) <- names(map)
  }
  plot(0, 0, type = "n", ylim = c(maxlen, minlen), xlim = thelim,
       xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome",
       axes = FALSE, ...)
  axis(side = 2,  ylim = c(maxlen, minlen))
  pins <- par()$plt
  for(i in 1:n.chr) {
      if(!is.null(marker.names)) {
          if(marker.names == "dist")
              alis <- list(x = chrpos[i] + 0.50, y =  mt[[i]], labels = as.character(round(map[[i]], 2)),
                           adj = c(0, 0.5), cex = m.cex)
          else
              alis <- list(x = chrpos[i] + 0.50, y = mt[[i]], labels = names(map[[i]]), adj = c(0, 0.5), cex = m.cex)
          do.call("text", c(alis, dots))
          segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 0.3, map[[i]])
          segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 0.4, mt[[i]])
          segments(chrpos[i] + 0.40, mt[[i]], chrpos[i] + 0.45, mt[[i]])
    }
    map[[i]] <- omap[[i]]
    barl <- chrpos[i]- 0.03
    barr <- chrpos[i]+ 0.03
    segments(barl, min(map[[i]]), barl, max(map[[i]]), lwd = 1)
    segments(barr, min(map[[i]]), barr, max(map[[i]]), lwd = 1)
    segments(barl - 0.17, map[[i]], barr + 0.17, map[[i]])
    # attempt to put curves at ends of chromosomes
    xseq <- seq(barl, barr, length = 20) - chrpos[i]
    yseq <- circ(xseq, xseq, ely = 1, elx = 0.07/maxlen)
    yseq <- yseq - max(yseq)
    lines(xseq + chrpos[i], min(map[[i]]) + yseq)
    lines(xseq + chrpos[i], max(map[[i]]) - yseq)
  }
  axis(side = 1, at = chrpos, labels = names(map), tick = FALSE)
  if(is.na(pmatch("main", names(dots))) & !as.logical(sys.parent()))
    title("Genetic Map")
  invisible(list(mt = mt, map = map, chrpos = chrpos))
}

linkMap.wgaim <- function (object, intervalObj, chr, chr.dist, marker.names = "markers", flanking = TRUE, list.col = list(q.col = "light blue", m.col = "red", t.col = "light blue"), list.cex = list(t.cex = 0.6, m.cex = 0.6), trait.labels = NULL, tick = FALSE, ...)
{
    dots <- list(...)
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (!length(wchr <- object$QTL$effects)) {
        warning("There are no significant QTL's. Plotting map only...")
        linkMap(intervalObj, chr, chr.dist, marker.names = marker.names,
            tick = tick, squash = FALSE, ...)
        return(invisible())
    }
    qtlm <- getQTL(object, intervalObj)
    wchr <- qtlm[,1]
    if(object$QTL$type == "interval")
        qtlm <- qtlm[,3:6]
    else
        qtlm <- cbind(qtlm[,3:4],qtlm[,3:4])
    if (is.null(list.col$q.col))
        list.col$q.col <- "light blue"
    if (is.null(list.col$t.col))
        list.col$t.col <- list.col$q.col
    if (missing(chr))
        chr <- unique(wchr)[order(unique(wchr))]
    if(is.null(list.cex$m.cex))
        list.cex$m.cex <- 0.6
    if(is.null(list.cex$t.cex))
        list.cex$t.cex <- 0.6
    if(flanking)
        attr(intervalObj, "flanking") <- unique(c(qtlm[,1],qtlm[,3]))
    lmap <- linkMap(intervalObj, chr, chr.dist, marker.names = marker.names,
        tick = tick, squash = TRUE, m.cex = list.cex$m.cex, ...)
    map <- lmap$map
    if (is.null(trait <- trait.labels))
        trait <- rep(as.character(object$call$fixed[[2]]), length(wchr))
    if (length(trait.labels) == 1)
        trait <- rep(trait.labels, length(wchr))
    qtrait <- unique(trait)
    qtlm <- cbind.data.frame(qtlm, trait = factor(trait, levels = unique(trait)))
    qtlm[, 2] <- as.numeric(as.character(qtlm[, 2]))
    qtlm[, 4] <- as.numeric(as.character(qtlm[, 4]))
    qtlList <- lapply(split(qtlm, wchr), function(el) el[order(el[,
        2]), ])
    qtlm <- do.call("rbind", qtlList)
    wchr <- wchr[order(wchr)]
    chr <- names(map)
    if (!missing(chr)) {
        oc <- wchr %in% chr
        if(!all(oc)) {
            warning("Some QTL's exist outside chromosome(s) subset, Omitting QTL's....")
            qtlm <- qtlm[oc,]
            wchr <- wchr[oc]
        }
    }
    if (!missing(chr.dist)) {
        mins <- unlist(lapply(map, min))
        maxs <- unlist(lapply(map, max))
        mins <- mins[wchr]; maxs <- maxs[wchr]
        om <- qtlm[, 4] < maxs & qtlm[,2] > mins
        if (!all(om)) {
            warning("Some QTL regions outside distances specified. Omitting QTL's....")
            qtlm <- qtlm[om, ]
            wchr <- wchr[om]
        }
    }
    n.chr <- length(map)
    maxlen <- max(unlist(lapply(map, max)))
    chrpos <- lmap$chrpos
    mt <- lmap$mt
    if (!is.na(cind <- pmatch("col", names(dots))))
        dots <- dots[-cind]
    if (is.null(dim(qtlm)))
        qtlm <- matrix(qtlm, nrow = 1, byrow = FALSE)
    tlis <- list()
    if (!is.na(pmatch("cex", names(dots))))
        p.cex <- dots$cex
    else p.cex <- par("cex")
    dots$cex <- NULL
    for (i in 1:n.chr) {
        conv <- par("pin")[2]/maxlen
        ind <- wchr %in% names(map)[i]
#        if (as.logical(length(ind <- wchr %in% names(map)[i]))) {
         if(any(ind)){
            tlis[[i]] <- as.vector(qtlm[ind, 2] + qtlm[ind, 4])/2
            names(tlis[[i]]) <- as.character(qtlm[ind, 5])
            if (length(tlis[[i]]) > 1) {
                for (j in 1:(length(tlis[[i]]) - 1)) {
                  ch <- tlis[[i]][j + 1] * conv - (tlis[[i]][j] *
                    conv + 10 * par("csi") * list.cex$t.cex/9)
                  if (ch < 0) {
                    temp <- tlis[[i]][j + 1] * conv + abs(ch)
                    tlis[[i]][j + 1] <- temp/conv
                  }
                }
            }
        }
    }
    qtld <- qtlm[, 1:4]
    nodup <- !duplicated(do.call("paste", qtld))
    qtls <- qtld[nodup, ]
    whd <- pmatch(do.call("paste", qtld), do.call("paste", qtls),
        duplicates.ok = TRUE)
    dlis <- split(as.character(qtlm[, 5]), whd)
    qtlm <- qtls
    wchr <- wchr[nodup]
    for (i in 1:n.chr) {
        ind <- wchr %in% names(map)[i]
        if(any(ind)){
            ind <- (1:length(wchr))[ind]
            for (j in ind) {
                if (!is.null(marker.names)) {
                    wh <- mt[[i]][c(as.character(qtlm[j, 1]), as.character(qtlm[j, 3]))]
                  if (marker.names == "dist") {
                    dist <- map[[i]][c(as.character(qtlm[j, 1]), as.character(qtlm[j, 3]))]
                    dlabs <- as.character(round(as.numeric(unlist(dist)), 2))
                    alis <- list(x = chrpos[i] + 0.5, y = wh,
                      labels = dlabs, adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                  }
                  else alis <- list(x = chrpos[i] + 0.5, y = wh,
                    labels = names(wh), adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                  do.call("text", c(alis, dots))
                }
                yv <- c(qtlm[j, 2], qtlm[j, 4])
                yv <- c(yv, rev(yv))
                dind <- dlis[[j]]
                q.cols <- list.col$q.col[pmatch(dind, qtrait)]
                qind <- 1:length(dind)
                if (length(dlis[[j]]) > 1) {
                  int <- seq(chrpos[i] - 0.2, chrpos[i] + 0.2,
                    length = length(dind) + 1)
                  for (k in 1:length(dind)) {
                    xv <- c(rep(int[k], 2), rep(int[k + 1], 2))
                    if(object$QTL$type == "interval")
                        polygon(xv, y = yv, border = NA, col = q.cols[k])
                    else {
                        plis <- list(x = (xv[1] + xv[3])/2,y = yv[1], col = q.cols[k], cex = p.cex)
                        do.call("points", c(plis, dots))
                    }
                  }
                }
                else {
                  xv <- c(rep(chrpos[i] - 0.2, 2), rep(chrpos[i] +
                    0.2, 2))
                  if(object$QTL$type == "interval")
                      polygon(xv, y = yv, border = NA, col = q.cols)
                  else {
                      plis <- list(x = chrpos[i], y = yv[1], col = q.cols, cex = p.cex)
                      do.call("points", c(plis,dots))
                  }
                }
                segments(chrpos[i] - 0.25, yv[1], chrpos[i] -
                  0.25, yv[2])
                segments(chrpos[i] - 0.25, sum(yv[1:2])/2, chrpos[i] -
                  0.3, sum(yv[1:2])/2)
                segments(chrpos[i] - 0.3, sum(yv[1:2])/2, chrpos[i] -
                  0.4, tlis[[i]][qind])
                segments(chrpos[i] - 0.4, tlis[[i]][qind], chrpos[i] -
                  0.45, tlis[[i]][qind])
                if (length(list.col$t.col) > 1)
                  t.cols <- list.col$t.col[pmatch(dind, qtrait)]
                else t.cols <- list.col$t.col
                text(chrpos[i] - 0.5, tlis[[i]][qind], names(tlis[[i]][qind]),
                  adj = c(1, 0.3), col = t.cols, cex = list.cex$t.cex)
                tlis[[i]] <- tlis[[i]][-qind]
            }
        }
        segments(chrpos[i] - 0.2, map[[i]], chrpos[i] + 0.2,
            map[[i]])
    }
    if (is.na(pmatch("main", names(dots))))
        title("Genetic Map with QTL")
}

linkMap.default <- function (object, intervalObj, chr, chr.dist, marker.names = "markers", flanking = TRUE, list.col = list(q.col = rainbow(length(object)), m.col = "red", t.col = rainbow(length(object))),
    list.cex = list(m.cex = 0.6, t.cex = 0.6), trait.labels = NULL, tick = FALSE, ...)
{
    old.par <- par(no.readonly = TRUE)
    par(mar = c(5, 4, 5, 2) + 0.1)
    par(xpd = TRUE)
    on.exit(par(old.par))
    dlist <- list()
    lclass <- lapply(object, function(el) {
        if (!inherits(el, "wgaim"))
            stop("list objects need to inherit from class \"wgaim\"")
        class(el)[1]
    })
    lclass <- unique(unlist(lclass))
    if (any(is.na(pmatch(lclass, "wgaim"))))
        stop("link.map method is only for list objects of class \"wgaim\"")
    type <- unique(unlist(lapply(object, function(el) el$QTL$type)))
    if(length(type) > 1)
        stop("Models need to have same analyses \"type\".")
    dlist$QTL$type <- type
    effects <- lapply(object, function(el) el$QTL$effects)
    len <- lapply(effects, length)
    if (!is.null(trait.labels)) {
        if (length(trait.labels) != length(object))
            stop("Length of trait labels does not equal number of models specified.")
    }
    else trait.labels <- unlist(lapply(object, function(el) as.character(el$call$fixed[[2]])))
    tt <- table(trait.labels)
    tt <- tt[tt > 1]
    if(length(tt)){
        tn <- names(tt)
        for(i in length(tt)) {
            whl <- trait.labels %in% tn[i]
            trait.labels[whl] <- paste(trait.labels[whl], 1:tt[i], sep = "")
        }
    }
    trait <- rep(trait.labels, times = len)
    dlist$QTL$effects <- unlist(effects)
    class(dlist) <- "wgaim"
    if (length(list.col$q.col) != length(object)) {
        warning("QTL colours not of the same length as the number of traits,\n Choosing \"q.col = rainbow(",
            length(object), ")\".")
        list.col$q.col <- rainbow(length(object))
    }
    if (length(list.col$t.col) != length(object)) {
        if (length(list.col$t.col) != 1) {
            warning("Inappropriate length for trait name colours, using QTL colours")
            list.col$t.col <- list.col$q.col
        }
    }
    linkMap(dlist, intervalObj, chr, chr.dist, marker.names = marker.names,
        list.col = list.col, list.cex = list.cex, trait.labels = trait, tick = tick, flanking = flanking, ...)
}

theme_scatter <- function (base_size = 11, base_family = "") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
        legend.position = "none",
        panel.grid.minor = element_line(colour = "grey90", size = 0.4),
        panel.grid.major = element_line(colour = "grey90", size = 0.8),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = base_size),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = base_size),
        strip.text = element_text(size = base_size))
    }

outStat <- function (object, intervalObj, iter = NULL, chr = NULL, statistic = "outlier", plot.chr = FALSE, chr.lines = FALSE)
{
    if (missing(object))
        stop("model object is a required argument.")
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if(plot.chr & is.null(object$QTL$diag$ochr[[1]]))
        stop("There are no chromosome outlier statistics to display.")
    if(plot.chr & !(statistic == "outlier"))
        warning("Ignoring statistic argument and plotting chromosome outlier statistics.")
    if(!is.null(iter)){
        if(any(iter > length(object$QTL$diag$oint)))
            stop("iter argument contains integers greater than number of analysis iterations.")
    } else iter <- 1:length(object$QTL$diag$oint)
    c.iter <- paste("Iteration: ", iter, sep = "")
    if(!is.null(chr)){
        if(!all(chr %in% names(nmar(intervalObj))))
            stop("Some chromosome names do not exist in intervalobj.")
        intervalObj <- subset(intervalObj, chr = chr)
    } else chr <- names(nmar(intervalObj))
    ann.labels <- TRUE
    if(plot.chr){
        ochr <- object$QTL$diag$ochr[iter]
        names(ochr) <- c.iter
        ochrl <- lapply(ochr, function(el, chr){
            values <- el[names(el) %in% chr]
            cbind.data.frame(values = values, chr = names(values))
        }, chr)
        char.iter <- rep(c.iter, times = sapply(ochrl, nrow))
        char.iter <- factor(char.iter, levels = unique(char.iter))
        ochrd <- cbind(do.call("rbind.data.frame", ochrl), iteration = char.iter)
        gp <- ggplot(ochrd, aes_string(x = "chr", y = "values")) +
            facet_wrap( ~ iteration, ncol = 1) + geom_col(fill = "light blue", colour = "grey50") +
            scale_y_continuous(expand = c(0,0)) + ylab("Outlier Statistic") + xlab("Chromosome") +
            theme_scatter()
        ann.labels <- FALSE
    } else {
        if(statistic == "outlier"){
            y.lab <- "Outlier Statistic"
            oint <- object$QTL$diag$oint[iter]
        } else {
            oint <- object$QTL$diag$blups[iter]
            y.lab <- "Scaled BLUPs"
        }
        names(oint) <- c.iter
        ointl <- lapply(oint, function(el, chr){
            echr <- sapply(strsplit(names(el), "\\."), "[", 2)
            whc <- echr %in% chr
            values <- el[whc]
            cbind.data.frame(values = values, chr = echr[whc])
        }, chr)
        char.iter <- rep(c.iter, times = sapply(ointl, nrow))
        char.iter <- factor(char.iter, levels = unique(char.iter))
        distl <- lapply(intervalObj$geno, function(el, type){
            if(length(el$map) == 1) return(0.05)
            if(type %in% "interval")
                el$dist
            else c(0.05, el$dist)
        }, type = object$QTL$type)
        dist <- cumsum(unlist(distl))
        if(object$QTL$type == "interval")
            dist <- dist - unlist(distl)/2
        dist <- rep(dist, length(iter))
        ointd <- cbind(do.call("rbind.data.frame", ointl), dist = dist, iteration = char.iter)
        gp <- ggplot(ointd, aes_string(x = "dist", y = "values", colour = "chr")) +
            facet_wrap( ~ iteration, ncol = 1) + geom_line() + geom_rug(sides = "b", length = unit(0.02, "npc")) +
            ylab(y.lab) + xlab("") + theme_scatter()
        if(chr.lines){
            if(length(chr) > 1){
                ci <- sapply(distl, length)[1:(length(chr) - 1)]
                ci <- c(ci[1], cumsum(ci))
            }
        gp <- gp + geom_vline(xintercept = dist[ci], colour = "grey80", size = 1)
        }
    }
    if(!is.null(qtl <- object$QTL$qtl)){
        iterq <- iter[iter < (length(qtl) + 1)]
        qtls <- gsub("Chr.", "", qtl[iterq])
        qchr <- sapply(strsplit(qtls, "\\."),"[", 1)
        if(any(wchr <- qchr %in% chr)){
            iterq <- iterq[wchr]
            char.iter <- paste("Iteration: ", iterq, sep = "")
            qlabel <- qtls[wchr]
            chrd <- cbind.data.frame(chr = qchr[wchr], qlabel = qlabel, iteration = char.iter)
            valq <- c()
            if(plot.chr){
                for(i in 1:nrow(chrd)){
                    it <- as.character(chrd$iteration[i])
                    valq[i] <- ochr[[it]][names(ochr[[it]]) %in% chrd$chr[i]]
                }
            } else {
                wint <- sapply(strsplit(as.character(chrd$qlabel), "\\."),"[", 2)
                distq <- c()
                cmar <- cumsum(sapply(distl, length))
                for(i in 1:nrow(chrd)){
                    it <- as.character(chrd$iteration[i])
                    valq[i] <- oint[[it]][names(oint[[it]]) %in% qtl[iterq[i]]]
                    cint <- (1:length(cmar))[names(cmar) %in% chrd$chr[i]]
                    cint <- ifelse(cint > 1, cmar[cint - 1], 0)
                    iint <- (i - 1)*cmar[length(cmar)]
                    distq[i] <- dist[iint + cint + as.numeric(wint[i])]
                }
                chrd$dist <- distq
            }
            yrc <- ggplot_build(gp)$layout$panel_scales_y[[1]]$range$range
            ylc <- ifelse(statistic == "outlier", 0, yrc[1] - abs(diff(yrc))/12)
            chrd$values <- valq + sign(valq)*(diff(yrc)/15)
            mx <- max(chrd$values)
            gp <- gp + geom_text(data = chrd, aes_string(label = "qlabel"))
               #
        }
    }
    if(object$QTL$iterations == 1){
            yrc <- ggplot_build(gp)$layout$panel_scales_y[[1]]$range$range
            ylc <- ifelse(statistic == "outlier", 0, yrc[1] + abs(diff(yrc))/40)
    }
    if(ann.labels){
        ann.iter <- paste("Iteration: ", iter[length(iter)], sep = "")
        annd <- ointd[ointd$iteration %in% ann.iter,]
        dist.label <- sapply(split(annd$dist, annd$chr), function(el)
            min(el) + (max(el) - min(el))/2 )
        gp <- gp + scale_x_continuous(breaks = dist.label, labels = chr) + coord_cartesian(ylim = c(ylc, yrc[2] + diff(yrc)/12))
    }
    gp
}

######################## end functions







