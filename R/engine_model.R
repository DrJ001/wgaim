# =============================================================================
# engine_model.R
# Internal functions for genome-wide model construction and management.
# Includes: relationship matrix, model building, covObj rebuilding,
#           variance parameter reset, and environment cleanup.
# =============================================================================

#' @keywords internal
.constructCM <- function(genoData, scale.method = "diag") {
    tg <- t(genoData)
    relm <- crossprod(tg)
    scale <- mean(diag(relm))
    relm <- relm / scale
    sv <- eigen(relm)
    if (any(sv$values < 0)) {
        sv$values[sv$values < 0] <- 0
        eps <- 1e-12 * scale
        relm <- sv$vectors %*% diag(sv$values) %*% t(sv$vectors) + eps * diag(nrow(relm))
    }
    attr(relm, "rowNames") <- dimnames(relm)[[1]] <- dimnames(relm)[[2]] <- rownames(genoData)
    ch <- chol(relm)
    chol.inv <- chol2inv(ch)
    rm.env <- new.env()
    rm.env$trans <- (tg %*% chol.inv) / scale
    rm.env$relm <- relm
    rm.env$scale <- scale
    rm.env
}

#' @keywords internal
.buildGenomeModel <- function(baseModel, genoData, phenoData, merge.by,
                               intervalObj, force, rterms, caller.env, ...) {
    qtlModel <- baseModel
    vm <- FALSE
    vmterms <- NULL
    cov.env <- NULL
    if ((ncol(genoData) > nrow(genoData)) & !force) {
        cov.env <- .constructCM(genoData)
        covObj <- cov.env$relm
        vmterms <- c(paste0("vm(", merge.by, ", covObj)"), merge.by)
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
    assign("covObj", covObj, envir = caller.env)
    cat("\nRandom Effects Interval/Marker Model Iteration (1):\n")
    cat("============================================\n")
    qtlModel$call$data <- quote(phenoData)
    qtlModel <- update(qtlModel, random. = ran.form, ...)
    list(qtlModel = qtlModel, intervalObj = intervalObj,
         cov.env = cov.env, vm = vm, vmterms = vmterms)
}

#' @keywords internal
.rebuildCovObj <- function(genoData, state, merge.by, intervalObj,
                            force, vm, vmterms, qtlModel, caller.env) {
    mout <- (1:ncol(genoData))[!as.logical(state)]
    genoSub <- genoData[, -mout, drop = FALSE]
    cov.env <- NULL
    if ((ncol(genoSub) > nrow(genoSub)) & !force) {
        cov.env <- .constructCM(genoSub)
        covObj <- cov.env$relm
        attr(intervalObj, "env") <- cov.env
    } else {
        covObj <- cbind.data.frame(rownames(genoSub), genoSub)
        names(covObj)[1] <- merge.by
        if (is.null(qtlModel$call$mbf$ints) & vm) {
            attr(intervalObj, "env") <- NULL
            rterms <- unlist(strsplit(deparse(qtlModel$call$random[[2]]), " \\+ "))
            rterms <- rterms[!(rterms %in% vmterms)]
            qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
            qtlModel$call$mbf$ints$cov <- "covObj"
            ran.form <- as.formula(paste(c("~ mbf('ints')", merge.by, rterms), collapse = " + "))
            qtlModel$call$random <- ran.form
        }
    }
    assign("covObj", covObj, envir = caller.env)
    list(cov.env = cov.env, intervalObj = intervalObj, qtlModel = qtlModel)
}

#' @keywords internal
.vModify <- function(model, merge.by) {
    namg <- names(model$G.param)
    terms <- paste("mbf*.ints*", paste0("vm\\(", merge.by, "*"), "X\\.", sep = "|")
    nterm <- grep(terms, namg)
    if (length(nterm)) {
        for (i in nterm) {
            con.term <- model$G.param[[i]][[1]]$con == "B"
            if (any(con.term)) {
                model$G.param[[i]][[1]]$con[con.term] <- "P"
                model$G.param[[i]][[1]]$initial[con.term] <- 0.1
            }
        }
    }
    model
}

#' @keywords internal
.envFix <- function(model, asremlEnv) {
    for (i in names(asremlEnv)) {
        attr(model$formulae[[i]], ".Environment") <- asremlEnv[[i]]
        if (i %in% names(model$call))
            environment(model$call[[i]]) <- asremlEnv[[i]]
    }
    for (i in names(attributes(model$mf)$model.terms))
        attr(attributes(model$mf)$model.terms[[i]]$Terms.obj, ".Environment") <- NULL
    attributes(model$mf)$mbf.env <- attributes(model$mf)$points.env <- NULL
    model
}
