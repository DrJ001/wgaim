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
                               intervalObj, force, rterms, caller.env,
                               Trait = NULL, n.fa = 0L, ...) {
    qtlModel <- baseModel
    vm <- FALSE
    vmterms <- NULL
    cov.env <- NULL

    # Build the Trait variance structure suffix for the genome-wide term.
    # ntrait == 1 (Trait = NULL): no suffix -- univariate, behaviour unchanged.
    # ntrait == 2: start with diag(Trait), upgrade to corh(Trait) if n.fa > 0.
    # ntrait >= 3: start with diag(Trait), upgrade to fa(Trait, k) for k = 1..n.fa.
    trait.suffix <- if (!is.null(Trait)) paste0(":diag(", Trait, ")") else ""

    if ((ncol(genoData) > nrow(genoData)) & !force) {
        cov.env <- .constructCM(genoData)
        covObj  <- cov.env$relm
        gterm   <- paste0("vm(", merge.by, ", covObj)", trait.suffix)
        vmterms <- c(gterm, merge.by)
        ran.form <- as.formula(paste(c("~", vmterms, rterms), collapse = " + "))
        attr(intervalObj, "env") <- cov.env
        vm <- TRUE
    } else {
        covObj <- cbind.data.frame(rownames(genoData), genoData)
        names(covObj)[1] <- merge.by
        qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
        qtlModel$call$mbf$ints$cov <- "covObj"
        gterm    <- paste0("mbf('ints')", trait.suffix)
        ran.form <- as.formula(paste(c(paste("~", gterm), merge.by, rterms), collapse = " + "))
    }
    assign("covObj", covObj, envir = caller.env)
    cat("\nRandom Effects Interval/Marker Model Iteration (1):\n")
    cat("============================================\n")
    qtlModel$call$data <- quote(phenoData)
    qtlModel <- update(qtlModel, random. = ran.form, ...)

    # Multivariate only: upgrade variance structure from diag to corh / fa(k)
    if (!is.null(Trait) && n.fa > 0L) {
        ntrait <- length(levels(phenoData[[Trait]]))
        rterms.cur <- attr(terms.formula(qtlModel$call$random), "term.labels")
        gterm.cur  <- rterms.cur[grep("vm.*covObj|mbf.*ints", rterms.cur)]
        if (ntrait == 2L) {
            gterm.new  <- gsub("diag", "corh", gterm.cur)
            ran.form   <- as.formula(paste(
                c("~", gterm.new, rterms.cur[rterms.cur != gterm.cur]), collapse = " + "))
            message("\nQTL x ", Trait, " Bivariate (corh) Random Effects Model.")
            cat("===============================================\n")
            qtlModel <- update(qtlModel, random. = ran.form, ...)
        } else {
            for (k in seq_len(n.fa)) {
                gterm.new <- if (k == 1L)
                    gsub("diag\\(", paste0("fa("), gsub("\\)$", paste0(",", 1L, ")"), gterm.cur))
                else
                    gsub(paste0(",", k - 1L, "\\)"), paste0(",", k, ")"), gterm.cur)
                ran.form  <- as.formula(paste(
                    c("~", gterm.new, rterms.cur[rterms.cur != gterm.cur]), collapse = " + "))
                message("\nQTL x ", Trait, " Factor Analytic(", k, ") Random Effects Model.")
                cat("===================================================\n")
                qtlModel  <- update(qtlModel, random. = ran.form, ...)
                gterm.cur <- gterm.new
            }
        }
        # Keep vmterms in sync with the upgraded term name
        if (vm) vmterms <- c(gterm.cur, merge.by)
    }

    list(qtlModel = qtlModel, intervalObj = intervalObj,
         cov.env = cov.env, vm = vm, vmterms = vmterms)
}

#' @keywords internal
.rebuildCovObj <- function(genoData, state, merge.by, intervalObj,
                            force, vm, vmterms, qtlModel, caller.env) {
    mout    <- (1:ncol(genoData))[!as.logical(state)]
    genoSub <- genoData[, -mout, drop = FALSE]
    cov.env <- NULL
    if ((ncol(genoSub) > nrow(genoSub)) & !force) {
        cov.env <- .constructCM(genoSub)
        covObj  <- cov.env$relm
        attr(intervalObj, "env") <- cov.env
    } else {
        covObj <- cbind.data.frame(rownames(genoSub), genoSub)
        names(covObj)[1] <- merge.by
        if (is.null(qtlModel$call$mbf$ints) & vm) {
            # vm -> mbf switch: reconstruct the random formula preserving any
            # Trait variance-structure suffix (e.g. ":diag(Trait)", ":corh(Trait)",
            # ":fa(Trait,k)") that was appended to the old vm term.
            attr(intervalObj, "env") <- NULL
            all.rterms <- unlist(strsplit(deparse(qtlModel$call$random[[2]]), " \\+ "))
            # Identify the genome-wide vm term specifically (starts with "vm(")
            vm.term      <- all.rterms[grep("^vm\\(", all.rterms)][1L]
            other.rterms <- all.rterms[all.rterms != vm.term]
            # Extract any Trait suffix from the old vm term and attach it to mbf
            trait.suffix <- sub(paste0("^vm\\(", merge.by, ",\\s*covObj\\)"), "", vm.term)
            mbf.term     <- paste0("mbf('ints')", trait.suffix)
            qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
            qtlModel$call$mbf$ints$cov <- "covObj"
            ran.form <- as.formula(paste("~",
                paste(c(mbf.term, merge.by, other.rterms), collapse = " + ")))
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
