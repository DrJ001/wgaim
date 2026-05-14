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

    # Build the Trait variance structure for the genome-wide and residual polygenic terms.
    #
    # The correct ASReml formula for a multivariate (MET) model is:
    #   ~ diag(Trait):vm(merge.by, covObj) + diag(Trait):merge.by + <other rterms>
    # where:
    #   diag(Trait):vm(merge.by, covObj) -- composite additive G x E term
    #   diag(Trait):merge.by             -- residual polygenic G x E term
    #
    # The Trait variance structure is a PREFIX to both the genomic and
    # residual terms. When Trait = NULL (univariate) no prefix is applied
    # and the formula reduces to the existing ~ vm(merge.by, covObj) + merge.by.
    trait.prefix  <- if (!is.null(Trait)) paste0("diag(", Trait, "):") else ""
    resid.term    <- paste0(trait.prefix, merge.by)   # residual polygenic term

    if ((ncol(genoData) > nrow(genoData)) & !force) {
        cov.env <- .constructCM(genoData)
        covObj  <- cov.env$relm
        gterm   <- paste0(trait.prefix, "vm(", merge.by, ", covObj)")
        vmterms <- c(gterm, resid.term)
        ran.form <- as.formula(paste(c("~", vmterms, rterms), collapse = " + "))
        attr(intervalObj, "env") <- cov.env
        vm <- TRUE
    } else {
        covObj <- cbind.data.frame(rownames(genoData), genoData)
        names(covObj)[1] <- merge.by
        qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
        qtlModel$call$mbf$ints$cov <- "covObj"
        gterm    <- paste0(trait.prefix, "mbf('ints')")
        ran.form <- as.formula(paste(c(paste("~", gterm), resid.term, rterms), collapse = " + "))
    }
    assign("covObj", covObj, envir = caller.env)
    cat("\nRandom Effects Interval/Marker Model Iteration (1):\n")
    cat("============================================\n")
    qtlModel$call$data <- quote(phenoData)
    qtlModel <- update(qtlModel, random. = ran.form, ...)

    # Multivariate only: upgrade the genome-wide term's variance structure
    # from diag to corh (ntrait=2) or fa(k) (ntrait>=3).
    # Only the genomic term is upgraded; the residual diag(Trait):merge.by stays.
    if (!is.null(Trait) && n.fa > 0L) {
        ntrait <- length(levels(phenoData[[Trait]]))
        rterms.cur <- attr(terms.formula(qtlModel$call$random), "term.labels")
        gterm.cur  <- rterms.cur[grep("vm.*covObj|mbf.*ints", rterms.cur)]
        if (ntrait == 2L) {
            # diag(Trait):vm() --> corh(Trait):vm()
            gterm.new  <- sub(paste0("^diag\\(", Trait, "\\):"),
                              paste0("corh(", Trait, "):"), gterm.cur)
            ran.form   <- as.formula(paste(
                c("~", gterm.new, rterms.cur[rterms.cur != gterm.cur]), collapse = " + "))
            message("\nQTL x ", Trait, " Bivariate (corh) Random Effects Model.")
            cat("===============================================\n")
            qtlModel <- update(qtlModel, random. = ran.form, ...)
            gterm.cur <- gterm.new
        } else {
            for (k in seq_len(n.fa)) {
                # diag(Trait): --> fa(Trait,k):  targeting the exact prefix
                old.struct <- if (k == 1L)
                    paste0("^diag\\(", Trait, "\\):")
                else
                    paste0("^fa\\(", Trait, ",\\s*", k - 1L, "\\):")
                new.struct <- paste0("fa(", Trait, ",", k, "):")
                gterm.new  <- sub(old.struct, new.struct, gterm.cur)
                ran.form   <- as.formula(paste(
                    c("~", gterm.new, rterms.cur[rterms.cur != gterm.cur]), collapse = " + "))
                message("\nQTL x ", Trait, " Factor Analytic(", k, ") Random Effects Model.")
                cat("===================================================\n")
                qtlModel  <- update(qtlModel, random. = ran.form, ...)
                gterm.cur <- gterm.new
            }
        }
        # Keep vmterms[1] in sync with the upgraded genomic term name;
        # vmterms[2] (residual polygenic) is unchanged.
        if (vm) vmterms[1L] <- gterm.cur
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
            # vm -> mbf switch: reconstruct the random formula, preserving the
            # Trait variance-structure PREFIX (e.g. "diag(Trait):", "corh(Trait):",
            # "fa(Trait,k):") that precedes the old vm term.
            attr(intervalObj, "env") <- NULL
            all.rterms <- unlist(strsplit(deparse(qtlModel$call$random[[2]]), " \\+ "))
            # Find the genome-wide vm term -- may now start with a variance structure
            # prefix such as "diag(Site):vm(...)", so match on "vm.*covObj" not "^vm("
            vm.term      <- all.rterms[grep("vm.*covObj", all.rterms)][1L]
            # Everything else is kept as-is (includes residual diag(Trait):merge.by)
            other.rterms <- all.rterms[!(all.rterms %in% vmterms)]
            # Extract the Trait prefix (e.g. "diag(Site):") from the vm term
            # and transplant it onto the mbf term.
            vm.core      <- paste0("vm\\(", merge.by, ",\\s*covObj\\)$")
            trait.prefix <- sub(vm.core, "", vm.term)
            mbf.term     <- paste0(trait.prefix, "mbf('ints')")
            qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
            qtlModel$call$mbf$ints$cov <- "covObj"
            ran.form <- as.formula(paste("~",
                paste(c(mbf.term, vmterms[2L], other.rterms), collapse = " + ")))
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
