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
    # Clamp negative eigenvalues (numerical noise from floating-point arithmetic).
    # Also apply a small unconditional ridge: as markers are progressively
    # excluded across iterations, the remaining genoData can become rank-deficient
    # (e.g. nearly collinear intervals near detected QTL), causing ASReml to
    # abort with "singularities in covObj". The ridge 1e-6 * I is negligible
    # relative to the mean diagonal of 1 but sufficient to keep K positive
    # definite throughout all iterations.
    sv$values[sv$values < 0] <- 0
    eps <- max(1e-6, 1e-10 * max(sv$values))
    relm <- sv$vectors %*% diag(sv$values + eps) %*% t(sv$vectors)
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

    # -------------------------------------------------------------------------
    # Extract the residual genetic term directly from the base model's random
    # formula.  The term containing merge.by is the genetic line term; all
    # other terms have already been removed into rterms by .fixLines().
    #
    # Examples of what this term may look like:
    #   Univariate  : "Variety"                  (bare line factor)
    #   Multivariate: "diag(Trial):Variety"
    #                 "corh(Trial):Variety"
    #                 "corgh(Trial):Variety"
    #                 "us(Trial):Variety"
    #                 "fa(Trial,2):Variety"
    #
    # (After fix.lines, Variety may have been replaced by "Gsave", but the
    # structure prefix, if any, is preserved by .fixLines.)
    # -------------------------------------------------------------------------
    all.rand  <- unlist(strsplit(deparse(baseModel$call$random[[2]]), " \\+ "))
    resid.term <- all.rand[grep(merge.by, all.rand)][1L]

    # Isolate the variance structure prefix by stripping ":merge.by" from the end.
    # Use a regex that escapes any special characters in merge.by.
    mb.esc        <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", merge.by)
    mb.pattern    <- paste0(":", mb.esc, "$")
    struct.prefix <- sub(mb.pattern, "", resid.term)
    has.struct    <- struct.prefix != resid.term   # FALSE for univariate bare term

    # Classify the residual structure to determine if/how to upgrade the
    # additive genomic term after the initial diag(Trait) model is fitted.
    #
    # corh / corgh / us  -> upgrade additive to corgh (heterogeneous corr/var)
    # fa(Trait, k)       -> upgrade additive to fa(Trait, n.fa)
    # diag / bare        -> no upgrade
    #
    # us(Trial) (unstructured) maps to corgh on the additive term: both span
    # the same space of positive-definite matrices but corgh is the natural
    # ASReml parameterisation for vm() / mbf() terms.
    upgrade.to.corgh <- has.struct &&
        grepl("^corh\\(|^corgh\\(|^us\\(", struct.prefix)
    upgrade.to.fa    <- has.struct && grepl("^fa\\(", struct.prefix) && n.fa > 0L

    # The additive genomic term always starts with diag(Trait): when multivariate
    # (regardless of the residual structure) and has no prefix for univariate.
    diag.prefix <- if (!is.null(Trait)) paste0("diag(", Trait, "):") else ""

    # -------------------------------------------------------------------------
    # Build and fit the initial genome-wide model.
    #   vm path  (markers > lines): vm(merge.by, covObj)
    #   mbf path (lines >= markers): mbf('ints')
    # Both get the diag(Trait): prefix for multivariate; no prefix for univariate.
    # -------------------------------------------------------------------------
    if ((ncol(genoData) > nrow(genoData)) & !force) {
        cov.env  <- .constructCM(genoData)
        covObj   <- cov.env$relm
        gterm    <- paste0(diag.prefix, "vm(", merge.by, ", covObj)")
        vmterms  <- c(gterm, resid.term)
        ran.form <- as.formula(paste(c("~", vmterms, rterms), collapse = " + "))
        attr(intervalObj, "env") <- cov.env
        vm <- TRUE
    } else {
        covObj <- cbind.data.frame(rownames(genoData), genoData)
        names(covObj)[1] <- merge.by
        qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
        qtlModel$call$mbf$ints$cov <- "covObj"
        gterm   <- paste0(diag.prefix, "mbf('ints')")
        vmterms <- c(gterm, resid.term)
        ran.form <- as.formula(paste(c(paste("~", gterm), resid.term, rterms),
                                     collapse = " + "))
    }
    assign("covObj", covObj, envir = caller.env)
    cat("\nRandom Effects Interval/Marker Model Iteration (1):\n")
    cat("============================================\n")
    qtlModel$call$data <- quote(phenoData)
    qtlModel <- update(qtlModel, random. = ran.form, ...)

    # -------------------------------------------------------------------------
    # Upgrade the additive genomic term's variance structure if required.
    # The residual genetic term (resid.term / vmterms[2]) is NEVER changed --
    # it is the user's input structure and stays as-is throughout.
    #
    #   corh / corgh / us in residual -> upgrade additive to corgh(Trait):vm/mbf
    #   fa(Trait,k)   in residual    -> upgrade additive to fa(Trait,n.fa):vm/mbf
    #                                   (step through fa(1)...fa(n.fa))
    #   diag / bare                  -> no upgrade; stay at diag(Trait): or bare
    #
    # Note: we carry the current additive term name in gterm (set above) and
    # build ran.form directly from it; we do NOT re-read from
    # qtlModel$call$random, which avoids any dependency on how ASReml stores
    # the formula internally after update().
    # -------------------------------------------------------------------------
    if (upgrade.to.corgh || upgrade.to.fa) {
        # other.rterms: all current random terms except the additive genomic term
        # and the residual genetic term (those are in vmterms).  At this point
        # only rterms (the non-genetic other terms) qualify.
        other.rterms <- rterms   # rterms are the non-genetic terms from .fixLines()
        gterm.cur    <- gterm    # carry the additive term name built above

        if (upgrade.to.corgh) {
            # diag(Trait):vm/mbf --> corgh(Trait):vm/mbf
            gterm.new <- sub(paste0("^diag\\(", Trait, "\\):"),
                             paste0("corgh(", Trait, "):"), gterm.cur)
            ran.form  <- as.formula(paste(c("~", gterm.new, resid.term, other.rterms),
                                          collapse = " + "))
            message("\nQTL x ", Trait, " corgh Random Effects Model.")
            cat("===============================================\n")
            qtlModel  <- update(qtlModel, random. = ran.form, ...)
            gterm.cur <- gterm.new

        } else {
            # diag(Trait):vm/mbf --> fa(Trait,1):vm/mbf --> ... --> fa(Trait,n.fa):vm/mbf
            for (k in seq_len(n.fa)) {
                old.struct <- if (k == 1L)
                    paste0("^diag\\(", Trait, "\\):")
                else
                    paste0("^fa\\(", Trait, ",\\s*", k - 1L, "\\):")
                new.struct <- paste0("fa(", Trait, ",", k, "):")
                gterm.new  <- sub(old.struct, new.struct, gterm.cur)
                ran.form   <- as.formula(paste(c("~", gterm.new, resid.term, other.rterms),
                                               collapse = " + "))
                message("\nQTL x ", Trait, " Factor Analytic(", k,
                        ") Random Effects Model.")
                cat("===================================================\n")
                qtlModel  <- update(qtlModel, random. = ran.form, ...)
                gterm.cur <- gterm.new
            }
        }

        # Keep vmterms[1] in sync with the upgraded additive term.
        # vmterms[2] (the residual genetic term) is never changed.
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
            # variance structure prefix (e.g. "diag(Trait):", "corgh(Trait):",
            # "fa(Trait,k):") that precedes the vm term in the additive genomic
            # term (vmterms[1]).  vmterms[2] is the residual genetic term
            # as extracted from the original base model; it is passed through
            # unchanged.
            attr(intervalObj, "env") <- NULL
            all.rterms <- unlist(strsplit(deparse(qtlModel$call$random[[2]]), " \\+ "))
            # Find the additive vm term in the current formula.  It may carry
            # any variance structure prefix so match on "vm.*covObj".
            vm.term      <- all.rterms[grep("vm.*covObj", all.rterms)][1L]
            # Collect all terms except both vmterms (additive + residual genetic).
            # These are the QTL effects and other random terms added during the loop.
            other.rterms <- all.rterms[!(all.rterms %in% vmterms)]
            # Extract the variance structure prefix from the vm term and
            # transplant it onto the new mbf term.
            vm.core      <- paste0("vm\\(", merge.by, ",\\s*covObj\\)$")
            struct.prefix <- sub(vm.core, "", vm.term)
            mbf.term     <- paste0(struct.prefix, "mbf('ints')")
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
