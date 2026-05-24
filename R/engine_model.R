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
                               Trait = NULL, str = NULL, ...) {
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
    all.rand   <- .rhs_terms(baseModel$call$random)
    resid.term <- all.rand[grep(merge.by, all.rand)][1L]

    # Isolate the variance structure prefix by stripping ":merge.by" from the end.
    mb.esc        <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", merge.by)
    mb.pattern    <- paste0(":", mb.esc, "$")
    struct.prefix <- sub(mb.pattern, "", resid.term)
    has.struct    <- struct.prefix != resid.term   # FALSE for univariate bare term

    # -------------------------------------------------------------------------
    # Determine upgrade flags and n.fa from str (explicit) or residual (mirror).
    #
    # str = NULL  -> mirror the residual structure exactly:
    #   corh / corgh / us in residual  -> upgrade additive to corgh
    #   fa(Trait, k) in residual       -> upgrade additive to fa(Trait, k)
    #   diag / bare                    -> no upgrade
    #
    # str non-NULL -> override, regardless of residual structure:
    #   "corh" / "corgh" / "us"  -> upgrade additive to corgh
    #   "fa2", "fa3", ...         -> upgrade additive to fa(Trait, k)
    #   "diag"                    -> no upgrade (independent per-trait)
    #   "fa" (no number)          -> error: number of factors must be specified
    #
    # us(Trial) maps to corgh on the additive term: both span the same space of
    # positive-definite matrices but corgh is the natural ASReml parameterisation
    # for vm() / mbf() terms.
    # -------------------------------------------------------------------------
    # vm.struct: the exact variance structure prefix to apply to the additive
    # genomic term.  NULL = no upgrade (diag or bare univariate).
    # upgrade.to.fa: TRUE when vm.struct == "fa".
    n.fa         <- 0L
    vm.struct    <- NULL   # NULL -> no upgrade
    upgrade.to.fa <- FALSE

    if (is.null(str)) {
        # Mirror from residual: extract the structure name from struct.prefix.
        # For us() in the residual we mirror to "corgh" -- us() is not supported
        # as a variance structure on vm()/mbf() terms in ASReml.
        if (has.struct) {
            if (grepl("^corh\\(", struct.prefix))
                vm.struct <- "corh"
            else if (grepl("^corgh\\(", struct.prefix))
                vm.struct <- "corgh"
            else if (grepl("^us\\(", struct.prefix))
                vm.struct <- "corgh"   # us not supported on vm; fall back to corgh
            else if (grepl("^fa\\(", struct.prefix)) {
                vm.struct     <- "fa"
                n.fa          <- as.integer(sub(".*,\\s*(\\d+)\\).*", "\\1", struct.prefix))
                upgrade.to.fa <- TRUE
            }
            # diag() -> vm.struct stays NULL (no upgrade)
        }
    } else {
        # Explicit override: use exactly what the user asked for.
        str.l <- tolower(trimws(str))
        if (str.l == "corh") {
            vm.struct <- "corh"
        } else if (str.l == "corgh") {
            vm.struct <- "corgh"
        } else if (str.l == "us") {
            vm.struct <- "us"
        } else if (str.l == "diag") {
            vm.struct <- NULL   # no upgrade
        } else if (grepl("^fa[0-9]+$", str.l)) {
            n.fa          <- as.integer(sub("^fa", "", str.l))
            vm.struct     <- "fa"
            upgrade.to.fa <- TRUE
        } else if (str.l == "fa") {
            stop("str = \"fa\" requires a number of factors, e.g. str = \"fa2\".")
        } else {
            stop("str must be NULL, \"diag\", \"corh\", \"corgh\", \"us\", ",
                 "or \"fa1\", \"fa2\", etc. Got: \"", str, "\".")
        }
    }

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

    # The additive genomic term is returned with diag(Trait): structure only.
    # Upgrading to corh/corgh/fa is deferred to .upgradeVmStructure(), called
    # by the analysis function after an initial LRT confirms there is signal.
    list(qtlModel = qtlModel, intervalObj = intervalObj,
         cov.env = cov.env, vm = vm, vmterms = vmterms, n.fa = n.fa,
         vm.struct = vm.struct, upgrade.to.fa = upgrade.to.fa,
         gterm = gterm, resid.term = resid.term)
}

# =============================================================================
# .upgradeVmStructure
#
# Upgrades the additive genomic term's variance structure from the initial
# diag(Trait): form to the requested corh/corgh/fa(k) structure.
#
# Called by qtlAim.asreml() and gwasAim.asreml() AFTER the initial LRT has
# confirmed significant genome-wide additive variance, so that complex
# unstructured models are only fitted when there is evidence of signal.
#
# For ntrait == 1 or vm.struct == NULL (diag or bare univariate): no-op.
#
# Arguments:
#   qtlModel     : the diag(Trait):vm model returned by .buildGenomeModel()
#   vm.struct    : character upgrade target ("corh","corgh","us","fa") or NULL
#   n.fa         : integer, number of FA factors (0 if not FA)
#   upgrade.to.fa: logical, TRUE when vm.struct == "fa"
#   Trait        : character, name of the Trait factor column
#   rterms       : character vector of other (non-genetic) random terms
#   resid.term   : character, the residual genetic term (vmterms[2L])
#   gterm        : character, the current additive term name (diag prefix)
#   vm           : logical, TRUE = vm path; FALSE = mbf path
#   vmterms      : character(2) — updated in place and returned
#   phenoData    : phenotypic data frame (must be in scope for ASReml update())
#   ...          : further args passed to update()
#
# Returns: list(qtlModel, vmterms)
# =============================================================================
#' @keywords internal
.upgradeVmStructure <- function(qtlModel, vm.struct, n.fa, upgrade.to.fa,
                                 Trait, rterms, resid.term, gterm,
                                 vm, vmterms, phenoData, ...) {
    if (is.null(vm.struct))
        return(list(qtlModel = qtlModel, vmterms = vmterms))

    other.rterms <- rterms
    gterm.cur    <- gterm   # starts as diag(Trait):vm or diag(Trait):mbf

    if (!upgrade.to.fa) {
        # diag(Trait):vm/mbf --> vm.struct(Trait):vm/mbf
        gterm.new <- sub(paste0("^diag\\(", Trait, "\\):"),
                         paste0(vm.struct, "(", Trait, "):"), gterm.cur)
        ran.form  <- as.formula(paste(c("~", gterm.new, resid.term, other.rterms),
                                      collapse = " + "))
        message("\nQTL x ", Trait, " ", vm.struct, " Random Effects Model.")
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
    if (vm) vmterms[1L] <- gterm.cur
    list(qtlModel = qtlModel, vmterms = vmterms)
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
            all.rterms <- .rhs_terms(qtlModel$call$random)
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

# =============================================================================
# .rhs_terms() -- robustly split a random formula RHS into individual terms.
#
# Problem: deparse() wraps long expressions into a character vector of length
# > 1 when the expression exceeds width.cutoff (default 60).  Continuation
# lines are indented with leading spaces.  Naively splitting on " + " then
# fails because:
#   - The wrap boundary produces " +\n    " not " + ", so the pattern misses
#     the split entirely and a mangled term containing the newline survives.
#   - Or the wrap produces a trailing empty element at end of line 1.
#   - Leading spaces on continuation lines leave " Day" instead of "Day".
#
# This affects any model whose random formula exceeds ~60 characters, which
# is common in field trials with spatial terms (at(Day):Plate:PCol etc.).
#
# Fix: use width.cutoff = 500L (no practical formula will exceed this), join
# the elements with a space, then split on whitespace-tolerant \\s*+\\s* and
# drop any empty strings that survive.
# =============================================================================
#' @keywords internal
.rhs_terms <- function(formula.obj) {
    raw   <- paste(deparse(formula.obj[[2L]], width.cutoff = 500L), collapse = " ")
    terms <- trimws(unlist(strsplit(raw, "\\s*\\+\\s*")))
    terms[nchar(terms) > 0L]
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
