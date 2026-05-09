# =============================================================================
# GWASAim.R
# Whole Genome GWAS Analyses via Integrated Modelling.
#
# Forward-selection marker association analysis using the shared wgAim engine.
# Differences from QTLAim:
#   - Takes a 'panel' object (from makePanel()) instead of an 'interval' object
#   - gen.type is always 'marker' (no interval midpoints in GWAS)
#   - Bonferroni correction applied to TypeI threshold by default
#   - Returns class c('GWASAim', 'asreml')
# =============================================================================

GWASAim <- function(baseModel, ...)
    UseMethod("GWASAim")

GWASAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

GWASAim.asreml <- function(baseModel, panelObj, merge.by = NULL,
                            fix.lines = TRUE, method = "fixed",
                            selection = "interval", force = FALSE,
                            exclusion.window = 20, breakout = -1,
                            TypeI = 0.05, bonferroni = TRUE,
                            trace = TRUE, verboseLev = 0, ...) {

    # Capture calling environment early — needed for assign() in engine
    caller.env <- parent.frame()

    # Trace/sink setup must live here so on.exit fires at the right level
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }

    # -------------------------------------------------------------------------
    # Phase 1: Shared model validation
    # -------------------------------------------------------------------------
    vd        <- .validateModel(baseModel, merge.by, method, selection, breakout)
    baseModel <- vd$baseModel
    asremlEnv <- vd$asremlEnv
    phenoData <- vd$phenoData

    # GWAS-specific: panelObj validation and line matching
    if (missing(panelObj))
        stop("panelObj is a required argument. Use makePanel() to create one.")
    if (!inherits(panelObj, "panel"))
        stop("panelObj must be of class \"panel\". Use makePanel() to create one.")
    glines <- panelObj$pheno[, merge.by]
    if (is.null(glines))
        stop("Panel data does not contain column \"", merge.by, "\".")
    plines <- phenoData[, merge.by]
    if (is.null(plines))
        stop("Phenotypic data does not contain column \"", merge.by, "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in panel \"", merge.by, "\" column do not match any names ",
             "in phenotypic \"", merge.by, "\" column.")

    # -------------------------------------------------------------------------
    # Phase 2a: Build genotype data matrix (marker mode forced for GWAS)
    # -------------------------------------------------------------------------
    gd       <- .buildGenoData(panelObj, "marker", glines, plines)
    genoData <- gd$genoData
    mnams    <- gd$mnams
    state    <- gd$state

    # Apply Bonferroni correction to significance threshold
    n.markers <- ncol(genoData)
    TypeI.eff <- if (bonferroni) TypeI / n.markers else TypeI
    if (bonferroni) {
        cat("\nBonferroni-adjusted significance threshold:",
            formatC(TypeI.eff, format = "e", digits = 3),
            sprintf("(TypeI=%.3f / %d markers)\n", TypeI, n.markers))
    }

    # Phase 2b: Handle lines present in phenotypic but absent from panel
    fl           <- .fixLines(baseModel, phenoData, genoData, merge.by, plines, fix.lines, ...)
    baseModel    <- fl$baseModel
    phenoData    <- fl$phenoData
    merge.by     <- fl$merge.by
    rterms       <- fl$rterms
    genetic.term <- fl$genetic.term

    # -------------------------------------------------------------------------
    # Phase 3: Build and fit initial genome-wide marker model (vm or mbf path)
    # -------------------------------------------------------------------------
    gm        <- .buildGenomeModel(baseModel, genoData, phenoData, merge.by,
                                   panelObj, force, rterms, caller.env, ...)
    qtlModel  <- gm$qtlModel
    panelObj  <- gm$intervalObj   # may have env attribute set (vm path)
    cov.env   <- gm$cov.env
    vm        <- gm$vm
    vmterms   <- gm$vmterms

    # -------------------------------------------------------------------------
    # Phase 4: Iterative forward-selection loop
    # -------------------------------------------------------------------------
    ldiag <- coef.list <- vcoef.list <- list()
    qtl   <- c()
    iter  <- 1

    repeat {
        # Compute outlier statistics and select best marker
        selq               <- .qtlSelect(qtlModel, phenoData, panelObj, "marker",
                                         selection, exclusion.window, state, verboseLev)
        state              <- selq$state
        ldiag$oint[[iter]] <- selq$oint
        ldiag$ochr[[iter]] <- selq$ochr
        ldiag$blups[[iter]] <- selq$blups

        # Likelihood ratio test with Bonferroni-adjusted threshold
        lrt               <- .lrtTest(qtlModel, baseModel, TypeI.eff)
        ldiag$lik[[iter]] <- c(lrt$baseLogL, qtlModel$loglik, lrt$stat, lrt$pvalue)
        if (!lrt$pass | breakout == iter) break

        # Record significant marker and report
        qtl[iter] <- selq$qtl
        cqtl <- strsplit(qtl[iter], "\\.")
        message("Found significant marker on chromosome ", sapply(cqtl, "[", 2),
                " marker ", sapply(cqtl, "[", 3))

        # Merge marker genotype column into phenoData
        me        <- .mergeEffect(phenoData, genoData, qtl[iter], merge.by)
        phenoData <- me$phenoData
        qtl.x     <- me$qtl.x

        # Rebuild covariance object with selected marker excluded
        rc        <- .rebuildCovObj(genoData, state, merge.by, panelObj,
                                    force, vm, vmterms, qtlModel, caller.env)
        cov.env   <- rc$cov.env
        panelObj  <- rc$intervalObj
        qtlModel  <- rc$qtlModel

        qtlModel$call$data <- baseModel$call$data <- quote(phenoData)

        # Add selected marker effect (fixed or random) to both models
        ae                 <- .addEffect(baseModel, qtlModel, phenoData, merge.by,
                                         qtl.x, method, iter, ...)
        baseModel          <- ae$baseModel
        qtlModel           <- ae$qtlModel
        coef.list[[iter]]  <- ae$coefs
        vcoef.list[[iter]] <- ae$vcoefs

        iter <- iter + 1
    }

    # -------------------------------------------------------------------------
    # Phase 5: Package results and clean up
    # -------------------------------------------------------------------------
    qtl.list            <- .packResults(qtl, coef.list, vcoef.list, ldiag, state,
                                         iter, breakout, cov.env, genetic.term,
                                         method, "marker", selection)
    qtl.list$bonferroni <- bonferroni
    qtl.list$TypeI.eff  <- TypeI.eff
    qtl.list$n.markers  <- n.markers

    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    qtlModel <- .envFix(qtlModel, asremlEnv)
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("GWASAim", "asreml")
    qtlModel
}
