# =============================================================================
# QTLAim.R
# Whole Genome QTL Analyses via Integrated Modelling.
#
# User-facing orchestrator for QTL detection using a forward-selection
# iterative approach. All heavy lifting is delegated to the shared engine
# pieces in engine_*.R. This function is responsible only for:
#   - QTL-specific argument validation
#   - trace/sink setup (on.exit must be in the outermost frame)
#   - Calling engine pieces in the correct order
#   - Assigning results back to the calling environment
# =============================================================================

QTLAim <- function(baseModel, ...)
    UseMethod("QTLAim")

QTLAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

QTLAim.asreml <- function(baseModel, intervalObj, merge.by = NULL, fix.lines = TRUE,
                           gen.type = "interval", method = "fixed",
                           selection = "interval", force = FALSE,
                           exclusion.window = 20, breakout = -1,
                           TypeI = 0.05, trace = TRUE, verboseLev = 0, ...) {

    # Capture calling environment early — needed for assign() calls in engine
    caller.env <- parent.frame()

    # Trace/sink setup must live here so on.exit fires when QTLAim.asreml() returns
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }

    # -------------------------------------------------------------------------
    # Phase 1: Shared model validation
    # -------------------------------------------------------------------------
    vd <- .validateModel(baseModel, merge.by, method, selection, breakout)
    baseModel <- vd$baseModel
    asremlEnv <- vd$asremlEnv
    phenoData <- vd$phenoData

    # QTL-specific: intervalObj validation and line matching
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj must be of class \"interval\".")
    glines <- intervalObj$pheno[, merge.by]
    if (is.null(glines))
        stop("Genotypic data does not contain column \"", merge.by, "\".")
    plines <- phenoData[, merge.by]
    if (is.null(plines))
        stop("Phenotypic data does not contain column \"", merge.by, "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in genotypic \"", merge.by, "\" column do not match any names ",
             "in phenotypic \"", merge.by, "\" column.")

    # -------------------------------------------------------------------------
    # Phase 2a: Build genotype data matrix
    # -------------------------------------------------------------------------
    gd       <- .buildGenoData(intervalObj, gen.type, glines, plines)
    genoData <- gd$genoData
    mnams    <- gd$mnams
    state    <- gd$state

    # Phase 2b: Handle lines present in phenotypic but absent from genotypic data
    fl           <- .fixLines(baseModel, phenoData, genoData, merge.by, plines, fix.lines, ...)
    baseModel    <- fl$baseModel
    phenoData    <- fl$phenoData
    merge.by     <- fl$merge.by
    rterms       <- fl$rterms
    genetic.term <- fl$genetic.term

    # -------------------------------------------------------------------------
    # Phase 3: Build and fit initial genome-wide model (vm or mbf path)
    # -------------------------------------------------------------------------
    gm          <- .buildGenomeModel(baseModel, genoData, phenoData, merge.by,
                                     intervalObj, force, rterms, caller.env, ...)
    qtlModel    <- gm$qtlModel
    intervalObj <- gm$intervalObj
    cov.env     <- gm$cov.env
    vm          <- gm$vm
    vmterms     <- gm$vmterms

    # -------------------------------------------------------------------------
    # Phase 4: Iterative forward-selection loop
    # -------------------------------------------------------------------------
    ldiag <- coef.list <- vcoef.list <- list()
    qtl   <- c()
    iter  <- 1

    repeat {
        # Compute outlier statistics and select best interval/marker
        selq               <- .qtlSelect(qtlModel, phenoData, intervalObj, gen.type,
                                         selection, exclusion.window, state, verboseLev)
        state              <- selq$state
        ldiag$oint[[iter]] <- selq$oint
        ldiag$ochr[[iter]] <- selq$ochr
        ldiag$blups[[iter]] <- selq$blups

        # Likelihood ratio test against base model
        lrt               <- .lrtTest(qtlModel, baseModel, TypeI)
        ldiag$lik[[iter]] <- c(lrt$baseLogL, qtlModel$loglik, lrt$stat, lrt$pvalue)
        if (!lrt$pass | breakout == iter) break

        # Record selected QTL and report
        qtl[iter] <- selq$qtl
        cqtl <- strsplit(qtl[iter], "\\.")
        message("Found QTL on chromosome ", sapply(cqtl, "[", 2),
                " ", gen.type, " ", sapply(cqtl, "[", 3))

        # Merge selected genotype column into phenoData
        me        <- .mergeEffect(phenoData, genoData, qtl[iter], merge.by)
        phenoData <- me$phenoData
        qtl.x     <- me$qtl.x

        # Rebuild covariance object with selected interval excluded
        rc          <- .rebuildCovObj(genoData, state, merge.by, intervalObj,
                                      force, vm, vmterms, qtlModel, caller.env)
        cov.env     <- rc$cov.env
        intervalObj <- rc$intervalObj
        qtlModel    <- rc$qtlModel

        qtlModel$call$data <- baseModel$call$data <- quote(phenoData)

        # Add selected effect (fixed or random) to both models
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
    qtl.list <- .packResults(qtl, coef.list, vcoef.list, ldiag, state, iter,
                              breakout, cov.env, genetic.term, method, gen.type, selection)

    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    qtlModel <- .envFix(qtlModel, asremlEnv)
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("QTLAim", "asreml")
    qtlModel
}
