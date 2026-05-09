# =============================================================================
# GPAim.R
# Genomic Prediction via Integrated Modelling.
#
# Fits a genomic best linear unbiased prediction (G-BLUP) model using the
# shared wgAim engine. Unlike QTLAim and GWASAim there is no iterative
# selection loop — GPAim fits a single genome-wide model using the genomic
# relationship matrix (G matrix) as the covariance structure, then extracts
# genomic estimated breeding values (GEBVs) for all genotyped lines.
#
# The vm path (relationship matrix) is always used regardless of the
# markers:lines ratio — the G matrix IS the genomic prediction model.
#
# Arguments intentionally absent (vs QTLAim/GWASAim):
#   method          - always 'random' (GEBVs are random effects by definition)
#   selection       - no selection; no loop
#   force           - vm path is always taken; not a user decision
#   exclusion.window- no iterative selection to exclude around
#   breakout        - no loop to break out of
#   TypeI           - no significance testing
#   bonferroni      - removed from GWASAim; not applicable here either
# =============================================================================

GPAim <- function(baseModel, ...)
    UseMethod("GPAim")

GPAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

GPAim.asreml <- function(baseModel, genoObj, merge.by = NULL,
                          fix.lines = TRUE, gen.type = "marker",
                          trace = TRUE, ...) {

    caller.env <- parent.frame()

    # Trace/sink must live here so on.exit fires at the right level
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }

    # -------------------------------------------------------------------------
    # Phase 1: Validation
    # Inline only the shared parts — GPAim has no method/selection/breakout
    # -------------------------------------------------------------------------
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
        if (!baseModel$converge)
            stop("Base model not converged: Check base model before proceeding.")
    }
    asremlEnv <- lapply(baseModel$formulae, function(el) attr(el, ".Environment"))
    phenoData <- eval(baseModel$call$data)

    if (is.null(merge.by))
        stop("merge.by: name of the column linking phenotypic and genotypic data is required.")
    if (!(gen.type %in% c("marker", "interval")))
        stop("gen.type must be \"marker\" or \"interval\".")

    # genoObj: accept either 'interval' (from cross2int) or 'panel' (from makePanel)
    if (missing(genoObj))
        stop("genoObj is a required argument. Supply an 'interval' or 'panel' object.")
    if (!inherits(genoObj, "interval"))
        stop("genoObj must be of class \"interval\" or \"panel\".")
    if (gen.type == "interval" && inherits(genoObj, "panel"))
        stop("gen.type = \"interval\" requires an 'interval' object from cross2int(), ",
             "not a 'panel' object. Use gen.type = \"marker\" with makePanel() output.")

    glines <- genoObj$pheno[, merge.by]
    if (is.null(glines))
        stop("Genotypic data does not contain column \"", merge.by, "\".")
    plines <- phenoData[, merge.by]
    if (is.null(plines))
        stop("Phenotypic data does not contain column \"", merge.by, "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in genotypic \"", merge.by, "\" column do not match any names ",
             "in phenotypic \"", merge.by, "\" column.")

    # -------------------------------------------------------------------------
    # Phase 2: Build genotype data matrix
    # -------------------------------------------------------------------------
    gd       <- .buildGenoData(genoObj, gen.type, glines, plines)
    genoData <- gd$genoData

    fl           <- .fixLines(baseModel, phenoData, genoData, merge.by, plines, fix.lines, ...)
    baseModel    <- fl$baseModel
    phenoData    <- fl$phenoData
    merge.by     <- fl$merge.by
    rterms       <- fl$rterms
    genetic.term <- fl$genetic.term

    # -------------------------------------------------------------------------
    # Phase 3: Build genomic relationship matrix and fit GP model
    # Always uses the vm path — the G matrix is the GP covariance structure.
    # -------------------------------------------------------------------------
    n.markers <- ncol(genoData)
    cat(sprintf("\nBuilding genomic relationship matrix  (%d markers)...\n", n.markers))
    cov.env <- .constructCM(genoData)
    covObj  <- cov.env$relm
    assign("covObj", covObj, envir = caller.env)

    vmterm   <- paste0("vm(", merge.by, ", covObj)")
    ran.form <- as.formula(paste(c("~", vmterm, rterms), collapse = " + "))

    cat("Fitting Genomic Prediction model...\n")
    gpModel           <- baseModel
    gpModel$call$data <- quote(phenoData)
    gpModel <- update(gpModel, random. = ran.form, ...)

    # -------------------------------------------------------------------------
    # Extract GEBVs via predict() on the vm term
    # -------------------------------------------------------------------------
    cat("Extracting GEBVs...\n")
    pv   <- predict(gpModel, classify = merge.by, only = vmterm,
                    vcov = FALSE, data = phenoData)
    pvdf <- pv$pvals
    gebv <- data.frame(
        pvdf[[merge.by]],
        pvdf[["predicted.value"]],
        pvdf[["std.error"]],
        stringsAsFactors = FALSE
    )
    names(gebv) <- c(genetic.term, "GEBV", "SE")

    # -------------------------------------------------------------------------
    # Variance components and narrow-sense heritability
    # -------------------------------------------------------------------------
    sigma2 <- gpModel$sigma2
    if (gpModel$vparameters.con[length(gpModel$vparameters.con)] == 4)
        sigma2 <- 1
    var.genetic <- sigma2 *
        gpModel$vparameters[grep("vm.*covObj", names(gpModel$vparameters))]
    var.resid   <- sigma2   # residual variance in ASReml parameterisation
    h2          <- as.numeric(var.genetic / (var.genetic + var.resid))

    # -------------------------------------------------------------------------
    # Package results and clean up
    # -------------------------------------------------------------------------
    gp.list <- list(
        gebv         = gebv,
        gen.type     = gen.type,
        var.genetic  = as.numeric(var.genetic),
        var.resid    = as.numeric(var.resid),
        heritability = h2,
        n.markers    = n.markers,
        rel.scale    = cov.env$scale,
        genetic.term = genetic.term
    )

    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    gpModel <- .envFix(gpModel, asremlEnv)
    gpModel$GP <- gp.list
    class(gpModel) <- c("GPAim", "asreml")
    gpModel
}
