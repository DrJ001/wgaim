# =============================================================================
# GPAim.R
# Genomic Prediction via Integrated Modelling.
#
# Fits a genomic best linear unbiased prediction (G-BLUP) model and extracts
# GEBVs for all genotyped lines. Two computational paths are used:
#
#   vm path  (markers > lines): builds the genomic relationship matrix G = XX'
#            and fits vm(line, G). GEBVs extracted directly via predict().
#
#   mbf path (lines >= markers): models marker effects q directly as random
#            effects via ASReml's mbf() facility. GEBVs computed as M %*% q̂,
#            where M is the marker matrix and q̂ are the marker effect BLUPs.
#            The mbf path avoids the singularity of G when lines >= markers.
#
# Arguments intentionally absent (vs QTLAim/GWASAim):
#   method          - always 'random' (GEBVs are random effects by definition)
#   selection       - no selection; no loop
#   exclusion.window- no iterative selection to exclude around
#   breakout        - no loop to break out of
#   TypeI           - no significance testing
# =============================================================================

GPAim <- function(baseModel, ...)
    UseMethod("GPAim")

GPAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

GPAim.asreml <- function(baseModel, genoObj, merge.by = NULL,
                          fix.lines = TRUE, gen.type = "marker",
                          force = FALSE, trace = TRUE, ...) {

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
    # Phase 3: Fit GP model — vm or mbf path
    #
    # vm  path (markers > lines): G = XX' is full rank; fit vm(line, G).
    # mbf path (lines >= markers): G is singular; fit marker effects directly
    #                              via mbf(), then GEBVs = M %*% q̂.
    # -------------------------------------------------------------------------
    n.markers    <- ncol(genoData)
    n.lines.geno <- nrow(genoData)
    use.vm       <- (n.markers > n.lines.geno) & !force

    gpModel           <- baseModel
    gpModel$call$data <- quote(phenoData)

    if (use.vm) {
        cat(sprintf("\nvm path: building relationship matrix (%d markers > %d lines)...\n",
                    n.markers, n.lines.geno))
        cov.env  <- .constructCM(genoData)
        covObj   <- cov.env$relm
        assign("covObj", covObj, envir = caller.env)
        vmterm   <- paste0("vm(", merge.by, ", covObj)")
        ran.form <- as.formula(paste(c("~", vmterm, rterms), collapse = " + "))
    } else {
        cat(sprintf("\nmbf path: fitting marker effects directly (%d lines >= %d markers)...\n",
                    n.lines.geno, n.markers))
        cov.env  <- NULL
        covObj   <- cbind.data.frame(rownames(genoData), genoData)
        names(covObj)[1] <- merge.by
        gpModel$call$mbf$markers$key <- rep(merge.by, 2)
        gpModel$call$mbf$markers$cov <- "covObj"
        assign("covObj", covObj, envir = caller.env)
        ran.form <- as.formula(paste(c("~ mbf('markers')", rterms), collapse = " + "))
    }

    cat("Fitting Genomic Prediction model...\n")
    gpModel <- update(gpModel, random. = ran.form, ...)

    # -------------------------------------------------------------------------
    # Extract GEBVs and compute variance components
    # -------------------------------------------------------------------------
    cat("Extracting GEBVs...\n")
    sigma2 <- gpModel$sigma2
    if (gpModel$vparameters.con[length(gpModel$vparameters.con)] == 4)
        sigma2 <- 1

    if (use.vm) {
        # vm path: predict genetic values directly from the vm term
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

        var.genetic <- sigma2 *
            gpModel$vparameters[grep("vm.*covObj", names(gpModel$vparameters))]
        var.resid   <- sigma2

    } else {
        # mbf path: GEBVs = M %*% q̂
        # Extract marker effect BLUPs (q̂) from the mbf random coefficients
        mbf.rows <- grep("mbf", rownames(gpModel$coefficients$random))
        q.hat    <- gpModel$coefficients$random[mbf.rows, 1]
        gebvs    <- as.numeric(genoData %*% q.hat)

        # Approximate SE: sqrt( M^2 %*% PEV(q̂) )
        # Treats marker effects as independent — ignores their covariance.
        pev     <- sigma2 * gpModel$vcoeff$random[mbf.rows]
        se.gebv <- sqrt(as.numeric(genoData^2 %*% pev))

        gebv <- data.frame(
            rownames(genoData),
            gebvs,
            se.gebv,
            stringsAsFactors = FALSE
        )
        names(gebv) <- c(genetic.term, "GEBV", "SE")

        # Genetic variance from marker effects: Var(Mq) = sigma2 * vpar * mean(sum_j m_ij^2)
        vpar.mbf    <- gpModel$vparameters[grep("mbf.*markers", names(gpModel$vparameters))]
        var.genetic <- sigma2 * vpar.mbf * mean(rowSums(genoData^2), na.rm = TRUE)
        var.resid   <- sigma2
    }

    h2 <- as.numeric(var.genetic / (var.genetic + var.resid))

    # -------------------------------------------------------------------------
    # Package results and clean up
    # -------------------------------------------------------------------------
    gp.list <- list(
        gebv         = gebv,
        gen.type     = gen.type,
        path         = ifelse(use.vm, "vm", "mbf"),
        var.genetic  = as.numeric(var.genetic),
        var.resid    = as.numeric(var.resid),
        heritability = h2,
        n.markers    = n.markers,
        rel.scale    = if (use.vm) cov.env$scale else 1,
        genetic.term = genetic.term
    )

    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    gpModel <- .envFix(gpModel, asremlEnv)
    gpModel$GP <- gp.list
    class(gpModel) <- c("GPAim", "asreml")
    gpModel
}
