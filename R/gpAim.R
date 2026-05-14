# =============================================================================
# gpAim.R
# Genomic Prediction via Integrated Modelling.
#
# Fits a genomic best linear unbiased prediction (G-BLUP) model and extracts
# GEBVs for all genotyped lines. Two computational paths are used:
#
#   vm path  (markers > lines): builds the genomic relationship matrix G = XX'
#            and fits vm(line, G). GEBVs extracted directly via predict().
#
#   mbf path (lines >= markers): models marker effects q directly as random
#            effects via ASReml's mbf() facility. GEBVs computed as M %*% q.hat,
#            where M is the marker matrix and q.hat are the marker effect BLUPs.
#            The mbf path avoids the singularity of G when lines >= markers.
#
# genObj must be of class "wgCross" (from primeCross()) or "wgPanel"
# (from primePanel()).
#
# Arguments intentionally absent (vs qtlAim/gwasAim):
#   method          - always 'random' (GEBVs are random effects by definition)
#   selection       - no selection; no loop
#   exclusion.window- no iterative selection to exclude around
#   breakout        - no loop to break out of
#   TypeI           - no significance testing
# =============================================================================

#' Genomic Prediction via Integrated Modelling
#'
#' @description
#' Fits a genomic best linear unbiased prediction (G-BLUP) model and extracts
#' genomic estimated breeding values (GEBVs) for all genotyped lines. Two
#' computational paths are used depending on the relationship between the
#' number of markers and lines:
#'
#' \describe{
#'   \item{\strong{vm path} (markers > lines)}{Constructs the genomic
#'     relationship matrix \eqn{G = XX'/s} (where \eqn{X} is the marker
#'     matrix and \eqn{s} is a scaling constant) and fits
#'     \code{vm(line, G)} as the genetic random effect in ASReml. GEBVs
#'     are the predicted values of this random effect, obtained via
#'     \code{\link[asreml]{predict.asreml}}.}
#'   \item{\strong{mbf path} (lines \eqn{\ge} markers)}{Models marker effects
#'     \eqn{\mathbf{q}} directly as i.i.d. random effects via ASReml's
#'     \code{mbf()} facility. GEBVs are then computed as
#'     \eqn{\hat{\mathbf{g}} = M\hat{\mathbf{q}}}, where \eqn{M} is the
#'     \eqn{(\text{lines} \times \text{markers})} genotype matrix and
#'     \eqn{\hat{\mathbf{q}}} are the marker effect BLUPs. This avoids the
#'     singularity of \eqn{G} when lines equal or exceed markers.}
#' }
#'
#' @param baseModel A converged \code{asreml} model object representing the
#'   null model without any genomic genetic term. Should capture experimental
#'   design effects, spatial terms, and fixed covariates. The genetic line
#'   term should be included in the random formula; it is removed internally
#'   and replaced by the genomic term. If the model has not converged, one
#'   update attempt is made automatically.
#' @param genObj A genotypic data object of class \code{"wgCross"} (from
#'   \code{\link[wgAim]{primeCross}}, for biparental populations) or class
#'   \code{"wgPanel"} (from \code{\link{primePanel}}, for diversity panels).
#'   Must contain per-chromosome genotype matrices accessible at
#'   \code{$geno[[chr]]$imputed.data} (marker mode) or
#'   \code{$geno[[chr]]$interval.data} (interval mode).
#' @param merge.by Character string naming the column present in both the
#'   phenotypic data and \code{genObj} linking lines across datasets.
#' @param fix.lines Logical. If \code{TRUE} (default), lines in the
#'   phenotypic data that are absent from \code{genObj} are accommodated
#'   by adding a fixed \code{Gomit} factor and restricting the genomic
#'   random effect to genotyped lines only. See \code{\link{qtlAim}} for
#'   full details.
#' @param gen.type Character string, either \code{"marker"} (default) or
#'   \code{"interval"}. Determines which genotype scores are used to build
#'   the marker matrix \eqn{M}. \code{"marker"} uses imputed marker
#'   scores (\code{$imputed.data}); \code{"interval"} uses interval
#'   midpoint scores (\code{$interval.data}) and requires a
#'   \code{"wgCross"} class \code{genObj}. For diversity panels
#'   (\code{"wgPanel"} class), \code{gen.type = "marker"} is required.
#' @param force Logical. If \code{FALSE} (default), the vm path is used when
#'   markers exceed lines and the mbf path otherwise (automatic selection).
#'   Set \code{force = TRUE} to always use the mbf path regardless of the
#'   markers-to-lines ratio. The two paths are statistically equivalent.
#' @param trace Logical or character string. If \code{TRUE} (default), ASReml
#'   output is printed to the console. If a file path character string, output
#'   is redirected there.
#' @param \dots Additional arguments passed to \code{update.asreml}, such as
#'   \code{na.action = na.method(x = "include")}.
#'
#' @details
#' \strong{GEBV standard errors under the mbf path:}
#' For the mbf path, GEBV standard errors are computed using a diagonal
#' approximation: \eqn{SE(\hat{g}_i) \approx \sqrt{M_i^2 \cdot
#' PEV(\hat{q})}}, where \eqn{M_i} is row \eqn{i} of the marker matrix
#' and \eqn{PEV(\hat{q})} is the vector of prediction error variances for
#' the marker effects. This approximation treats marker effects as
#' independent, ignoring their covariances, and may underestimate true
#' standard errors somewhat.
#'
#' \strong{Heritability:}
#' Narrow-sense heritability is estimated as
#' \eqn{h^2 = V_g / (V_g + V_e)}, where \eqn{V_g} is the estimated
#' additive genetic variance and \eqn{V_e} is the residual variance. For
#' the vm path, \eqn{V_g = \sigma^2 \cdot \lambda_{vm}} where
#' \eqn{\lambda_{vm}} is the vm variance parameter. For the mbf path,
#' \eqn{V_g = \sigma^2 \cdot \lambda_{mbf} \cdot \overline{\sum_j m_{ij}^2}},
#' where \eqn{\overline{\sum_j m_{ij}^2}} is the mean sum of squared marker
#' scores per line.
#'
#' \strong{Genomic relationship matrix:}
#' The relationship matrix \eqn{G} is stored in the returned object
#' (\code{$GP$rel.matrix}) for both paths and is used by
#' \code{\link{plot.gpAim}} to produce the relatedness heatmap. For the
#' vm path this is the matrix used in model fitting; for the mbf path it
#' is computed as \eqn{MM'/s} from the fitted marker matrix.
#'
#' @return An object of class \code{c("gpAim","asreml")} -- the fitted ASReml
#'   model augmented with a \code{$GP} list containing:
#' \describe{
#'   \item{\code{$gebv}}{A \code{data.frame} with columns for the line
#'     identifier, \code{GEBV} (genomic estimated breeding value), and
#'     \code{SE} (standard error of the GEBV).}
#'   \item{\code{$gen.type}}{The genotype data type used (\code{"marker"} or
#'     \code{"interval"}).}
#'   \item{\code{$path}}{The computational path used (\code{"vm"} or
#'     \code{"mbf"}).}
#'   \item{\code{$var.genetic}}{Estimated additive genetic variance \eqn{V_g}.}
#'   \item{\code{$var.resid}}{Estimated residual variance \eqn{V_e}.}
#'   \item{\code{$heritability}}{Estimated narrow-sense heritability
#'     \eqn{h^2 = V_g/(V_g + V_e)}.}
#'   \item{\code{$n.markers}}{Number of markers used.}
#'   \item{\code{$rel.matrix}}{The genomic relationship matrix \eqn{G},
#'     stored for use in \code{\link{plot.gpAim}}.}
#'   \item{\code{$genetic.term}}{The name of the genetic line identifier
#'     column (as named in the returned \code{$gebv} data frame).}
#' }
#' Additionally, the (possibly modified) phenotypic data frame is assigned
#' to \code{<response>.data} in the calling environment.
#'
#' @seealso \code{\link{print.gpAim}}, \code{\link{summary.gpAim}},
#'   \code{\link{plot.gpAim}}, \code{\link{primePanel}},
#'   \code{\link[wgAim]{primeCross}}, \code{\link{qtlAim}},
#'   \code{\link{gwasAim}}
#'
#' @examples
#' \dontrun{
#' library(asreml)
#'
#' # --- Diversity panel (vm path, markers > lines) ---
#' panel <- primePanel(geno = geno.mat, map = map.df, encoding = "012")
#' panel$pheno$line <- factor(line.ids)
#' base.mod <- asreml(yield ~ rep, random = ~ line, data = pheno.df)
#'
#' gp.fit <- gpAim(base.mod, genObj = panel, merge.by = "line",
#'                 trace = "trace.txt",
#'                 na.action = na.method(x = "include"))
#'
#' print(gp.fit)
#' summary(gp.fit)
#' plot(gp.fit, type = "caterpillar")
#' plot(gp.fit, type = "heatmap")
#' plot(gp.fit, type = "density", prop.select = 0.10)
#'
#' # --- Biparental population (interval mode) ---
#' data(genoRxK, package = "wgAim")
#' genoRxK <- primeCross(genoRxK, impute = "Martinez", id = "Genotype")
#' base.mod2 <- asreml(yld ~ Type, random = ~ Genotype, data = phenoRxK)
#'
#' gp.fit2 <- gpAim(base.mod2, genObj = genoRxK,
#'                  merge.by = "Genotype", gen.type = "interval")
#' }
#'
#' @name gpAim
#' @export
gpAim <- function(baseModel, ...)
    UseMethod("gpAim")

#' @rdname gpAim
#' @exportS3Method
gpAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

#' @rdname gpAim
#' @exportS3Method
gpAim.asreml <- function(baseModel, genObj, merge.by = NULL,
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
    # Inline only the shared parts -- gpAim has no method/selection/breakout
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

    # genObj: accept either 'wgCross' (from primeCross) or 'wgPanel' (from primePanel)
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, c("wgCross", "wgPanel")))
        stop("genObj must be of class \"wgCross\" (from primeCross()) or ",
             "\"wgPanel\" (from primePanel()).")
    if (gen.type == "interval" && inherits(genObj, "wgPanel"))
        stop("gen.type = \"interval\" requires a 'wgCross' object from primeCross(), ",
             "not a 'wgPanel' object. Use gen.type = \"marker\" with primePanel() output.")
    if (gen.type == "interval" && genObj$type == "marker")
        stop("gen.type = \"interval\" requested but genObj contains no interval data.\n",
             "Re-run primeCross() with type = \"interval\".")

    glines <- genObj$pheno[, merge.by]
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
    gd       <- .buildGenoData(genObj, gen.type, glines, plines)
    genoData <- gd$genoData

    fl           <- .fixLines(baseModel, phenoData, genoData, merge.by, plines, fix.lines, ...)
    baseModel    <- fl$baseModel
    phenoData    <- fl$phenoData
    merge.by     <- fl$merge.by
    rterms       <- fl$rterms
    genetic.term <- fl$genetic.term

    # -------------------------------------------------------------------------
    # Phase 3: Fit GP model -- vm or mbf path
    #
    # vm  path (markers > lines): G = XX' is full rank; fit vm(line, G).
    # mbf path (lines >= markers): G is singular; fit marker effects directly
    #                              via mbf(), then GEBVs = M %*% q.hat.
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
        # mbf path: GEBVs = M %*% q.hat
        # Extract marker effect BLUPs (q.hat) from the mbf random coefficients
        mbf.rows <- grep("mbf", rownames(gpModel$coefficients$random))
        q.hat    <- gpModel$coefficients$random[mbf.rows, 1]
        gebvs    <- as.numeric(genoData %*% q.hat)

        # Approximate SE: sqrt( M^2 %*% PEV(q.hat) )
        # Treats marker effects as independent -- ignores their covariance.
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
    # Marker effects for blups plot
    # vm path:  q.hat = trans %*% g.hat  (back-computed from GEBVs)
    # mbf path: q.hat directly from model coefficients
    # -------------------------------------------------------------------------
    if (use.vm) {
        g.hat  <- gebv$GEBV
        q.hat  <- as.numeric(cov.env$trans %*% g.hat)
    }
    # q.hat already set for mbf path above; marker names from genoData cols
    marker.effects <- data.frame(
        marker = colnames(genoData),
        effect = as.numeric(q.hat),
        stringsAsFactors = FALSE
    )

    # -------------------------------------------------------------------------
    # Package results and clean up
    # -------------------------------------------------------------------------
    # Genomic relationship matrix -- stored for heatmap plot.
    # vm path: already computed; mbf path: compute M %*% M' / scale on demand.
    if (use.vm) {
        rel.matrix <- cov.env$relm
    } else {
        tg         <- t(genoData)
        relm.raw   <- crossprod(tg)
        rel.scale2 <- mean(diag(relm.raw))
        rel.matrix <- relm.raw / rel.scale2
    }

    gp.list <- list(
        gebv           = gebv,
        marker.effects = marker.effects,
        gen.type       = gen.type,
        path           = ifelse(use.vm, "vm", "mbf"),
        var.genetic    = as.numeric(var.genetic),
        var.resid      = as.numeric(var.resid),
        heritability   = h2,
        n.markers      = n.markers,
        rel.scale      = if (use.vm) cov.env$scale else rel.scale2,
        rel.matrix     = rel.matrix,
        genetic.term   = genetic.term
    )

    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    gpModel <- .envFix(gpModel, asremlEnv)
    gpModel$GP <- gp.list
    class(gpModel) <- c("gpAim", "asreml")
    gpModel
}
