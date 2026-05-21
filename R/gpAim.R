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
#'   \code{\link{primeCross}}, for biparental populations) or class
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
#' @param Trait Character string naming a factor column in the phenotypic data
#'   identifying the environment or trial (e.g. \code{"Trial"}). When supplied,
#'   a multivariate G-BLUP model is fitted and per-environment GEBVs are
#'   returned. See \code{str} for how the variance structure of the additive
#'   genomic term is controlled. Default \code{NULL} (univariate).
#' @param str Character string controlling the variance structure applied to
#'   the additive genomic term (\code{vm()} or \code{mbf()}) when
#'   \code{Trait} is non-\code{NULL}. Options:
#'   \describe{
#'     \item{\code{NULL} (default)}{Mirror the structure on the genetic line
#'       term in \code{baseModel}'s random formula: \code{corh}/\code{corgh}/
#'       \code{us} maps to \code{corgh(Trait):vm()}; \code{fa(Trait,k)} maps
#'       to \code{fa(Trait,k):vm()}; \code{diag} stays as
#'       \code{diag(Trait):vm()} (independent per-environment GEBVs).}
#'     \item{\code{"corh"} / \code{"corgh"} / \code{"us"}}{Force
#'       \code{corgh(Trait):vm()}: heterogeneous variances with a common
#'       genetic correlation across environments.}
#'     \item{\code{"fa1"}, \code{"fa2"}, \dots}{Force
#'       \code{fa(Trait,k):vm()} with \eqn{k} factors. The number must be
#'       explicit; \code{str = "fa"} alone returns an error.}
#'     \item{\code{"diag"}}{Force \code{diag(Trait):vm()}: independent
#'       per-environment GEBVs, no borrowing of information across
#'       environments.}
#'   }
#'   Ignored when \code{Trait = NULL}. When \code{str} is supplied, the
#'   structure on the residual genetic term in \code{baseModel} is retained
#'   unchanged; only the additive genomic term is overridden.
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
#'   \code{\link{primeCross}}, \code{\link{qtlAim}},
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
                          force = FALSE, trace = TRUE,
                          Trait = NULL, str = NULL, ...) {

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

    # ------------------------------------------------------------------
    # Multivariate setup
    # ------------------------------------------------------------------
    is.mv        <- !is.null(Trait)
    ntrait       <- 1L
    trait.levels <- NULL
    Ga           <- NULL
    Gcor         <- NULL

    if (is.mv) {
        if (!Trait %in% names(phenoData))
            stop("Trait '", Trait, "' not found in phenotypic data.")
        if (!is.factor(phenoData[[Trait]]))
            phenoData[[Trait]] <- factor(phenoData[[Trait]])
        trait.levels <- levels(phenoData[[Trait]])
        ntrait       <- length(trait.levels)
        if (ntrait < 2L)
            stop("Trait '", Trait, "' must have at least 2 levels for multivariate GP.")
        # Validate str: parse n.fa from "faK" and check it is feasible
        if (!is.null(str)) {
            str.l <- tolower(trimws(str))
            if (str.l == "fa")
                stop("str = \"fa\" requires a number of factors, e.g. str = \"fa2\".")
            if (grepl("^fa[0-9]+$", str.l)) {
                n.fa.str <- as.integer(sub("^fa", "", str.l))
                n.par.fa <- (n.fa.str + 1L) * ntrait - n.fa.str * (n.fa.str - 1L) %/% 2L
                n.par.us <- ntrait * (ntrait + 1L) %/% 2L
                if (n.par.fa > n.par.us)
                    stop("str = \"", str, "\" is too large for ", ntrait,
                         " traits (exceeds unstructured): reduce the number of factors.")
            }
        }
    }

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
    # Phase 3: Fit GP model via .buildGenomeModel()
    #
    # .buildGenomeModel() matches the qtlAim/gwasAim pattern:
    #   - Reads the residual genetic term (containing merge.by) from the base
    #     model random formula and retains it alongside the new genomic vm/mbf
    #     term, so the fitted model is:
    #       random = ~ [prefix]:vm(merge.by, covObj) + [prefix]:merge.by
    #     The residual genetic term captures line-level effects not explained
    #     by the marker panel (polygenic background + ungenotyped lines).
    #   - Handles vm vs mbf path selection (vm when markers > lines, mbf
    #     otherwise, overridden by force = TRUE).
    #   - For MV: the genomic vm term mirrors the variance structure on the
    #     residual genetic term in the base model (corh -> corgh, fa(k) ->
    #     fa(k), diag -> diag).  The user controls the structure by specifying
    #     it on the genetic line term in the base model's random formula.
    #     n.fa controls how many factor-analytic steps to fit when the base
    #     model carries an fa() structure.
    #   - Assigns covObj to caller.env.
    # -------------------------------------------------------------------------
    n.markers    <- ncol(genoData)
    n.lines.geno <- nrow(genoData)
    cat("Fitting Genomic Prediction model...\n")
    gm      <- .buildGenomeModel(baseModel, genoData, phenoData, merge.by,
                                  genObj, force, rterms, caller.env,
                                  Trait = Trait, str = str, ...)
    gpModel <- gm$qtlModel
    cov.env <- gm$cov.env     # NULL for mbf path
    use.vm  <- gm$vm
    vmterm  <- gm$vmterms[1L] # genomic term label; [2] = residual genetic term

    # -------------------------------------------------------------------------
    # Extract GEBVs and compute variance components
    # -------------------------------------------------------------------------
    cat("Extracting GEBVs...\n")
    sigma2 <- gpModel$sigma2
    if (gpModel$vparameters.con[length(gpModel$vparameters.con)] == 4)
        sigma2 <- 1

    if (use.vm) {
        # Internal G.param name for the vm term (works for any variance prefix)
        vm_only <- names(gpModel$G.param)[grep("vm.*covObj", names(gpModel$G.param))]
        if (!is.mv) {
            # Univariate vm: predict GEBVs directly, with SED for Cullis H²
            pv   <- predict(gpModel, classify = vm_only,
                            only = vm_only, vcov = FALSE, sed = TRUE,
                            data = phenoData)
            pvdf    <- pv$pvals
            SE_vec  <- pvdf[["std.error"]]
            var.genetic <- sigma2 *
                gpModel$vparameters[grep("vm.*covObj", names(gpModel$vparameters))]
            var.resid   <- sigma2
            # Mrode (2014) per-BLUP accuracy: sqrt(1 - PEV / Vg)
            acc     <- sqrt(pmax(0, 1 - SE_vec^2 / var.genetic))
            # Cullis (2006) generalised H² from SED upper triangle
            gen.H2  <- .cullis_H2(pv$sed, var.genetic)
            gebv <- data.frame(
                pvdf[[merge.by]],
                pvdf[["predicted.value"]],
                SE_vec,
                acc,
                stringsAsFactors = FALSE
            )
            names(gebv) <- c(genetic.term, "GEBV", "SE", "Accuracy")
            gebv$gen.H2 <- gen.H2
        } else {
            # MV vm: predict per-trial GEBVs via classify = "Trait:merge.by",
            # with SED for per-trial Cullis H²
            classify.term <- paste(Trait, merge.by, sep = ":")
            pv   <- predict(gpModel, classify = classify.term,
                            only = vm_only, vcov = FALSE, sed = TRUE,
                            data = phenoData)
            ord  <- order(pv$pvals[[Trait]], pv$pvals[[merge.by]])
            pvdf <- pv$pvals[ord, ]

            # Extract Ga (ntrait x ntrait genetic covariance matrix)
            Ga          <- .extract_Ga_gp(gpModel, sigma2, ntrait, trait.levels)
            var.genetic <- setNames(diag(Ga), trait.levels)
            var.resid   <- setNames(rep(sigma2, ntrait), trait.levels)

            # Per-(line, trial) Mrode accuracy: sqrt(1 - PEV_ij / G_jj)
            nlines  <- nrow(pvdf) %/% ntrait
            G_diag  <- setNames(diag(Ga), trait.levels)
            SE_vec  <- pvdf[["std.error"]]
            G_jj_vec <- G_diag[as.character(pvdf[[Trait]])]
            acc     <- sqrt(pmax(0, 1 - SE_vec^2 / G_jj_vec))

            # Per-trial Cullis H² from diagonal SED blocks (Trait-major ordering)
            # Rows in pvdf are Trait-major after ordering: trial 1 rows 1:nlines,
            # trial 2 rows (nlines+1):(2*nlines), etc.
            sed_mat <- pv$sed[ord, ord]
            gen.H2_vec <- vapply(seq_len(ntrait), function(j) {
                idx_j <- seq((j - 1L) * nlines + 1L, j * nlines)
                .cullis_H2(sed_mat[idx_j, idx_j, drop = FALSE], G_diag[j])
            }, numeric(1L))
            gen.H2_named <- setNames(gen.H2_vec, trait.levels)

            gebv <- data.frame(
                as.character(pvdf[[merge.by]]),
                factor(as.character(pvdf[[Trait]]), levels = trait.levels),
                pvdf[["predicted.value"]],
                SE_vec,
                acc,
                stringsAsFactors = FALSE
            )
            names(gebv) <- c(genetic.term, Trait, "GEBV", "SE", "Accuracy")
            # Broadcast per-trial gen.H2 to each line row
            gebv$gen.H2 <- gen.H2_named[as.character(pvdf[[Trait]])]

            # Genetic correlation matrix for display
            sds        <- sqrt(pmax(diag(Ga), 0))
            Gcor       <- Ga / outer(sds, sds)
            diag(Gcor) <- 1
            dimnames(Gcor) <- dimnames(Ga)
        }
    } else {
        # Univariate mbf path: GEBVs = M %*% q.hat
        mbf.rows <- grep("mbf", rownames(gpModel$coefficients$random))
        q.hat    <- gpModel$coefficients$random[mbf.rows, 1]
        gebvs    <- as.numeric(genoData %*% q.hat)
        pev      <- sigma2 * gpModel$vcoeff$random[mbf.rows]
        se.gebv  <- sqrt(as.numeric(genoData^2 %*% pev))
        vpar.mbf    <- gpModel$vparameters[grep("mbf", names(gpModel$vparameters))]
        var.genetic <- sigma2 * vpar.mbf * mean(rowSums(genoData^2), na.rm = TRUE)
        var.resid   <- sigma2
        # Mrode accuracy from diagonal PEV approximation (se.gebv is approximate)
        acc.mbf <- sqrt(pmax(0, 1 - se.gebv^2 / var.genetic))
        # Cullis H² requires a SED matrix which is not available on the mbf path
        gen.H2  <- NA_real_
        gebv <- data.frame(
            rownames(genoData), gebvs, se.gebv, acc.mbf,
            stringsAsFactors = FALSE
        )
        names(gebv) <- c(genetic.term, "GEBV", "SE", "Accuracy")
        gebv$gen.H2 <- gen.H2
    }

    h2 <- if (!is.mv) {
        as.numeric(var.genetic / (var.genetic + var.resid))
    } else {
        var.genetic / (var.genetic + var.resid)
    }

    # -------------------------------------------------------------------------
    # Marker effects for blups plot (univariate only)
    # -------------------------------------------------------------------------
    if (use.vm && !is.mv) {
        g.hat  <- gebv$GEBV
        q.hat  <- as.numeric(cov.env$trans %*% g.hat)
    }
    if (!is.mv) {
        marker.effects <- data.frame(
            marker = colnames(genoData),
            effect = as.numeric(q.hat),
            stringsAsFactors = FALSE
        )
    } else {
        marker.effects <- NULL  # per-trial marker effects not computed for MV
    }

    # -------------------------------------------------------------------------
    # Package results and clean up
    # -------------------------------------------------------------------------
    # Genomic relationship matrix -- stored for heatmap plot.
    # vm path: already computed; mbf path: compute M %*% M' / scale on demand.
    if (use.vm) {
        rel.matrix <- cov.env$relm
        rel.scale  <- cov.env$scale
    } else {
        tg         <- t(genoData)
        relm.raw   <- crossprod(tg)
        rel.scale  <- mean(diag(relm.raw))
        rel.matrix <- relm.raw / rel.scale
    }

    # Collect gen.H2: scalar (or NA) for univariate, named vector for MV.
    # gen.H2_named is only defined on the MV vm path; fall back to NULL otherwise.
    gen.H2.stored <- if (is.mv && exists("gen.H2_named")) gen.H2_named
                     else if (exists("gen.H2"))           gen.H2
                     else                                 NA_real_

    gp.list <- list(
        gebv           = gebv,
        marker.effects = marker.effects,
        gen.type       = gen.type,
        path           = ifelse(use.vm, "vm", "mbf"),
        var.genetic    = if (!is.mv) as.numeric(var.genetic) else var.genetic,
        var.resid      = if (!is.mv) as.numeric(var.resid)   else var.resid,
        heritability   = h2,
        gen.H2         = gen.H2.stored,
        n.markers      = n.markers,
        rel.scale      = rel.scale,
        rel.matrix     = rel.matrix,
        genetic.term   = genetic.term,
        # Multivariate slots (NULL for univariate)
        Trait          = Trait,
        trait.levels   = trait.levels,
        Ga             = Ga,
        Gcor           = Gcor
    )

    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    gpModel <- .envFix(gpModel, asremlEnv)
    gpModel$GP <- gp.list
    class(gpModel) <- c("gpAim", "asreml")
    gpModel
}

# =============================================================================
# Helper: Cullis (2006) generalised H² from a symmetric SED sub-matrix.
#
# H² = 1 - mean(upper-triangle SED)² / (2 * G_jj)
#
# G_jj  : genetic variance for the group (scalar, must be > 0).
# sed_sub: square SED matrix for that group's predictions.
#
# Only finite, positive SED values contribute.  Returns NA if insufficient
# data, and clamps to [0, 1].
# =============================================================================
#' @keywords internal
.cullis_H2 <- function(sed_sub, G_jj) {
    if (is.null(sed_sub) || is.na(G_jj) || G_jj <= 0) return(NA_real_)
    ut <- sed_sub[upper.tri(sed_sub)]
    ut <- ut[is.finite(ut) & ut > 0]
    if (!length(ut)) return(NA_real_)
    max(0, min(1, 1 - mean(ut)^2 / (2 * G_jj)))
}

# =============================================================================
# Helper: extract ntrait x ntrait genetic covariance matrix Ga from a fitted
# gpAim model.  Logic mirrors .qtlSelect() in engine_select.R so that Ga
# is consistent regardless of variance structure (diag / corh / fa).
# =============================================================================
#' @keywords internal
.extract_Ga_gp <- function(gpModel, sigma2, ntrait, trait.levels) {
    sterms    <- "vm.*covObj|mbf.*ints"
    gp_names  <- names(gpModel$G.param)
    mterm     <- gp_names[grep(sterms, gp_names)][1L]
    vpar_all  <- gpModel$vparameters
    vpar_nms  <- names(vpar_all)
    vm_idx    <- grep(sterms, vpar_nms)
    vpars.raw <- vpar_all[vm_idx]
    vpar.nms  <- vpar_nms[vm_idx]

    Ga <- if (grepl("diag", mterm)) {
        diag(pmax(vpars.raw * sigma2, 0), ntrait)
    } else if (grepl("corh|corgh", mterm)) {
        is.cor   <- grepl("\\.cor$", vpar.nms)
        var.pars <- pmax(vpars.raw[!is.cor] * sigma2, 0)
        cor.pars <- vpars.raw[is.cor]
        sds      <- sqrt(var.pars)
        G        <- diag(var.pars)
        idx <- 0L
        for (col in seq_len(ntrait - 1L)) {
            for (row in (col + 1L):ntrait) {
                idx <- idx + 1L
                G[row, col] <- G[col, row] <- cor.pars[idx] * sds[row] * sds[col]
            }
        }
        G
    } else if (grepl("fa", mterm)) {
        n.fa.fit   <- as.integer(sub(".*,\\s*(\\d+)\\).*", "\\1", mterm))
        is.loading <- grepl("!fa[0-9]+$", vpar.nms)
        psi <- pmax(vpars.raw[!is.loading] * sigma2, 0)
        Lam <- matrix(vpars.raw[is.loading], ncol = n.fa.fit)
        Lam %*% t(Lam) + diag(as.numeric(psi))
    } else {
        diag(pmax(vpars.raw[seq_len(ntrait)] * sigma2, 0))
    }
    dimnames(Ga) <- list(trait.levels, trait.levels)
    Ga
}
