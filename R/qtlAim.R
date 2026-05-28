# =============================================================================
# qtlAim.R
# Whole Genome QTL Analyses via Integrated Modelling.
#
# User-facing orchestrator for QTL detection using a forward-selection
# iterative approach. All heavy lifting is delegated to the shared engine
# pieces in engine_*.R. This function is responsible only for:
#   - QTL-specific argument validation (genObj must be class "wgCross" from primeCross())
#   - trace/sink setup (on.exit must be in the outermost frame)
#   - Calling engine pieces in the correct order
#   - Assigning results back to the calling environment
# =============================================================================

#' QTL Analyses via Integrated Modelling
#'
#' @description
#' Performs forward-selection QTL detection in biparental mapping populations
#' using the whole-genome average interval mapping (wgAim) algorithm of
#' Verbyla et al. (2007, 2012), implemented through the ASReml-R linear mixed
#' modelling engine. All marker intervals or markers are included simultaneously
#' as a composite genome-wide random effect; at each iteration the interval
#' with the highest outlier statistic is nominated as a candidate QTL and
#' tested for significance via a likelihood ratio test (LRT). Significant QTL
#' are added as fixed or random effects and the process repeats until no
#' further significant QTL are detected.
#'
#' @param baseModel A converged \code{asreml} model object representing the
#'   null model (no QTL effects). Must include the genetic line term in the
#'   random formula. If the model has not converged, a single update attempt
#'   is made automatically before proceeding.
#' @param genObj An object of class \code{"wgCross"} produced by
#'   \code{\link[wgAim]{primeCross}}. Contains per-chromosome interval midpoint
#'   data (\code{$interval.data}) and imputed marker data
#'   (\code{$imputed.data}) used according to \code{gen.type}.
#' @param merge.by Character string naming the column present in both the
#'   phenotypic data (extracted from \code{baseModel}) and \code{genObj}
#'   that links lines across the two datasets. This is typically the genetic
#'   line identifier (e.g. \code{"Genotype"}).
#' @param fix.lines Logical. If \code{TRUE} (default), lines present in the
#'   phenotypic data but absent from \code{genObj} are accommodated by
#'   adding a fixed factor (\code{Gomit}) to the model and restricting the
#'   genetic random effect to genotyped lines only (\code{Gsave}). This
#'   ensures that additive and non-additive genetic variance components are
#'   estimated from genotyped lines only, avoiding confounding with
#'   ungenotyped lines. A harmless warning from ASReml about zero-df terms
#'   may appear and can be ignored.
#' @param gen.type Character string, either \code{"interval"} (default) or
#'   \code{"marker"}. Controls whether interval midpoint scores
#'   (\code{$interval.data}) or imputed marker scores (\code{$imputed.data})
#'   are used to construct the genome-wide covariate matrix. Interval midpoints
#'   are recommended for biparental populations as they better represent the
#'   genetic content between flanking markers.
#' @param method Character string, either \code{"fixed"} (default) or
#'   \code{"random"}. Determines how detected QTL are incorporated into
#'   subsequent iterations. Under \code{"fixed"} each QTL is added as an
#'   additive fixed effect (Verbyla et al. 2007). Under \code{"random"} each
#'   QTL is added as an additive random effect with its own variance component
#'   (Verbyla et al. 2012). The \code{"random"} approach is recommended for
#'   populations with many QTL and limited replication as it borrows strength
#'   across QTL through shrinkage.
#' @param selection Character string, either \code{"interval"} (default) or
#'   \code{"chromosome"}. Under \code{"interval"} the interval with the
#'   globally highest outlier statistic is selected at each iteration. Under
#'   \code{"chromosome"} a chromosome-level aggregated outlier statistic is
#'   first computed and the best chromosome is selected; the interval with the
#'   highest statistic within that chromosome is then chosen. The
#'   \code{"chromosome"} option is known to perform poorly when chromosomes
#'   have few markers and is not recommended in that case.
#' @param force Logical. If \code{FALSE} (default), the algorithm
#'   automatically selects the \strong{vm} path (genomic relationship matrix,
#'   via \code{vm()} in ASReml) when the number of intervals exceeds the
#'   number of genotyped lines, and the \strong{mbf} path (marker-by-file,
#'   via \code{mbf()} in ASReml) otherwise. If \code{TRUE}, the mbf path is
#'   forced regardless of the markers-to-lines ratio. The vm and mbf paths
#'   are statistically equivalent but differ numerically; the mbf path is
#'   more efficient when lines exceed markers.
#' @param exclusion.window Numeric scalar giving the distance in centimorgans
#'   on each side of a detected QTL within which further intervals are excluded
#'   from selection in subsequent iterations. Default is \code{20} cM. This
#'   prevents selecting multiple intervals within a single QTL region in the
#'   same iteration step. Increasing this value leads to more conservative
#'   spacing between detected QTL.
#' @param breakout Integer. If set to a positive integer \eqn{n}, the
#'   algorithm terminates after \eqn{n} iterations regardless of significance,
#'   returning diagnostic information for those \eqn{n} iterations but
#'   \emph{without} adding the \eqn{n}-th QTL as a fixed or random effect.
#'   This is useful for inspecting outlier statistics or BLUPs from early
#'   iterations. Default is \code{-1} (run to completion).
#' @param TypeI Numeric scalar giving the significance level for the
#'   likelihood ratio test used to determine whether a nominated interval is
#'   a significant QTL. Default is \code{0.05}. The LRT tests the additive
#'   variance parameter of the composite genome-wide term (not an individual
#'   interval effect), so this threshold operates as a \emph{family-wise}
#'   error rate without requiring additional multiple-testing correction.
#' @param trace Logical or character string. If \code{TRUE} (default), all
#'   ASReml output is printed to the console. If a character string giving a
#'   file path, output from all ASReml model fits is redirected to that file
#'   (errors, warnings and QTL detection messages still appear on screen).
#'   Using \code{trace = "trace.txt"} is recommended to reduce console clutter
#'   during analyses with many iterations.
#' @param verboseLev Integer, either \code{0} (default) or \code{1}. At
#'   \code{verboseLev = 1} a table of per-interval and per-chromosome outlier
#'   statistics is printed to the console at each iteration, in addition to
#'   the standard ASReml fitting output. Useful for diagnosing the selection
#'   process.
#' @param Trait Character string naming a factor column in the phenotypic data
#'   that defines the grouping structure for multivariate analysis.  The column
#'   may represent trials, environments, sites, treatments, or measured traits
#'   — any factor whose levels define the multivariate response
#'   (e.g. \code{"Trial"}, \code{"Environment"}, \code{"Treatment"}).  When
#'   supplied, a multivariate QTL model is fitted using a mixture chi-squared
#'   LRT and Wald test pruning to separate main-effect from
#'   group-by-genotype interaction QTL.  Default \code{NULL} (univariate).
#' @param str Character string controlling the variance structure of the
#'   additive genomic term in the multivariate model.  \code{NULL} (default)
#'   mirrors the structure on the residual genetic term in \code{baseModel}.
#'   Explicit options: \code{"corh"}, \code{"corgh"}, \code{"us"},
#'   \code{"diag"}, or \code{"fa1"}, \code{"fa2"}, etc. for factor-analytic
#'   structures.  The number of factors must be stated explicitly;
#'   \code{str = "fa"} alone is an error.  Only relevant when
#'   \code{Trait} is non-\code{NULL}.
#' @param \dots Additional arguments passed to all internal
#'   \code{\link[asreml]{update.asreml}} calls, such as
#'   \code{na.action = na.method(x = "include")} or \code{maxit}.
#'
#' @details
#' \strong{Algorithm overview:}
#' \enumerate{
#'   \item All marker intervals are assembled into a genome-wide covariate
#'     matrix \eqn{M} (lines x intervals).
#'   \item An initial genome-wide model is fitted in which all intervals
#'     collectively contribute a composite random genetic term (either
#'     \code{vm(line, G)} or \code{mbf('ints')}).
#'   \item At each iteration, an \emph{outlier statistic}
#'     \eqn{\tilde{q}_i^2 / \tilde{v}_i} is computed for every active interval
#'     \eqn{i}, where \eqn{\tilde{q}_i} is the back-transformed conditional
#'     BLUP of the interval effect and \eqn{\tilde{v}_i} is its posterior
#'     variance. The interval with the maximum outlier statistic is nominated.
#'   \item A likelihood ratio test (LRT) is performed: \eqn{2(\ell_1 -
#'     \ell_0) \sim \chi^2_1}, testing whether the composite genome-wide
#'     additive variance parameter is zero. If the one-sided p-value
#'     \eqn{< TypeI}, the nominated interval is a significant QTL.
#'   \item The selected QTL is added to the model as a fixed or random effect
#'     (according to \code{method}), the composite genome-wide term is
#'     rebuilt excluding the selected interval and any intervals within
#'     \code{exclusion.window} cM, and the process repeats.
#' }
#'
#' \strong{vm vs mbf path:}
#' The vm path builds a genomic relationship matrix \eqn{G = XX'/s} (where
#' \eqn{s} is a scaling constant) and fits \code{vm(line, G)} in ASReml.
#' The mbf path passes the interval matrix directly to ASReml via
#' \code{mbf('ints')} and is equivalent but more efficient when lines \eqn{\ge}
#' markers. Both paths produce identical QTL selections.
#'
#' @return An object of class \code{c("qtlAim","asreml")} -- the final fitted
#'   ASReml model augmented with a \code{$QTL} list containing:
#' \describe{
#'   \item{\code{$qtl}}{Character vector of internal key names of detected QTL
#'     (e.g. \code{"Chr.3B.14"}), in detection order.}
#'   \item{\code{$effects}}{Named numeric vector of estimated QTL additive
#'     effects from the final model.}
#'   \item{\code{$veffects}}{Named numeric vector of estimated variances of
#'     those QTL effects.}
#'   \item{\code{$method}}{The \code{method} argument used (\code{"fixed"} or
#'     \code{"random"}).}
#'   \item{\code{$type}}{The \code{gen.type} used (\code{"interval"} or
#'     \code{"marker"}).}
#'   \item{\code{$selection}}{The \code{selection} argument used.}
#'   \item{\code{$TypeI}}{The significance threshold used for the LRT.}
#'   \item{\code{$iterations}}{Total number of forward-selection iterations
#'     performed.}
#'   \item{\code{$breakout}}{Integer; the \code{breakout} value passed by the
#'     user (\code{-1} means the algorithm ran to completion).}
#'   \item{\code{$Trait}}{The \code{Trait} column name for multivariate models;
#'     \code{NULL} for univariate.}
#'   \item{\code{$trait.levels}}{Character vector of factor levels of
#'     \code{Trait}; \code{NULL} for univariate.}
#'   \item{\code{$is.interaction}}{Logical vector (one entry per detected QTL)
#'     indicating whether each QTL is a G\eqn{\times}E interaction effect
#'     (\code{TRUE}) or a main effect (\code{FALSE}).  \code{NULL} for
#'     univariate models.}
#'   \item{\code{$diag}}{A diagnostic list containing:
#'     \describe{
#'       \item{\code{$oint}}{Per-iteration list of named numeric vectors of
#'         outlier statistics \eqn{\tilde{q}_i^2/\tilde{v}_i} for all active
#'         intervals.}
#'       \item{\code{$ochr}}{Per-iteration list of chromosome-level aggregated
#'         statistics (used when \code{selection = "chromosome"}).}
#'       \item{\code{$blups}}{Per-iteration list of named numeric vectors of
#'         scaled BLUPs \eqn{\tilde{q}_i/\sqrt{|\tilde{v}_i|}}.}
#'       \item{\code{$lik}}{Per-iteration list of LRT components
#'         (\code{$baseLogL}, \code{$stat}, \code{$pvalue}, \code{$pass}).}
#'       \item{\code{$lik.mat}}{Matrix of LRT statistics across iterations,
#'         used by \code{\link{aimTrace}}.}
#'       \item{\code{$coef.list}}{Per-iteration list of named numeric vectors
#'         of cumulative QTL fixed (or random) effect estimates.}
#'       \item{\code{$vcoef.list}}{Per-iteration list of named numeric vectors
#'         of cumulative QTL effect variances.}
#'       \item{\code{$state}}{Named integer vector recording which intervals
#'         have been selected (1) or are still active (0) at the end of the
#'         last iteration.}
#'       \item{\code{$genetic.term}}{Name of the genetic line identifier column
#'         used in the model.}
#'       \item{\code{$rel.scale}}{Numeric scaling constant \eqn{s} used to
#'         normalise the genomic relationship matrix \eqn{G = XX'/s}.}
#'     }
#'   }
#' }
#' Additionally, the updated phenotypic data frame (with QTL genotype columns
#' appended) is assigned to \code{<response>.data} in the calling environment.
#'
#' @references
#' Verbyla, A.P., Cullis, B.R. and Thompson, R. (2007). The analysis of QTL
#' by simultaneous use of the full linkage map. \emph{Theoretical and Applied
#' Genetics}, \bold{116}, 95--111.
#'
#' Verbyla, A.P., Taylor, J.D. and Verbyla, K.L. (2012). RWGAIM: An efficient
#' high-dimensional random whole genome average (QTL) interval mapping approach.
#' \emph{Genetics Research}, \bold{94}, 291--306.
#'
#' @seealso \code{\link{summary.qtlAim}}, \code{\link{print.qtlAim}},
#'   \code{\link{aimTrace}}, \code{\link{linkMap.qtlAim}},
#'   \code{\link{primeCross}},
#'   \code{\link{gwasAim}}, \code{\link{gpAim}}
#'
#' @examples
#' \dontrun{
#' library(asreml)
#' data(phenoRxK, package = "wgAim")
#' data(genoRxK,  package = "wgAim")
#'
#' # Prepare wgCross object
#' genoRxK <- subset(genoRxK, chr = c("1A","2D1","2D2","3B"))
#' genoRxK <- primeCross(genoRxK, impute = "Martinez", id = "Genotype")
#'
#' # Fit null base model
#' base.mod <- asreml(yld ~ Type + lrow,
#'                    random   = ~ Genotype + Range,
#'                    residual = ~ ar1(Range):ar1(Row),
#'                    data     = phenoRxK)
#'
#' # Forward-selection QTL analysis
#' qtl.fit <- qtlAim(base.mod, genObj = genoRxK,
#'                   merge.by = "Genotype",
#'                   trace    = "trace.txt",
#'                   na.action = na.method(x = "include"))
#'
#' # Results
#' print(qtl.fit,   genObj = genoRxK)
#' summary(qtl.fit, genObj = genoRxK)
#' aimTrace(qtl.fit)
#' linkMap(qtl.fit, genObj = genoRxK)
#' }
#'
#' @name qtlAim
#' @export
qtlAim <- function(baseModel, ...)
    UseMethod("qtlAim")

#' @rdname qtlAim
#' @exportS3Method
qtlAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

#' @rdname qtlAim
#' @exportS3Method
qtlAim.asreml <- function(baseModel, genObj, merge.by = NULL, fix.lines = TRUE,
                           gen.type = NULL, method = "fixed",
                           selection = "interval", force = FALSE,
                           exclusion.window = 20, breakout = -1,
                           TypeI = 0.05, trace = TRUE, verboseLev = 0,
                           Trait = NULL, str = NULL, ...) {

    # Capture calling environment early -- needed for assign() calls in engine
    caller.env <- parent.frame()

    # Trace/sink setup must live here so on.exit fires when qtlAim.asreml() returns
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

    # QTL-specific: genObj validation and line matching
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "wgCross"))
        stop("genObj must be of class \"wgCross\" produced by primeCross().")
    glines <- genObj$pheno[, merge.by]
    if (is.null(glines))
        stop("Genotypic data does not contain column \"", merge.by, "\".")
    plines <- phenoData[, merge.by]
    if (is.null(plines))
        stop("Phenotypic data does not contain column \"", merge.by, "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in genotypic \"", merge.by, "\" column do not match any names ",
             "in phenotypic \"", merge.by, "\" column.")

    # Multivariate argument guards
    ntrait <- 1L
    if (!is.null(Trait)) {
        if (!Trait %in% names(phenoData))
            stop("Trait column \"", Trait, "\" not found in the phenotypic data.")
        if (!is.factor(phenoData[[Trait]]))
            stop("Trait column \"", Trait, "\" must be a factor.")
        ntrait <- length(levels(phenoData[[Trait]]))
        if (ntrait < 2L)
            stop("Trait column \"", Trait, "\" must have at least 2 levels.")
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

    # -------------------------------------------------------------------------
    # Phase 2a: Build genotype data matrix
    # -------------------------------------------------------------------------
    if (missing(gen.type) || is.null(gen.type))
        gen.type <- genObj$type
    if (gen.type == "interval" && genObj$type == "marker")
        stop("gen.type = \"interval\" requested but genObj contains no interval data.\n",
             "Re-run primeCross() with type = \"interval\".")
    gd       <- .buildGenoData(genObj, gen.type, glines, plines)
    genoData <- gd$genoData
    mnams    <- gd$mnams
    state    <- gd$state

    # Phase 2b: Handle lines present in phenotypic but absent from genotypic data
    fl           <- .fixLines(baseModel, phenoData, genoData, merge.by, plines, fix.lines, Trait = Trait, ...)
    baseModel    <- fl$baseModel
    phenoData    <- fl$phenoData
    merge.by     <- fl$merge.by
    rterms       <- fl$rterms
    genetic.term <- fl$genetic.term

    # -------------------------------------------------------------------------
    # Phase 3: Build and fit initial genome-wide model (vm or mbf path)
    # -------------------------------------------------------------------------
    gm            <- .buildGenomeModel(baseModel, genoData, phenoData, merge.by,
                                       genObj, force, rterms, caller.env,
                                       Trait = Trait, str = str, ...)
    qtlModel      <- gm$qtlModel
    genObj        <- gm$intervalObj
    cov.env       <- gm$cov.env
    vm            <- gm$vm
    vmterms       <- gm$vmterms
    n.fa          <- gm$n.fa          # FA factors on residual (0 if none)
    vm.struct     <- gm$vm.struct     # upgrade target for additive term (NULL if diag/UV)
    upgrade.to.fa <- gm$upgrade.to.fa # TRUE when vm.struct == "fa"
    gterm         <- gm$gterm         # current additive term name (diag prefix)
    resid.term    <- gm$resid.term    # residual genetic term

    # -------------------------------------------------------------------------
    # Phase 4: Initial LRT + optional upgrade + iterative forward-selection loop
    # -------------------------------------------------------------------------
    ldiag <- coef.list <- vcoef.list <- list()
    qtl   <- c()
    iter  <- 1L

    # Initial LRT on the diag(Trait):vm model: tests whether there is any
    # significant genome-wide additive variance before committing to building
    # an expensive corh/fa structure.  For UV the diag model IS the full model.
    # .buildLRTModels() downgrades the FA residual to diag for MV.
    cat("\nLikelihood Ratio Test of Additive Variance (initial):\n")
    cat("======================================================\n")
    lrt.mods        <- .buildLRTModels(qtlModel, baseModel, ntrait, Trait,
                                        merge.by, n.fa, vmterms, phenoData, ...)
    lrt             <- .lrtTest(lrt.mods$qtlForLRT, lrt.mods$baseForLRT,
                                 TypeI, ntrait = ntrait)
    ldiag$lik[[1L]] <- c(lrt$baseLogL, lrt.mods$qtlForLRT$loglik,
                         lrt$stat, lrt$pvalue)
    message(sprintf("LRT: statistic = %.4f, p-value = %.4f", lrt$stat, lrt$pvalue))

    if (!lrt$pass) {
        message("No significant additive variance detected. No QTL found.")
    } else {
        # Upgrade additive genomic term to corh/fa structure when requested.
        # Only reached when the initial diag LRT has confirmed there is signal.
        if (!is.null(vm.struct)) {
            up       <- .upgradeVmStructure(qtlModel, vm.struct, n.fa,
                                             upgrade.to.fa, Trait, rterms,
                                             resid.term, gterm, vm, vmterms,
                                             phenoData, ...)
            qtlModel <- up$qtlModel
            vmterms  <- up$vmterms
        }

        repeat {
            # Compute outlier statistics and select best interval/marker
            selq                <- .qtlSelect(qtlModel, phenoData, genObj, gen.type,
                                              selection, exclusion.window, state, verboseLev,
                                              merge.by = merge.by, Trait = Trait, ntrait = ntrait)
            state               <- selq$state
            ldiag$oint[[iter]]  <- selq$oint
            ldiag$ochr[[iter]]  <- selq$ochr
            ldiag$blups[[iter]] <- selq$blups

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
            rc          <- .rebuildCovObj(genoData, state, merge.by, genObj,
                                          force, vm, vmterms, qtlModel, caller.env)
            cov.env     <- rc$cov.env
            genObj      <- rc$intervalObj
            qtlModel    <- rc$qtlModel

            qtlModel$call$data <- baseModel$call$data <- quote(phenoData)

            # Add selected effect (fixed or random) to both models
            ae                 <- .addEffect(baseModel, qtlModel, phenoData, merge.by,
                                             qtl.x, method, iter, Trait = Trait, ...)
            baseModel          <- ae$baseModel
            qtlModel           <- ae$qtlModel
            coef.list[[iter]]  <- ae$coefs
            vcoef.list[[iter]] <- ae$vcoefs

            # LRT at bottom: is there still significant genome-wide variance
            # after accounting for all detected QTL?
            cat(sprintf("\nLikelihood Ratio Test of Additive Variance - Iteration ( %d ):\n", iter))
            cat("=================================================================\n")
            lrt.mods               <- .buildLRTModels(qtlModel, baseModel, ntrait, Trait,
                                                       merge.by, n.fa, vmterms, phenoData, ...)
            lrt                    <- .lrtTest(lrt.mods$qtlForLRT, lrt.mods$baseForLRT,
                                               TypeI, ntrait = ntrait)
            ldiag$lik[[iter + 1L]] <- c(lrt$baseLogL, lrt.mods$qtlForLRT$loglik,
                                        lrt$stat, lrt$pvalue)
            lrt.stop <- !lrt$pass | breakout == iter
            message(sprintf("LRT: statistic = %.4f, p-value = %.4f%s",
                            lrt$stat, lrt$pvalue,
                            if (!lrt$pass)      "  -- No further QTL."
                            else if (lrt.stop)  "  -- Breakout."
                            else                ""))
            if (lrt.stop) break
            iter <- iter + 1L
        }
    }   # end else (lrt$pass)

    # -------------------------------------------------------------------------
    # Phase 5: Package results and clean up
    # -------------------------------------------------------------------------
    # Assign phenoData to the caller env NOW, before .packResults(), so that
    # the waldTest update() call inside .packResults() can find phenoData when
    # ASReml resolves quote(phenoData) from the model call.
    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    assign("phenoData",  phenoData, envir = caller.env)

    trait.levels <- if (!is.null(Trait)) levels(phenoData[[Trait]]) else NULL
    pr <- .packResults(qtl, coef.list, vcoef.list, ldiag, state, iter,
                       breakout, cov.env, genetic.term, method, gen.type,
                       selection, TypeI, Trait = Trait, qtlModel = qtlModel,
                       trait.levels = trait.levels, phenoData = phenoData)
    qtl.list <- pr$qtl.list
    qtlModel <- pr$qtlModel.pruned   # may be updated (interaction pruning) or unchanged
    qtlModel <- .envFix(qtlModel, asremlEnv)
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("qtlAim", "asreml")
    qtlModel
}
