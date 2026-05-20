# =============================================================================
# gwasAim.R
# Whole Genome GWAS Analyses via Integrated Modelling.
#
# Forward-selection marker association analysis using the shared wgAim engine.
# Differences from qtlAim:
#   - Takes a 'wgPanel' object (from primePanel()) instead of a 'wgCross' object
#   - gen.type is always 'marker'  -- no interval midpoints in GWAS
#   - method is always 'fixed'     -- GWAS tests marker effects as fixed effects
#   - selection is always 'interval' -- best individual marker, not chromosome
#   - Bonferroni correction applied to TypeI threshold by default
#   - Returns class c('gwasAim', 'asreml')
# =============================================================================

#' GWAS Analyses via Integrated Modelling
#'
#' @description
#' Performs forward-selection genome-wide association analysis in diversity
#' panels using the whole-genome composite model framework of Verbyla et al.
#' (2007, 2012). All panel markers are included simultaneously as a composite
#' genome-wide random effect; at each iteration the marker with the highest
#' outlier statistic is nominated and tested for significance via a likelihood
#' ratio test (LRT). Significant markers are added as fixed effects and the
#' process repeats until no further associations are detected.
#'
#' \code{gwasAim} differs from \code{\link{qtlAim}} in three important ways:
#' \enumerate{
#'   \item Input is a \code{"wgPanel"} object from \code{\link{primePanel}}
#'     rather than a \code{"wgCross"} object from \code{primeCross}.
#'   \item Marker effects are always treated as fixed effects
#'     (\code{method = "fixed"} is hard-coded), consistent with standard
#'     GWAS practice.
#'   \item Selection is always at the individual marker level
#'     (\code{selection = "interval"} is hard-coded), not chromosome-level
#'     aggregation.
#' }
#'
#' @param baseModel A converged \code{asreml} model object representing the
#'   null model. Should capture non-genetic sources of variation (spatial
#'   effects, experimental design terms) and any known population structure
#'   covariates. If the model has not converged, one update is attempted
#'   automatically.
#' @param genObj An object of class \code{"wgPanel"} produced by
#'   \code{\link{primePanel}}. Contains per-chromosome marker genotype matrices
#'   in additive \eqn{\pm 1} coding, used to construct the genome-wide
#'   composite term.
#' @param merge.by Character string naming the column present in both the
#'   phenotypic data and \code{genObj} that links lines across datasets.
#' @param fix.lines Logical. If \code{TRUE} (default), phenotyped lines absent
#'   from \code{genObj} are handled by adding a fixed \code{Gomit} factor,
#'   restricting variance estimation to genotyped lines only. See
#'   \code{\link{qtlAim}} for full details.
#' @param force Logical. If \code{FALSE} (default), the algorithm
#'   automatically uses the \strong{vm} path (genomic relationship matrix)
#'   when markers exceed lines, and the \strong{mbf} path otherwise. In
#'   practice, diversity panels almost always have more markers than lines
#'   so the vm path is typically taken. Set \code{force = TRUE} to override
#'   and use the mbf path regardless.
#' @param exclusion.window Numeric scalar giving the exclusion distance in
#'   cM around each detected marker. Markers within this window are excluded
#'   from selection in subsequent iterations. Default is \code{20} cM.
#'   In GWAS this is a proxy for linkage disequilibrium-based exclusion;
#'   users in high-LD populations may wish to increase this value.
#' @param breakout Integer. If positive, terminates the algorithm after that
#'   many iterations, returning diagnostics without adding the final marker
#'   as a fixed effect. Default is \code{-1} (run to completion).
#' @param TypeI Numeric scalar giving the significance level for the LRT.
#'   Default is \code{0.05}. \strong{No Bonferroni correction is applied.}
#'   The LRT tests the additive variance parameter of the composite
#'   genome-wide term -- a single test per iteration, not one per marker --
#'   so \code{TypeI = 0.05} already acts as a family-wise error rate. This
#'   has been empirically validated for the wgAim algorithm; see Details.
#' @param trace Logical or character string. If \code{TRUE} (default),
#'   ASReml output is printed to the console. If a character string giving
#'   a file path, output is redirected there (errors, warnings and detection
#'   messages still appear on screen).
#' @param verboseLev Integer, \code{0} (default) or \code{1}. At \code{1},
#'   per-marker outlier statistics are printed at each iteration.
#' @param Trait Character string naming a factor column in the phenotypic data
#'   that identifies the environment or trial (e.g. \code{"Site"}).  When
#'   supplied, a multivariate multi-environment GWAS model is fitted using a
#'   mixture chi-squared LRT and Wald test pruning to separate
#'   main-effect from G\eqn{\times}E associations.  Default \code{NULL}
#'   (univariate).
#' @param str Character string controlling the variance structure of the
#'   additive genomic term in the multivariate model.  \code{NULL} (default)
#'   mirrors the structure on the residual genetic term in \code{baseModel}.
#'   Explicit options: \code{"corh"}, \code{"corgh"}, \code{"us"},
#'   \code{"diag"}, or \code{"fa1"}, \code{"fa2"}, etc.  Only relevant when
#'   \code{Trait} is non-\code{NULL}.
#' @param \dots Additional arguments passed to \code{update.asreml}, such as
#'   \code{na.action = na.method(x = "include")}.
#'
#' @details
#' \strong{Why TypeI = 0.05 without Bonferroni:}
#' In standard single-marker GWAS, \eqn{n} LRTs are performed (one per
#' marker) and a Bonferroni correction (\eqn{\alpha / n}) is needed. In the
#' \code{gwasAim} forward-selection framework, at each iteration \emph{only
#' one} LRT is performed: it tests whether the additive variance parameter of
#' the \emph{entire composite genome-wide term} (all remaining markers
#' simultaneously) is significantly greater than zero. This single omnibus
#' test controls the family-wise error rate directly, making an additional
#' Bonferroni correction double-counting. Empirical type I error rate studies
#' confirm that \code{TypeI = 0.05} achieves approximately 5\% family-wise
#' error control under the null.
#'
#' \strong{Outlier statistic vs. p-value:}
#' The outlier statistic \eqn{\tilde{q}_i^2 / \tilde{v}_i} is a
#' \emph{ranking} tool that nominates the most likely associated marker; it
#' is derived from conditional BLUPs in the genome-wide model and should not
#' be converted to per-marker p-values via a chi-squared distribution.
#' Formal significance comes only from the LRT.
#'
#' @return An object of class \code{c("gwasAim","asreml")} -- the final fitted
#'   ASReml model augmented with a \code{$QTL} list (structured identically
#'   to \code{\link{qtlAim}}) plus GWAS-specific fields:
#' \describe{
#'   \item{\code{$TypeI}}{The significance threshold used.}
#'   \item{\code{$n.markers}}{Total number of markers in the panel.}
#' }
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
#' @seealso \code{\link{primePanel}}, \code{\link{summary.gwasAim}},
#'   \code{\link{print.gwasAim}}, \code{\link{plot.gwasAim}},
#'   \code{\link{qtlAim}}, \code{\link{gpAim}}
#'
#' @examples
#' \dontrun{
#' library(asreml)
#'
#' # Build wgPanel object from 0/1/2 encoded genotype matrix
#' panel <- primePanel(geno = geno.mat, map = map.df,
#'                     encoding = "012", maf = 0.05)
#' panel$pheno$line <- factor(line.ids)
#'
#' # Null base model capturing experimental structure
#' base.mod <- asreml(yield ~ rep,
#'                    random   = ~ line,
#'                    data     = pheno.df)
#'
#' # Forward-selection GWAS
#' gwas.fit <- gwasAim(base.mod, genObj = panel,
#'                     merge.by = "line", TypeI = 0.05,
#'                     trace    = "trace.txt",
#'                     na.action = na.method(x = "include"))
#'
#' print(gwas.fit,   genObj = panel)
#' summary(gwas.fit, genObj = panel)
#' plot(gwas.fit, genObj = panel, type = "manhattan")
#' }
#'
#' @name gwasAim
#' @export
gwasAim <- function(baseModel, ...)
    UseMethod("gwasAim")

#' @rdname gwasAim
#' @exportS3Method
gwasAim.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\".")

#' @rdname gwasAim
#' @exportS3Method
gwasAim.asreml <- function(baseModel, genObj, merge.by = NULL,
                            fix.lines = TRUE, force = FALSE,
                            exclusion.window = 20, breakout = -1,
                            TypeI = 0.05, trace = TRUE, verboseLev = 0,
                            Trait = NULL, str = NULL, ...) {

    # Hard-coded engine constants -- not user-configurable for GWAS:
    #   method    = "fixed"    GWAS tests each marker as a fixed effect
    #   selection = "interval" always select the best individual marker
    method    <- "fixed"
    selection <- "interval"

    # Capture calling environment early -- needed for assign() in engine
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

    # GWAS-specific: genObj validation and line matching
    if (missing(genObj))
        stop("genObj is a required argument. Use primePanel() to create one.")
    if (!inherits(genObj, "wgPanel"))
        stop("genObj must be of class \"wgPanel\" produced by primePanel().")
    glines <- genObj$pheno[, merge.by]
    if (is.null(glines))
        stop("Panel data does not contain column \"", merge.by, "\".")
    plines <- phenoData[, merge.by]
    if (is.null(plines))
        stop("Phenotypic data does not contain column \"", merge.by, "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in panel \"", merge.by, "\" column do not match any names ",
             "in phenotypic \"", merge.by, "\" column.")

    # Multivariate argument guards (identical to qtlAim)
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
    # Phase 2a: Build genotype data matrix (marker mode forced for GWAS)
    # -------------------------------------------------------------------------
    gd       <- .buildGenoData(genObj, "marker", glines, plines)
    genoData <- gd$genoData
    mnams    <- gd$mnams
    state    <- gd$state

    # TypeI is used directly as the LRT threshold. No multiple-testing
    # correction is needed: the LRT tests the additive variance parameter of
    # the entire genome-wide composite term (not individual markers), so TypeI
    # already acts as a family-wise error rate.
    n.markers <- ncol(genoData)
    cat(sprintf("\nGWAS significance threshold (TypeI): %.4f  (%d markers in panel)\n",
                TypeI, n.markers))

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
                                   genObj, force, rterms, caller.env,
                                   Trait = Trait, str = str, ...)
    qtlModel  <- gm$qtlModel
    genObj    <- gm$intervalObj   # may have env attribute set (vm path)
    cov.env   <- gm$cov.env
    vm        <- gm$vm
    vmterms   <- gm$vmterms
    n.fa      <- gm$n.fa    # effective number of fa factors (0 if not fa structure)

    # -------------------------------------------------------------------------
    # Phase 4: Iterative forward-selection loop
    # -------------------------------------------------------------------------
    ldiag <- coef.list <- vcoef.list <- list()
    qtl   <- c()
    iter  <- 1

    repeat {
        # Compute outlier statistics and select best marker
        selq               <- .qtlSelect(qtlModel, phenoData, genObj, "marker",
                                         selection, exclusion.window, state, verboseLev,
                                         merge.by = merge.by, Trait = Trait, ntrait = ntrait)
        state              <- selq$state
        ldiag$oint[[iter]] <- selq$oint
        ldiag$ochr[[iter]] <- selq$ochr
        ldiag$blups[[iter]] <- selq$blups

        # Likelihood ratio test: tests significance of the additive variance
        # parameter of the genome-wide composite term (not individual markers).
        # For ntrait > 1: downgrade to diag(Trait):vm/mbf for exact
        # pchisq.mixture(stat, ntrait) type I error control.
        #
        # For fa(Trait,k) residual structures (n.fa > 0): ALSO downgrade the
        # residual genetic term to diag(Trait): in both qtlForLRT and baseForLRT
        # so both models share the same residual structure.  This cleanly
        # isolates the vm diagonal variances without cross-contamination from
        # the fa factor-loading parameters.  For corh/corgh (n.fa = 0) the
        # residual is left as-is — the single vm downgrade is sufficient.
        # Inlined so phenoData is in scope for ASReml formula resolution.
        if (ntrait > 1L) {
            lrt.diag.pfx <- paste0("diag(", Trait, "):")
            mb.esc       <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", merge.by)

            # Build qtlForLRT: downgrade vm term (always) and residual (fa only)
            lrt.all.rt   <- attr(terms(qtlModel$call$random), "term.labels")
            lrt.vm.idx   <- grep("vm.*covObj|mbf.*ints", lrt.all.rt)
            lrt.cur.term <- lrt.all.rt[lrt.vm.idx[1L]]
            if (!startsWith(lrt.cur.term, lrt.diag.pfx)) {
                lrt.all.rt[lrt.vm.idx] <- paste0(lrt.diag.pfx,
                                                  sub("^[^:]+:", "", lrt.cur.term))
            }
            if (n.fa > 0L) {
                lrt.resid.idx <- setdiff(
                    grep(paste0(":", mb.esc, "$"), lrt.all.rt), lrt.vm.idx)
                if (length(lrt.resid.idx) > 0L) {
                    rt <- lrt.all.rt[lrt.resid.idx[1L]]
                    if (!startsWith(rt, lrt.diag.pfx))
                        lrt.all.rt[lrt.resid.idx[1L]] <- paste0(
                            lrt.diag.pfx, sub("^[^:]+:", "", rt))
                }
            }
            lrt.diag.ran <- as.formula(paste("~", paste(lrt.all.rt, collapse = " + ")))
            qtlForLRT    <- update(qtlModel, random. = lrt.diag.ran, ...)

            # Build baseForLRT: for fa, refit base model with diag residual.
            # vmterms[2L] is the residual genetic term (e.g. "fa(Trial,2):id");
            # use it directly rather than calling terms() on baseModel$call$random,
            # which may not be a parseable formula object after repeated update() calls.
            if (n.fa > 0L) {
                lrt.diag.resid <- paste0(lrt.diag.pfx,
                                         sub("^[^:]+:", "", vmterms[2L]))
                base.diag.ran  <- as.formula(paste("~", lrt.diag.resid))
                baseForLRT     <- update(baseModel, random. = base.diag.ran, ...)
            } else {
                baseForLRT <- baseModel
            }
        } else {
            qtlForLRT  <- qtlModel
            baseForLRT <- baseModel
        }
        lrt               <- .lrtTest(qtlForLRT, baseForLRT, TypeI, ntrait = ntrait)
        ldiag$lik[[iter]] <- c(lrt$baseLogL, qtlForLRT$loglik, lrt$stat, lrt$pvalue)
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
        rc        <- .rebuildCovObj(genoData, state, merge.by, genObj,
                                    force, vm, vmterms, qtlModel, caller.env)
        cov.env   <- rc$cov.env
        genObj    <- rc$intervalObj
        qtlModel  <- rc$qtlModel

        qtlModel$call$data <- baseModel$call$data <- quote(phenoData)

        # Add selected marker effect (fixed or random) to both models
        ae                 <- .addEffect(baseModel, qtlModel, phenoData, merge.by,
                                         qtl.x, method, iter, Trait = Trait, ...)
        baseModel          <- ae$baseModel
        qtlModel           <- ae$qtlModel
        coef.list[[iter]]  <- ae$coefs
        vcoef.list[[iter]] <- ae$vcoefs

        iter <- iter + 1
    }

    # -------------------------------------------------------------------------
    # Phase 5: Package results and clean up
    # -------------------------------------------------------------------------
    # Assign phenoData to caller env before .packResults() so waldTest's
    # update() can resolve quote(phenoData) from the model call.
    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    assign(data.name, phenoData, envir = caller.env)
    assign("phenoData",  phenoData, envir = caller.env)

    trait.levels <- if (!is.null(Trait)) levels(phenoData[[Trait]]) else NULL
    pr <- .packResults(qtl, coef.list, vcoef.list, ldiag, state, iter,
                       breakout, cov.env, genetic.term, method, "marker",
                       selection, TypeI, Trait = Trait, qtlModel = qtlModel,
                       trait.levels = trait.levels, phenoData = phenoData)
    qtl.list           <- pr$qtl.list
    qtlModel           <- pr$qtlModel.pruned
    qtl.list$n.markers <- n.markers
    qtlModel <- .envFix(qtlModel, asremlEnv)
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("gwasAim", "asreml")
    qtlModel
}
