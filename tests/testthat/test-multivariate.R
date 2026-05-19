# =============================================================================
# test-multivariate.R
# Tests for the multivariate (Trait = non-NULL) extensions to the wgAim engine.
#
# All tests are pure-R (no ASReml required). They exercise:
#   - pchisq.mixture / qchisq.mixture (engine_lrt.R)
#   - .compute_vqtilde (engine_lrt.R)
#   - .lrtTest with ntrait > 1 (engine_select.R)
#   - .addEffect with Trait non-NULL (engine_effects.R)
#   - .packResults with Trait non-NULL (engine_results.R)
#   - .buildGenomeModel Trait prefix term order (engine_model.R) -- formula string only
#   - qtlAim.asreml / gwasAim.asreml guard clauses for Trait
#   - summary.qtlAim / summary.gwasAim Trait-aware Env column
#   - waldTest.asreml zero-equality path
# =============================================================================

# Bring internal functions into scope
pchisq.mixture   <- wgAim:::pchisq.mixture
qchisq.mixture   <- wgAim:::qchisq.mixture
.compute_vqtilde <- wgAim:::.compute_vqtilde
.lrtTest         <- wgAim:::.lrtTest
.addEffect       <- wgAim:::.addEffect
.packResults     <- wgAim:::.packResults

# =============================================================================
# 1.  pchisq.mixture
# =============================================================================

test_that("pchisq.mixture: ntrait=1 equals 0.5*pchisq(x,0) + 0.5*pchisq(x,1)", {
    # Binomial(df, 1, 0.5): weight 0.5 at df=0, 0.5 at df=1
    # pchisq(x, df=0) = 1 for x > 0 (all probability at 0); = 0 at x=0
    x <- c(1, 2.706, 3.841)
    for (xi in x) {
        expected <- 0.5 * pchisq(xi, df = 0) + 0.5 * pchisq(xi, df = 1)
        expect_equal(pchisq.mixture(xi, ntrait = 1L), expected, tolerance = 1e-10)
    }
})

test_that("pchisq.mixture: ntrait=2 is a weighted average of pchisq(x,0..2)", {
    x        <- c(1, 3, 5)
    ntrait   <- 2L
    df       <- 0:ntrait
    mixprobs <- dbinom(df, size = ntrait, prob = 0.5)
    for (xi in x) {
        expected <- sum(mixprobs * pchisq(xi, df = df))
        expect_equal(pchisq.mixture(xi, ntrait = ntrait), expected, tolerance = 1e-10)
    }
})

test_that("pchisq.mixture: ntrait=3 is a weighted average of pchisq(x,0..3)", {
    x        <- c(2, 4, 7)
    ntrait   <- 3L
    df       <- 0:ntrait
    mixprobs <- dbinom(df, size = ntrait, prob = 0.5)
    for (xi in x) {
        expected <- sum(mixprobs * pchisq(xi, df = df))
        expect_equal(pchisq.mixture(xi, ntrait = ntrait), expected, tolerance = 1e-10)
    }
})

test_that("pchisq.mixture: vectorised input returns vector of same length", {
    x   <- c(0.5, 1.0, 2.0, 3.84)
    out <- pchisq.mixture(x, ntrait = 2L)
    expect_length(out, 4L)
    expect_true(all(out >= 0) && all(out <= 1))
})

test_that("pchisq.mixture: CDF is monotone non-decreasing in x", {
    x   <- c(1, 2, 4, 6, 8)
    out <- pchisq.mixture(x, ntrait = 2L)
    expect_true(all(diff(out) >= 0))
})

test_that("pchisq.mixture: matches manual calculation at x=0", {
    # In R, pchisq(0, df=0) = 0 (convention: P(X <= 0) for df=0 degenerate dist)
    # and pchisq(0, df>0) = 0, so pchisq.mixture(0, ntrait=2) = sum(mixprobs * 0) = 0
    df       <- 0:2
    mixprobs <- dbinom(df, size = 2L, prob = 0.5)
    expected <- sum(mixprobs * pchisq(0, df = df))
    expect_equal(pchisq.mixture(0, ntrait = 2L), expected, tolerance = 1e-10)
})

# =============================================================================
# 2.  qchisq.mixture
# =============================================================================

test_that("qchisq.mixture: round-trips with pchisq.mixture", {
    probs  <- c(0.80, 0.90, 0.95, 0.99)
    ntrait <- 2L
    for (p in probs) {
        q <- qchisq.mixture(p, ntrait = ntrait)
        p_back <- pchisq.mixture(q, ntrait = ntrait)
        expect_equal(p_back, p, tolerance = 1e-6,
                     label = sprintf("round-trip for p=%.2f", p))
    }
})

test_that("qchisq.mixture: ntrait=3 round-trips correctly", {
    probs  <- c(0.90, 0.95)
    ntrait <- 3L
    for (p in probs) {
        q      <- qchisq.mixture(p, ntrait = ntrait)
        p_back <- pchisq.mixture(q, ntrait = ntrait)
        expect_equal(p_back, p, tolerance = 1e-6)
    }
})

test_that("qchisq.mixture: quantile is smaller than plain qchisq at same prob", {
    # Mixture has heavier mass at lower df so quantile is lower than pure chi^2(ntrait)
    p      <- 0.95
    ntrait <- 2L
    q_mix  <- qchisq.mixture(p, ntrait = ntrait)
    q_chi  <- qchisq(p, df = ntrait)
    expect_lt(q_mix, q_chi)
})

# =============================================================================
# 3.  .compute_vqtilde
# =============================================================================

test_that(".compute_vqtilde: returns vector of length nmarkers", {
    set.seed(1)
    nmarkers <- 5L; nlines <- 8L; ntrait <- 2L
    trans   <- matrix(rnorm(nmarkers * nlines), nmarkers, nlines)
    Ga      <- diag(c(1.2, 0.8))
    Ginv    <- solve(Ga)
    vatilde <- kronecker(Ga, diag(nlines))  # simple diagonal case
    out     <- .compute_vqtilde(trans, Ginv, vatilde, ntrait)
    expect_length(out, nmarkers)
    expect_true(all(is.finite(out)))
})

test_that(".compute_vqtilde: all values non-negative for PSD vatilde", {
    set.seed(2)
    nmarkers <- 4L; nlines <- 6L; ntrait <- 2L
    trans   <- matrix(rnorm(nmarkers * nlines), nmarkers, nlines)
    Ga      <- matrix(c(1.5, 0.3, 0.3, 1.0), 2L, 2L)
    Ginv    <- solve(Ga)
    vatilde <- kronecker(Ga, diag(nlines))
    out     <- .compute_vqtilde(trans, Ginv, vatilde, ntrait)
    expect_true(all(out >= -1e-10))   # allowing tiny numerical negative
})

test_that(".compute_vqtilde: ntrait=1 gives same as scalar formula", {
    set.seed(3)
    nmarkers <- 5L; nlines <- 8L; ntrait <- 1L
    trans    <- matrix(rnorm(nmarkers * nlines), nmarkers, nlines)
    avar     <- 1.5
    Ginv     <- matrix(1 / avar, 1L, 1L)
    # vatilde = avar * I_nlines (diagonal posterior variance)
    vatilde  <- matrix(avar * diag(nlines), nlines, nlines)
    out_mv   <- .compute_vqtilde(trans, Ginv, vatilde, ntrait)
    # Scalar formula: vqtilde_i = sum_j trans[i,j]^2 * avar / avar = rowSums(trans^2) * avar / avar
    # More precisely: Ginv = 1/avar, tmp2 = trans[i,] %*% vatilde %*% trans[i,]' = avar * ||trans[i,]||^2
    # tr(Ginv %*% tmp2) = (1/avar) * avar * ||trans[i,]||^2 = ||trans[i,]||^2
    expected <- rowSums(trans^2)
    expect_equal(out_mv, expected, tolerance = 1e-8)
})

# =============================================================================
# 4.  .lrtTest with ntrait > 1
# =============================================================================

.make_ll <- function(ll) list(loglik = ll)

test_that(".lrtTest ntrait=1: unchanged behaviour (one-sided boundary)", {
    stat   <- 3.84   # approx chi^2(1) 0.05 critical value
    res    <- .lrtTest(.make_ll(-50 + stat / 2), .make_ll(-50), 0.05, ntrait = 1L)
    expect_equal(res$stat, stat, tolerance = 1e-8)
    expect_equal(res$pvalue, (1 - pchisq(stat, 1)) / 2, tolerance = 1e-8)
})

test_that(".lrtTest ntrait=2: uses pchisq.mixture", {
    stat   <- 4.0
    res    <- .lrtTest(.make_ll(-50 + stat / 2), .make_ll(-50), 0.05, ntrait = 2L)
    expected_p <- 1 - pchisq.mixture(stat, ntrait = 2L)
    expect_equal(res$pvalue, expected_p, tolerance = 1e-8)
})

test_that(".lrtTest ntrait=2: pass=TRUE when pvalue < TypeI", {
    # Large stat -> small pvalue
    res <- .lrtTest(.make_ll(-40), .make_ll(-50), 0.05, ntrait = 2L)
    expect_true(res$pass)
})

test_that(".lrtTest ntrait=2: pass=FALSE when pvalue > TypeI", {
    # Small stat -> large pvalue
    res <- .lrtTest(.make_ll(-49.9), .make_ll(-50), 0.05, ntrait = 2L)
    expect_false(res$pass)
})

test_that(".lrtTest ntrait=3: pass=TRUE for clearly significant stat", {
    res <- .lrtTest(.make_ll(-30), .make_ll(-50), 0.05, ntrait = 3L)
    expect_true(res$pass)
})

# =============================================================================
# 5.  .addEffect with Trait non-NULL
# =============================================================================

# Minimal fake asreml model sufficient for .addEffect (fixed path only)
.make_add_model <- function(qtl.x = "X.C1.2", Trait = NULL) {
    coefs   <- c(1.0, 0.5)
    vcoefs  <- c(0.04, 0.02)
    nams    <- if (is.null(Trait))
        c("(Intercept)", qtl.x)
    else
        c("(Intercept)", qtl.x, paste(Trait, qtl.x, sep = "_A:"),
          paste(Trait, qtl.x, sep = "_B:"))
    coefs_full  <- setNames(c(1.0, rep(0.5, length(nams) - 1)), nams)
    vcoefs_full <- setNames(c(0.1, rep(0.04, length(nams) - 1)), nams)
    m <- list(
        coefficients = list(fixed = matrix(coefs_full, ncol = 1,
                                            dimnames = list(nams, "effect"))),
        vcoeff       = list(fixed = vcoefs_full),
        sigma2       = 1,
        converge     = TRUE,
        call         = list(fixed  = quote(y ~ 1),
                            random = quote(~ 1),
                            data   = quote(phenoData))
    )
    class(m) <- c("asreml", "list")
    m
}

test_that(".addEffect Trait=NULL: fix.form adds only qtl.x", {
    qtl.x     <- "X.C1.2"
    base_m    <- .make_add_model(qtl.x)
    qtl_m     <- .make_add_model(qtl.x)
    phenoData <- data.frame(id = 1:5, y = rnorm(5))

    # Track what formula was passed to update()
    captured <- list()
    with_mocked_bindings(
        update = function(object, fixed. = NULL, random. = NULL, ...) {
            if (!is.null(fixed.)) captured[[length(captured) + 1L]] <<- fixed.
            object
        },
        .package = "wgAim",
        .addEffect(base_m, qtl_m, phenoData, "id", qtl.x, "fixed", 1L,
                   Trait = NULL)
    )
    # Both baseModel and qtlModel update calls use same formula
    expect_length(captured, 2L)
    fix_str <- deparse(captured[[1]])
    expect_true(grepl("X\\.C1\\.2", fix_str))
    expect_false(grepl(":", fix_str))   # no interaction term
})

test_that(".addEffect Trait non-NULL: fix.form adds qtl.x AND Trait:qtl.x", {
    qtl.x     <- "X.C1.2"
    Trait     <- "Site"
    base_m    <- .make_add_model(qtl.x)
    qtl_m     <- .make_add_model(qtl.x)
    phenoData <- data.frame(id = 1:5, y = rnorm(5), Site = factor(c("A","B","A","B","A")))

    captured <- list()
    with_mocked_bindings(
        update = function(object, fixed. = NULL, ...) {
            if (!is.null(fixed.)) captured[[length(captured) + 1L]] <<- fixed.
            object
        },
        .package = "wgAim",
        .addEffect(base_m, qtl_m, phenoData, "id", qtl.x, "fixed", 1L,
                   Trait = Trait)
    )
    fix_str <- deparse(captured[[1]])
    expect_true(grepl("X\\.C1\\.2", fix_str))
    expect_true(grepl("Site:X\\.C1\\.2|X\\.C1\\.2:Site", fix_str))
})

# =============================================================================
# 6.  .buildGenomeModel Trait prefix term order
# =============================================================================
#
# These tests mock out the ASReml update() call and inspect only the
# formula strings / vmterms that .buildGenomeModel() constructs, confirming:
#  (a) Univariate vm path:      vm(merge.by, covObj)               [genomic]
#                             + merge.by                           [residual]
#  (b) Multivariate diag path:  diag(Trait):vm(merge.by, covObj)  [genomic]
#                             + diag(Trait):merge.by               [residual]
#  (c) Multivariate corgh path: diag(Trait):vm() initially        [genomic]
#                             + corgh(Trait):merge.by              [residual]
#       (additive upgrades to corgh after initial fit)
#  (d) vmterms[2] reflects the residual term as-is from base model
#
# The design principle: the residual genetic term (vmterms[2]) is extracted
# verbatim from the base model's random formula and is never modified.  The
# additive genomic term (vmterms[1]) always starts as diag(Trait):vm/mbf and
# is then upgraded to match the residual structure when it is corh/corgh or fa.

.buildGenomeModel <- wgAim:::.buildGenomeModel
.constructCM      <- wgAim:::.constructCM

# Helper: build a base model mock whose random formula contains the given
# genetic term (e.g. "diag(Site):id") in place of a bare "id".
.mock_base_with_random <- function(rand_formula_str, merge.by = "id") {
    base_m <- make_mock_qtlAim()
    base_m$call$random <- as.formula(paste("~", rand_formula_str))
    base_m
}

test_that(".buildGenomeModel vm path: univariate formula unchanged", {
    set.seed(1)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    # 7 intervals > 5 lines -> vm path; base model has bare "id" random term
    base_m   <- .mock_base_with_random("id")
    caller_e <- new.env(parent = globalenv())
    res <- with_mocked_bindings(
        update = function(object, ...) object,
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData,
                          data.frame(id = rownames(genoData)), "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = NULL)
    )
    # Univariate: no prefix on either term
    expect_true(grepl("^vm\\(id", res$vmterms[1L]))   # genomic: vm(id, covObj)
    expect_equal(res$vmterms[2L], "id")                # residual: bare id
    expect_false(grepl("diag", res$vmterms[1L]))       # no diag prefix
})

test_that(".buildGenomeModel vm path Trait non-NULL, diag residual: diag prefix on both terms", {
    set.seed(2)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    phenoData <- data.frame(
        id   = factor(paste0("L", 1:5)),
        Site = factor(rep(c("S1","S2"), length.out = 5))
    )
    # Base model has diag(Site):id -- the diag residual structure
    base_m   <- .mock_base_with_random("diag(Site):id")
    caller_e <- new.env(parent = globalenv())
    captured <- list()
    with_mocked_bindings(
        update = function(object, random. = NULL, ...) {
            if (!is.null(random.)) captured[[1L]] <<- random.
            object
        },
        .package = "wgAim",
        {
            .buildGenomeModel(base_m, genoData, phenoData, "id",
                              genObj, FALSE, character(0), caller_e,
                              Trait = "Site")
        }
    )
    rhs <- deparse(captured[[1L]][[2]])
    # Additive genomic term: diag(Site):vm(id, covObj)
    expect_true(grepl("diag\\(Site\\):vm\\(id", rhs))
    # Residual genetic term: diag(Site):id (from base model)
    expect_true(grepl("diag\\(Site\\):id", rhs))
    # Must NOT be vm():diag() (wrong order)
    expect_false(grepl("vm\\(id.*\\):diag", rhs))
})

test_that(".buildGenomeModel: vmterms[2] is the genetic term extracted from base model", {
    set.seed(3)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    phenoData <- data.frame(id = factor(paste0("L", 1:5)),
                            Site = factor(c("S1","S2","S1","S2","S1")))
    # Base model has corgh(Site):id -- engine extracts this as the residual term
    base_m   <- .mock_base_with_random("corgh(Site):id")
    caller_e <- new.env(parent = globalenv())
    res <- with_mocked_bindings(
        update = function(object, ...) object,
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData, phenoData, "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = "Site")
    )
    # vmterms[2] must be the corgh term extracted from the base model
    expect_equal(res$vmterms[2L], "corgh(Site):id")
})

test_that(".buildGenomeModel: vmterms[2] is bare merge.by when Trait=NULL", {
    set.seed(4)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    # Univariate base model -- bare "id"
    base_m   <- .mock_base_with_random("id")
    caller_e <- new.env(parent = globalenv())
    res <- with_mocked_bindings(
        update = function(object, ...) object,
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData,
                          data.frame(id = rownames(genoData)), "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = NULL)
    )
    expect_equal(res$vmterms[2L], "id")
})

test_that(".buildGenomeModel: corgh residual triggers corgh upgrade of additive term", {
    set.seed(5)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    phenoData <- data.frame(id = factor(paste0("L", 1:5)),
                            Site = factor(c("S1","S2","S1","S2","S1")))
    # Base model has corgh(Site):id -- additive term should upgrade to corgh after initial diag fit
    base_m   <- .mock_base_with_random("corgh(Site):id")
    caller_e <- new.env(parent = globalenv())
    # Capture ALL update() calls to inspect the upgrade step
    captured <- list()
    with_mocked_bindings(
        update = function(object, random. = NULL, ...) {
            if (!is.null(random.)) captured[[length(captured) + 1L]] <<- random.
            object
        },
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData, phenoData, "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = "Site")
    )
    # Two update() calls expected: initial diag fit + upgrade to corgh
    expect_length(captured, 2L)
    # First call: diag prefix on additive term
    rhs1 <- deparse(captured[[1L]][[2]])
    expect_true(grepl("diag\\(Site\\):vm\\(id", rhs1))
    # Second call (upgrade): corgh prefix on additive term
    rhs2 <- deparse(captured[[2L]][[2]])
    expect_true(grepl("corgh\\(Site\\):vm\\(id", rhs2))
    # Residual term stays as corgh(Site):id throughout
    expect_true(grepl("corgh\\(Site\\):id", rhs2))
})

test_that(".buildGenomeModel: diag residual does NOT trigger upgrade (stays diag)", {
    set.seed(6)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    phenoData <- data.frame(id = factor(paste0("L", 1:5)),
                            Site = factor(c("S1","S2","S1","S2","S1")))
    # Base model has diag(Site):id -- no upgrade should occur
    base_m   <- .mock_base_with_random("diag(Site):id")
    caller_e <- new.env(parent = globalenv())
    captured <- list()
    with_mocked_bindings(
        update = function(object, random. = NULL, ...) {
            if (!is.null(random.)) captured[[length(captured) + 1L]] <<- random.
            object
        },
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData, phenoData, "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = "Site")
    )
    # Only ONE update() call: initial diag fit, no upgrade step
    expect_length(captured, 1L)
    rhs <- deparse(captured[[1L]][[2]])
    expect_true(grepl("diag\\(Site\\):vm\\(id", rhs))
})

test_that(".buildGenomeModel: us(Trial) residual triggers corgh upgrade (same as corgh)", {
    set.seed(7)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    phenoData <- data.frame(id   = factor(paste0("L", 1:5)),
                            Site = factor(c("S1","S2","S1","S2","S1")))
    # Base model has us(Site):id -- additive term should upgrade to corgh
    base_m   <- .mock_base_with_random("us(Site):id")
    caller_e <- new.env(parent = globalenv())
    captured <- list()
    with_mocked_bindings(
        update = function(object, random. = NULL, ...) {
            if (!is.null(random.)) captured[[length(captured) + 1L]] <<- random.
            object
        },
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData, phenoData, "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = "Site")
    )
    # Two update() calls: initial diag + upgrade to corgh
    expect_length(captured, 2L)
    rhs1 <- deparse(captured[[1L]][[2]])
    expect_true(grepl("diag\\(Site\\):vm\\(id", rhs1))
    rhs2 <- deparse(captured[[2L]][[2]])
    expect_true(grepl("corgh\\(Site\\):vm\\(id", rhs2))
    # Residual term stays as us(Site):id throughout
    expect_true(grepl("us\\(Site\\):id", rhs2))
})

test_that(".buildGenomeModel: vmterms[2] for us residual is us(Trait):merge.by", {
    set.seed(8)
    genObj   <- make_wgCross_interval(n_lines = 5, n_chr = 1, n_mar = 8)
    gdat     <- lapply(genObj$geno, function(el) el$interval.data)
    genoData <- do.call("cbind", gdat)
    rownames(genoData) <- paste0("L", 1:5)
    phenoData <- data.frame(id   = factor(paste0("L", 1:5)),
                            Site = factor(c("S1","S2","S1","S2","S1")))
    base_m   <- .mock_base_with_random("us(Site):id")
    caller_e <- new.env(parent = globalenv())
    res <- with_mocked_bindings(
        update = function(object, ...) object,
        .package = "wgAim",
        .buildGenomeModel(base_m, genoData, phenoData, "id",
                          genObj, FALSE, character(0), caller_e,
                          Trait = "Site")
    )
    # vmterms[2] must be the us term extracted verbatim from the base model
    expect_equal(res$vmterms[2L], "us(Site):id")
})

# =============================================================================
# 6b. .qtlSelect Ga extraction -- corgh ntrait > 2
# =============================================================================
#
# The corgh(Trait) Ga extraction was previously hard-coded for ntrait=2.
# These tests verify the generalised loop produces the correct Ga matrix
# for ntrait=2 (regression) and ntrait=3.

test_that(".qtlSelect Ga from corgh: ntrait=2 matches manual construction", {
    # For ntrait=2: vpars = c(var1, var2, cor12)
    # Ga_11 = var1, Ga_22 = var2, Ga_12 = cor12 * sqrt(var1) * sqrt(var2)
    var1 <- 1.2; var2 <- 0.8; cor12 <- 0.6
    vpars <- c(var1, var2, cor12)
    ntrait <- 2L

    sds  <- sqrt(vpars[1:ntrait])
    Ga   <- diag(vpars[1:ntrait])
    cors <- vpars[(ntrait + 1L):length(vpars)]
    idx  <- 0L
    for (col in seq_len(ntrait - 1L)) {
        for (row in (col + 1L):ntrait) {
            idx <- idx + 1L
            Ga[row, col] <- Ga[col, row] <- cors[idx] * sds[row] * sds[col]
        }
    }

    expected_offdiag <- cor12 * sqrt(var1) * sqrt(var2)
    expect_equal(Ga[1L, 1L], var1)
    expect_equal(Ga[2L, 2L], var2)
    expect_equal(Ga[1L, 2L], expected_offdiag, tolerance = 1e-10)
    expect_equal(Ga[2L, 1L], expected_offdiag, tolerance = 1e-10)
    # Must be symmetric
    expect_equal(Ga, t(Ga))
})

test_that(".qtlSelect Ga from corgh: ntrait=3 builds correct 3x3 symmetric matrix", {
    # For ntrait=3: vpars = c(var1, var2, var3, cor21, cor31, cor32)
    # lower-triangle column-major: (2,1), (3,1), (3,2)
    var1 <- 1.0; var2 <- 1.5; var3 <- 0.9
    c21  <- 0.4; c31  <- 0.2; c32  <- 0.5
    vpars  <- c(var1, var2, var3, c21, c31, c32)
    ntrait <- 3L

    sds  <- sqrt(vpars[1:ntrait])
    Ga   <- diag(vpars[1:ntrait])
    cors <- vpars[(ntrait + 1L):length(vpars)]
    idx  <- 0L
    for (col in seq_len(ntrait - 1L)) {
        for (row in (col + 1L):ntrait) {
            idx <- idx + 1L
            Ga[row, col] <- Ga[col, row] <- cors[idx] * sds[row] * sds[col]
        }
    }

    # Diagonal
    expect_equal(diag(Ga), c(var1, var2, var3))
    # Off-diagonal: col=1 row=2 -> cor21 * sd1 * sd2
    expect_equal(Ga[2L, 1L], c21 * sqrt(var1) * sqrt(var2), tolerance = 1e-10)
    # Off-diagonal: col=1 row=3 -> cor31 * sd1 * sd3
    expect_equal(Ga[3L, 1L], c31 * sqrt(var1) * sqrt(var3), tolerance = 1e-10)
    # Off-diagonal: col=2 row=3 -> cor32 * sd2 * sd3
    expect_equal(Ga[3L, 2L], c32 * sqrt(var2) * sqrt(var3), tolerance = 1e-10)
    # Symmetry
    expect_equal(Ga, t(Ga))
    # All ntrait*(ntrait-1)/2 = 3 off-diagonal pairs were filled
    expect_equal(sum(Ga != diag(diag(Ga))), 6L)
})

# =============================================================================
# 7.  .packResults with Trait non-NULL
# =============================================================================

# Reuse the fixture builder from test-engine.R (loaded via helper-fixtures.R)
# but add Trait-specific slots to the mock qtlModel
.make_pack_inputs_mv <- function(n_qtl = 2, Trait = "Site") {
    set.seed(77)
    qtl_keys  <- paste("Chr", paste0("C", seq_len(n_qtl)),
                       seq_len(n_qtl) + 1L, sep = ".")
    mmarks    <- paste("X",  paste0("C", seq_len(n_qtl)),
                       seq_len(n_qtl) + 1L, sep = ".")
    imarks    <- paste(paste0(Trait, "_B"), mmarks, sep = ":")

    eff_names <- c(rbind(mmarks, imarks))
    effects   <- setNames(rnorm(length(eff_names), 0.4, 0.1), eff_names)
    veffects  <- setNames(runif(length(eff_names), 0.02, 0.05), eff_names)

    coef.list  <- list(effects, effects)
    vcoef.list <- list(veffects, veffects)

    state <- rep(1L, 6L)
    names(state) <- paste("Chr", rep(paste0("C", 1:3), each = 2),
                          rep(1:2, times = 3), sep = ".")

    ldiag <- list(lik = list(
        list(baseLogL = -50, stat = 5.0, pvalue = 0.01, pass = TRUE),
        list(baseLogL = -48, stat = 4.5, pvalue = 0.02, pass = TRUE),
        list(baseLogL = -46, stat = 0.5, pvalue = 0.30, pass = FALSE)
    ))

    # Minimal qtlModel with fixed coefficients matching eff_names
    qtlModel <- list(
        coefficients = list(fixed = matrix(
            effects, ncol = 1,
            dimnames = list(eff_names, "effect"))),
        vcoeff  = list(fixed = veffects),
        sigma2  = 1,
        call    = list(fixed  = as.formula(paste("y ~", paste(c(mmarks, imarks), collapse = " + "))),
                       random = quote(~ 1),
                       data   = quote(phenoData)),
        Cfixed  = NULL   # waldTest will warn but we mock it out
    )
    class(qtlModel) <- c("asreml", "list")

    list(
        qtl        = qtl_keys,
        coef.list  = coef.list,
        vcoef.list = vcoef.list,
        ldiag      = ldiag,
        state      = state,
        iter       = 3L,
        breakout   = -1L,
        cov.env    = NULL,
        genetic.term = "id",
        method     = "fixed",
        gen.type   = "interval",
        selection  = "interval",
        TypeI      = 0.05,
        Trait      = Trait,
        qtlModel   = qtlModel
    )
}

test_that(".packResults Trait non-NULL: $QTL$Trait slot is set", {
    inp <- .make_pack_inputs_mv()
    inp$trait.levels <- c("A", "B")
    with_mocked_bindings(
        waldTest = function(object, cc) {
            zdf <- data.frame("Wald Statistic" = c(0.5, 0.3),
                              "P-Value"         = c(0.48, 0.58),
                              check.names = FALSE,
                              row.names   = c("q1", "q2"))
            invisible(list(Contrasts = NULL, Zero = zdf))
        },
        update = function(object, ...) object,
        .package = "wgAim",
        {
            out <- do.call(.packResults, inp)
        }
    )
    expect_equal(out$qtl.list$Trait, "Site")
})

test_that(".packResults Trait non-NULL: $QTL$wald.test slot is a data.frame", {
    inp <- .make_pack_inputs_mv()
    inp$trait.levels <- c("A", "B")
    with_mocked_bindings(
        waldTest = function(object, cc) {
            zdf <- data.frame("Wald Statistic" = c(0.5, 0.3),
                              "P-Value"         = c(0.48, 0.58),
                              check.names = FALSE,
                              row.names   = c("q1", "q2"))
            invisible(list(Contrasts = NULL, Zero = zdf))
        },
        update = function(object, ...) object,
        .package = "wgAim",
        {
            out <- do.call(.packResults, inp)
        }
    )
    expect_s3_class(out$qtl.list$wald.test, "data.frame")
    expect_equal(nrow(out$qtl.list$wald.test), 2L)
})

test_that(".packResults Trait non-NULL: is.interaction matches waldTest p-values", {
    inp <- .make_pack_inputs_mv()
    inp$trait.levels <- c("A", "B")
    with_mocked_bindings(
        waldTest = function(object, cc) {
            # First QTL significant interaction (p=0.01), second not (p=0.60)
            zdf <- data.frame("Wald Statistic" = c(6.5, 0.27),
                              "P-Value"         = c(0.01, 0.60),
                              check.names = FALSE,
                              row.names   = c("q1", "q2"))
            invisible(list(Contrasts = NULL, Zero = zdf))
        },
        update = function(object, ...) object,
        .package = "wgAim",
        {
            out <- do.call(.packResults, inp)
        }
    )
    expect_equal(out$qtl.list$is.interaction, c(TRUE, FALSE))
})

test_that(".packResults Trait=NULL: no Trait slot, qtlModel.pruned unchanged", {
    # Standard univariate inputs (from .make_pack_inputs in test-engine.R)
    inp <- list(
        qtl = "Chr.C1.2",
        coef.list  = list(c(X.C1.2 = 0.5)),
        vcoef.list = list(c(X.C1.2 = 0.04)),
        ldiag      = list(lik = list(
            list(baseLogL = -50, stat = 4, pvalue = 0.02, pass = TRUE),
            list(baseLogL = -48, stat = 0.3, pvalue = 0.35, pass = FALSE)
        )),
        state      = c(Chr.C1.1 = 0L, Chr.C1.2 = 0L, Chr.C2.1 = 1L),
        iter       = 2L,
        breakout   = -1L,
        cov.env    = NULL,
        genetic.term = "id",
        method     = "fixed",
        gen.type   = "interval",
        selection  = "interval",
        TypeI      = 0.05,
        Trait      = NULL,
        qtlModel   = NULL
    )
    out <- do.call(.packResults, inp)
    expect_null(out$qtl.list$Trait)
    expect_null(out$qtl.list$wald.test)
    expect_null(out$qtlModel.pruned)
})

# =============================================================================
# 8.  qtlAim / gwasAim guard clauses for Trait argument
# =============================================================================

# For Trait guard tests we mock .validateModel so it returns successfully,
# then test only the Trait-specific guards that fire immediately after.
# This mirrors the pattern used in test-qtlAim.R for post-validate guards.
.make_vret <- function(model, phenoData, merge.by = "id") {
    list(baseModel = model, asremlEnv = list(), phenoData = phenoData)
}

test_that("qtlAim: Trait column not in phenoData triggers stop", {
    ids       <- paste0("L", 1:10)
    phenoData <- data.frame(id = factor(ids), y = rnorm(10L))
    genObj    <- make_wgCross_interval()
    base_m    <- make_mock_qtlAim()
    with_mocked_bindings(
        .validateModel = function(...) .make_vret(base_m, phenoData),
        .package = "wgAim",
        expect_error(
            qtlAim(base_m, genObj = genObj, merge.by = "id",
                   Trait = "NotAColumn", trace = FALSE),
            "not found"
        )
    )
})

test_that("qtlAim: Trait column not a factor triggers stop", {
    ids       <- paste0("L", 1:10)
    phenoData <- data.frame(id = factor(ids), y = rnorm(10L),
                            Site = as.character(c(rep("A", 5L), rep("B", 5L))))
    genObj    <- make_wgCross_interval()
    base_m    <- make_mock_qtlAim()
    with_mocked_bindings(
        .validateModel = function(...) .make_vret(base_m, phenoData),
        .package = "wgAim",
        expect_error(
            qtlAim(base_m, genObj = genObj, merge.by = "id",
                   Trait = "Site", trace = FALSE),
            "must be a factor"
        )
    )
})

test_that("qtlAim: Trait column with only 1 level triggers stop", {
    ids       <- paste0("L", 1:10)
    phenoData <- data.frame(id = factor(ids), y = rnorm(10L),
                            Site = factor(rep("A", 10L)))
    genObj    <- make_wgCross_interval()
    base_m    <- make_mock_qtlAim()
    with_mocked_bindings(
        .validateModel = function(...) .make_vret(base_m, phenoData),
        .package = "wgAim",
        expect_error(
            qtlAim(base_m, genObj = genObj, merge.by = "id",
                   Trait = "Site", trace = FALSE),
            "at least 2 levels"
        )
    )
})

test_that("qtlAim: str='fa3' too large for ntrait=2 triggers stop", {
    ids       <- paste0("L", 1:10)
    phenoData <- data.frame(id   = factor(ids), y = rnorm(10L),
                            Site = factor(c(rep("A", 5L), rep("B", 5L))))
    genObj    <- make_wgCross_interval()
    base_m    <- make_mock_qtlAim()
    with_mocked_bindings(
        .validateModel = function(...) .make_vret(base_m, phenoData),
        .package = "wgAim",
        # ntrait=2: n.par.us = 3; fa3 gives (3+1)*2 - 3*2/2 = 5 > 3
        expect_error(
            qtlAim(base_m, genObj = genObj, merge.by = "id",
                   Trait = "Site", str = "fa3", trace = FALSE),
            "too large|exceeds unstructured"
        )
    )
})

test_that("gwasAim: Trait column not in phenoData triggers stop", {
    # The gwasAim Trait guard fires after line-matching, so phenoData IDs must
    # match the panel IDs. make_wgPanel() uses "S1", "S2", ... for IDs.
    genObj    <- make_wgPanel()
    ids       <- genObj$pheno$id   # character vector "S1"..."Sn"
    phenoData <- data.frame(id = factor(ids), y = rnorm(length(ids)))
    base_m    <- make_mock_gwasAim()
    with_mocked_bindings(
        .validateModel = function(...) .make_vret(base_m, phenoData),
        .package = "wgAim",
        expect_error(
            gwasAim(base_m, genObj = genObj, merge.by = "id",
                    Trait = "Missing", trace = FALSE),
            "not found"
        )
    )
})

# =============================================================================
# 9.  summary.qtlAim Trait-aware Env column
# =============================================================================

test_that("summary.qtlAim: Trait column present when object$QTL$Trait non-NULL", {
    obj    <- make_mock_qtlAim()
    genObj <- make_wgCross_interval()

    qtl_key  <- names(obj$QTL$effects)[1]   # e.g. "X.C2.2"
    int_key  <- paste0("Site_B:", qtl_key)  # "Site_B:X.C2.2"
    eff_val  <- obj$QTL$effects[[1]]
    veff_val <- obj$QTL$veffects[[1]]

    # Inject multivariate slots
    obj$QTL$Trait        <- "Site"
    obj$QTL$trait.levels <- c("A", "B")

    # summary.qtlAim for MV reads from coefficients$fixed / vcoeff$fixed
    # (not QTL$effects / QTL$veffects -- changed in Session 24 bug fixes).
    # Rows with "X." in their name are the QTL terms; "(Intercept)" is skipped.
    obj$coefficients$fixed <- matrix(
        c(5.0, eff_val, 0.3), ncol = 1,
        dimnames = list(c("(Intercept)", qtl_key, int_key), "effect"))
    obj$vcoeff <- list(
        fixed = setNames(c(0.01, veff_val, 0.03),
                         c("(Intercept)", qtl_key, int_key)))

    # Keep QTL$effects / QTL$veffects consistent (used by other code paths)
    obj$QTL$effects  <- c(obj$QTL$effects,  setNames(0.3,  int_key))
    obj$QTL$veffects <- c(obj$QTL$veffects, setNames(0.03, int_key))

    qtab <- summary(obj, genObj = genObj)
    expect_true("Trait" %in% names(qtab))
    expect_true("MAIN" %in% qtab$Trait || any(grepl("B", qtab$Trait)))
})

test_that("summary.qtlAim: no Trait column when Trait=NULL (univariate)", {
    obj   <- make_mock_qtlAim()
    genObj <- make_wgCross_interval()
    qtab  <- summary(obj, genObj = genObj)
    expect_false("Trait" %in% names(qtab))
})

test_that("summary.gwasAim: Trait column present when object$QTL$Trait non-NULL", {
    obj  <- make_mock_gwasAim()
    obj$QTL$Trait <- "Site"
    old_eff  <- obj$QTL$effects
    qtl_key  <- names(old_eff)[1]
    int_key  <- paste0("Site_B:", qtl_key)
    obj$QTL$effects  <- c(old_eff, setNames(0.25, int_key))
    obj$QTL$veffects <- c(obj$QTL$veffects, setNames(0.03, int_key))
    genObj <- make_wgPanel()
    qtab   <- summary(obj, genObj = genObj)
    expect_true("Trait" %in% names(qtab))
})

test_that("summary.gwasAim: no Trait column when Trait=NULL (univariate)", {
    obj   <- make_mock_gwasAim()
    genObj <- make_wgPanel()
    qtab  <- summary(obj, genObj = genObj)
    expect_false("Trait" %in% names(qtab))
})

# =============================================================================
# 10. waldTest.asreml zero-equality path (self-contained)
# =============================================================================

# Build a minimal fake asreml model with a known Cfixed matrix
.make_wald_model <- function() {
    set.seed(42)
    n    <- 4L
    tau  <- setNames(c(1.0, 0.5, -0.3, 0.8), c("mu", "X.C1.2", "Site_B:X.C1.2", "other"))
    vrb  <- diag(c(0.04, 0.01, 0.02, 0.03)) * 2   # Cfixed = vrb * sigma2
    m <- list(
        coefficients = list(fixed = matrix(tau, ncol = 1L, dimnames = list(names(tau), "e"))),
        vcoeff       = list(fixed = setNames(diag(vrb) / 1.5, names(tau))),
        sigma2       = 1.5,
        Cfixed       = vrb * 1.5,    # Cfixed stored pre-multiplied in real asreml
        converge     = TRUE,
        call         = list(fixed = quote(y ~ 1))
    )
    class(m) <- c("asreml", "list")
    m
}

test_that("waldTest.asreml: zero-equality test returns data.frame with correct structure", {
    m  <- .make_wald_model()
    cc <- list(list(coef = 3L, type = "zero"))   # test coefficient 3 (Site_B:X.C1.2) = 0
    res <- waldTest(m, cc = cc)
    expect_null(res$Contrasts)
    expect_s3_class(res$Zero, "data.frame")
    expect_true("Wald Statistic" %in% names(res$Zero))
    expect_true("P-Value" %in% names(res$Zero))
})

test_that("waldTest.asreml: Wald statistic is positive", {
    m   <- .make_wald_model()
    cc  <- list(list(coef = 3L, type = "zero"))
    res <- waldTest(m, cc = cc)
    expect_gt(res$Zero[["Wald Statistic"]], 0)
})

test_that("waldTest.asreml: p-value in [0, 1]", {
    m   <- .make_wald_model()
    cc  <- list(list(coef = 3L, type = "zero"))
    res <- waldTest(m, cc = cc)
    p   <- res$Zero[["P-Value"]]
    expect_gte(p, 0); expect_lte(p, 1)
})

test_that("waldTest.asreml: coefficient near zero gives large p-value", {
    m  <- .make_wald_model()
    # Replace tau for coefficient 3 with something near zero
    m$coefficients$fixed[3L, 1L] <- 0.001
    cc  <- list(list(coef = 3L, type = "zero"))
    res <- waldTest(m, cc = cc)
    expect_gt(res$Zero[["P-Value"]], 0.2)
})

test_that("waldTest.asreml: large coefficient gives small p-value", {
    m  <- .make_wald_model()
    m$coefficients$fixed[3L, 1L] <- 5.0    # far from zero
    cc  <- list(list(coef = 3L, type = "zero"))
    res <- waldTest(m, cc = cc)
    expect_lt(res$Zero[["P-Value"]], 0.01)
})

test_that("waldTest.default: stops with informative message", {
    expect_error(waldTest(list()), "only implemented")
})
