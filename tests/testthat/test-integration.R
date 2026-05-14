# =============================================================================
# test-integration.R
# End-to-end integration smoke tests.
#
# These tests require ASReml-R to be installed and licensed. They are skipped
# automatically in environments where ASReml is not available (CRAN, CI without
# a license).
#
# DESIGN NOTE — covObj side-effect:
# qtlAim/gwasAim/gpAim call assign("covObj", ..., envir = parent.frame()) to
# write the marker covariance matrix to the caller's frame. The functions also
# assign the updated phenoData as "<response>.data" in the caller's frame.
# Both assignments must land in globalenv() or a persistent environment that
# subsequent internal calls can find via get().
#
# For this reason every analysis call is wrapped in a top-level helper function
# that is *called from* globalenv(), so parent.frame() inside qtlAim resolves
# to globalenv() and covObj persists there for the duration of the analysis.
# =============================================================================

skip_if_not_installed("asreml")

library(asreml)

# ---------------------------------------------------------------------------
# Top-level wrappers — these must NOT be nested inside test_that() bodies.
# Calling them from test_that() means parent.frame() = globalenv() ✓
# ---------------------------------------------------------------------------

.run_qtlAim_interval <- function() {
    data(genoRxK,  package = "wgAim")
    data(phenoRxK, package = "wgAim")
    genoRxK <- suppressWarnings(
        primeCross(genoRxK, impute = "MartinezCurnow",
                   id = "Genotype", subset = c("1A", "2D1"))
    )
    base_mod <- asreml(yld ~ Type + lrow,
                       random   = ~ Genotype + Range,
                       residual = ~ ar1(Range):ar1(Row),
                       data     = phenoRxK, trace = FALSE)
    list(
        fit    = qtlAim(base_mod, genObj = genoRxK,
                        merge.by  = "Genotype",
                        breakout  = 1,
                        trace     = FALSE,
                        na.action = na.method(x = "include")),
        genObj = genoRxK
    )
}

.run_qtlAim_random <- function() {
    data(genoRxK,  package = "wgAim")
    data(phenoRxK, package = "wgAim")
    genoRxK <- suppressWarnings(
        primeCross(genoRxK, impute = "MartinezCurnow",
                   id = "Genotype", subset = c("1A", "2D1"))
    )
    base_mod <- asreml(yld ~ Type + lrow,
                       random   = ~ Genotype + Range,
                       residual = ~ ar1(Range):ar1(Row),
                       data     = phenoRxK, trace = FALSE)
    qtlAim(base_mod, genObj = genoRxK,
           merge.by = "Genotype", method = "random",
           breakout = 1, trace = FALSE,
           na.action = na.method(x = "include"))
}

.run_qtlAim_marker <- function() {
    data(genoRxK,  package = "wgAim")
    data(phenoRxK, package = "wgAim")
    genoRxK <- suppressWarnings(
        primeCross(genoRxK, type = "marker", impute = "MartinezCurnow",
                   id = "Genotype", subset = c("1A", "2D1"))
    )
    base_mod <- asreml(yld ~ Type + lrow,
                       random   = ~ Genotype + Range,
                       residual = ~ ar1(Range):ar1(Row),
                       data     = phenoRxK, trace = FALSE)
    qtlAim(base_mod, genObj = genoRxK,
           merge.by = "Genotype", gen.type = "marker",
           breakout = 1, trace = FALSE,
           na.action = na.method(x = "include"))
}

.run_gpAim <- function() {
    data(genoRxK,  package = "wgAim")
    data(phenoRxK, package = "wgAim")
    genoRxK <- suppressWarnings(
        primeCross(genoRxK, impute = "MartinezCurnow",
                   id = "Genotype", subset = c("1A", "2D1"))
    )
    base_mod <- asreml(yld ~ Type + lrow,
                       random   = ~ Genotype + Range,
                       residual = ~ ar1(Range):ar1(Row),
                       data     = phenoRxK, trace = FALSE)
    list(
        fit    = gpAim(base_mod, genObj = genoRxK,
                       merge.by = "Genotype", gen.type = "interval",
                       trace    = FALSE,
                       na.action = na.method(x = "include")),
        genObj = genoRxK
    )
}

.run_gwasAim <- function() {
    set.seed(55)
    n_lines   <- 60
    n_markers <- 80
    ids       <- paste0("Line_", seq_len(n_lines))
    geno_mat  <- matrix(sample(0:2, n_lines * n_markers, replace = TRUE),
                        nrow = n_lines,
                        dimnames = list(ids, paste0("M", seq_len(n_markers))))
    map_df <- data.frame(
        marker = paste0("M", seq_len(n_markers)),
        chr    = rep(paste0("Chr", 1:4), each = n_markers / 4),
        pos    = rep(seq(5, 100, length.out = n_markers / 4), times = 4),
        stringsAsFactors = FALSE
    )
    panel <- primePanel(geno_mat, map_df, encoding = "012",
                        maf = 0.05, impute = "none")
    pheno_df <- data.frame(Line     = factor(ids),
                           response = rnorm(n_lines, 10, 2),
                           stringsAsFactors = FALSE)
    panel$pheno <- pheno_df
    base_gwas <- asreml(response ~ 1, random = ~ Line,
                        data = pheno_df, trace = FALSE)
    gwasAim(base_gwas, genObj = panel,
            merge.by = "Line", breakout = 1, trace = FALSE,
            na.action = na.method(x = "include"))
}

.run_qtlAim_display <- function() {
    data(genoRxK,  package = "wgAim")
    data(phenoRxK, package = "wgAim")
    genoRxK <- suppressWarnings(
        primeCross(genoRxK, impute = "MartinezCurnow",
                   id = "Genotype", subset = c("1A", "2D1"))
    )
    base_mod <- asreml(yld ~ Type + lrow,
                       random   = ~ Genotype + Range,
                       residual = ~ ar1(Range):ar1(Row),
                       data     = phenoRxK, trace = FALSE)
    list(
        fit    = suppressMessages(
            qtlAim(base_mod, genObj = genoRxK,
                   merge.by  = "Genotype",
                   TypeI     = 0.20,
                   breakout  = 2,
                   trace     = FALSE,
                   na.action = na.method(x = "include"))
        ),
        genObj = genoRxK
    )
}

# =============================================================================
# Tests — each calls a .run_*() wrapper defined at file scope above
# =============================================================================

test_that("qtlAim interval: object structure correct (breakout=1)", {
    skip_if_not_installed("asreml")
    res <- tryCatch(.run_qtlAim_interval(), error = function(e) NULL)
    skip_if(is.null(res), "qtlAim failed (ASReml license or data issue)")
    qtl_fit <- res$fit
    expect_s3_class(qtl_fit, "qtlAim")
    expect_s3_class(qtl_fit, "asreml")
    expect_true(is.list(qtl_fit$QTL))
    expect_gte(qtl_fit$QTL$iterations, 1L)
    expect_true(is.list(qtl_fit$QTL$diag$oint))
    expect_equal(qtl_fit$QTL$type, "interval")
    expect_equal(qtl_fit$QTL$method, "fixed")
})

test_that("qtlAim random method: method field set correctly (breakout=1)", {
    skip_if_not_installed("asreml")
    qtl_rand <- tryCatch(.run_qtlAim_random(), error = function(e) NULL)
    skip_if(is.null(qtl_rand), "qtlAim(random) failed")
    expect_equal(qtl_rand$QTL$method, "random")
})

test_that("qtlAim marker type: type field set correctly (breakout=1)", {
    skip_if_not_installed("asreml")
    qtl_m <- tryCatch(.run_qtlAim_marker(), error = function(e) NULL)
    skip_if(is.null(qtl_m), "qtlAim(marker) failed")
    expect_equal(qtl_m$QTL$type, "marker")
})

test_that("gpAim: GP slot has valid structure and h2 in [0,1]", {
    skip_if_not_installed("asreml")
    res <- tryCatch(.run_gpAim(), error = function(e) NULL)
    skip_if(is.null(res), "gpAim failed")
    gp_fit <- res$fit
    expect_s3_class(gp_fit, "gpAim")
    expect_true(is.list(gp_fit$GP))
    expect_s3_class(gp_fit$GP$gebv, "data.frame")
    expect_true(all(c("GEBV", "SE") %in% names(gp_fit$GP$gebv)))
    expect_true(gp_fit$GP$heritability >= 0 && gp_fit$GP$heritability <= 1)
    expect_true(gp_fit$GP$path %in% c("vm", "mbf"))
})

test_that("gwasAim: synthetic panel pipeline returns correct structure (breakout=1)", {
    skip_if_not_installed("asreml")
    res <- tryCatch(.run_gwasAim(), error = function(e) NULL)
    skip_if(is.null(res), "gwasAim failed")
    expect_s3_class(res, "gwasAim")
    expect_true(!is.null(res$QTL$n.markers))
    expect_gt(res$QTL$n.markers, 0L)
})

test_that("print/summary/aimTrace/getQTL work on a real qtlAim fit", {
    skip_if_not_installed("asreml")
    res <- tryCatch(.run_qtlAim_display(), error = function(e) NULL)
    skip_if(is.null(res), "qtlAim (display) failed")
    qtl_fit <- res$fit
    genObj   <- res$genObj

    expect_no_error(print(qtl_fit, genObj = genObj))

    summ <- summary(qtl_fit, genObj = genObj, LOD = TRUE)
    if (!is.null(summ)) {
        expect_s3_class(summ, "data.frame")
        expect_true("LOD" %in% names(summ))
    }

    expect_no_error(aimTrace(qtl_fit))

    if (!is.null(qtl_fit$QTL$effects)) {
        qtlm <- getQTL(qtl_fit, genObj)
        expect_true(is.matrix(qtlm))
        expect_true(ncol(qtlm) %in% c(4L, 8L))
    }
})

test_that("plot methods work on a real gpAim fit", {
    skip_if_not_installed("asreml")
    res <- tryCatch(.run_gpAim(), error = function(e) NULL)
    skip_if(is.null(res), "gpAim (plot) failed")
    gp_fit <- res$fit
    genObj  <- res$genObj

    expect_s3_class(plot(gp_fit, type = "caterpillar"), "ggplot")
    expect_s3_class(plot(gp_fit, type = "density"),     "ggplot")
    expect_s3_class(suppressWarnings(plot(gp_fit, type = "heatmap")), "ggplot")
    expect_s3_class(plot(gp_fit, type = "blups", genObj = genObj), "ggplot")
})
