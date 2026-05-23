# =============================================================================
# test-aimTrace-extra.R
# Additional tests for aimTrace.qtlAim.R targeting uncovered branches:
#
#   - .build_stability_df: univariate (coefficient matching by name)
#   - .build_stability_df: multivariate (per-trial reconstruction,
#       main-effect reference, interaction term lookup)
#   - .plot_lrt: MV threshold (qchisq.mixture)
#   - .plot_stability: MV facet labels with MAIN/INT suffix
#   - aimTrace.qtlAim: multivariate object (ntrait > 1)
#   - aimTrace.gwasAim: multivariate object (ntrait > 1)
# =============================================================================

skip_if_not_installed("ggplot2")
skip_if_not_installed("ggrepel")

.build_stability_df <- wgAim:::.build_stability_df
.plot_lrt           <- wgAim:::.plot_lrt
.plot_stability     <- wgAim:::.plot_stability

# ---- helper: make a mock MV qtlAim object with Trait slot ------------------
.make_mv_aimTrace_qtlAim <- function(n_qtl = 2, trials = c("S1", "S2"),
                                      is.interaction = c(FALSE, TRUE)) {
    set.seed(77)
    genObj <- make_wgCross_interval(n_lines = 30, n_chr = 2, n_mar = 5)
    chrs   <- names(genObj$geno)
    n_lines <- 30L
    ids     <- paste0("L", seq_len(n_lines))
    nt      <- length(trials)

    chr_idx   <- chrs[seq_len(n_qtl) %% 2 + 1]
    int_idx   <- seq_len(n_qtl) + 1L
    qtl_keys  <- paste("Chr", chr_idx, int_idx, sep = ".")
    x_keys    <- paste("X",   chr_idx, int_idx, sep = ".")

    # Build coef.list: for each iteration k, include effects up to k QTL
    # For MV: main effect + Trial2 interaction for each QTL
    coef.list  <- lapply(seq_len(n_qtl), function(k) {
        mains <- setNames(rnorm(k, 0.4, 0.05), x_keys[seq_len(k)])
        ints  <- setNames(rnorm(k, 0.1, 0.02),
                          paste0(trials[2], ":", x_keys[seq_len(k)]))
        c(mains, ints)
    })
    vcoef.list <- lapply(coef.list, function(ef)
        setNames(runif(length(ef), 0.01, 0.03), names(ef)))

    n_ints_total <- 2 * 5
    all_keys <- paste("Chr", rep(chrs, each = 5),
                      rep(seq_len(5), times = 2), sep = ".")
    state <- rep(1, n_ints_total); names(state) <- all_keys

    oint_list  <- lapply(seq_len(n_qtl), function(k) {
        v <- runif(n_ints_total, 0, 2); names(v) <- all_keys; v
    })
    blup_list  <- lapply(seq_len(n_qtl), function(k) {
        # MV: matrix with ntrait columns
        m <- matrix(rnorm(n_ints_total * nt), n_ints_total, nt,
                    dimnames = list(all_keys, trials))
        m
    })
    lik_list <- lapply(seq_len(n_qtl), function(k)
        list(baseLogL = -50 + k, stat = 4 + k, pvalue = 0.02, pass = TRUE))
    lik.mat <- matrix(
        c(sapply(lik_list, function(l) c(l$baseLogL, l$baseLogL + l$stat/2, l$stat, l$pvalue))),
        ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue"))
    )

    QTL <- list(
        qtl         = qtl_keys,
        effects     = setNames(sapply(coef.list[[n_qtl]], `[`, seq_len(n_qtl)),
                               x_keys),
        veffects    = setNames(rep(0.02, n_qtl), x_keys),
        method      = "fixed",
        type        = "interval",
        selection   = "interval",
        TypeI       = 0.05,
        iterations  = n_qtl + 1L,
        breakout    = FALSE,
        Trait       = "Trial",
        trait.levels = trials,
        is.interaction = is.interaction[seq_len(n_qtl)],
        diag = list(
            oint       = oint_list,
            blups      = blup_list,
            lik        = lik_list,
            ochr       = NULL,
            lik.mat    = lik.mat,
            state      = state,
            genetic.term = "id",
            rel.scale  = 1,
            coef.list  = coef.list,
            vcoef.list = vcoef.list
        )
    )

    vpar <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
    obj <- list(
        converge        = TRUE,
        loglik          = -45,
        sigma2          = 1.0,
        vparameters     = vpar,
        vparameters.con = c(id = 0L, "R!variance" = 0L),
        coefficients    = list(
            fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
            random = matrix(rnorm(n_lines), n_lines, 1,
                            dimnames = list(paste0("id_", ids), "effect"))
        ),
        formulae = list(fixed = yld ~ 1, random = ~ Trial:id),
        mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
        call     = call("asreml"),
        QTL      = QTL
    )
    class(obj) <- c("qtlAim", "asreml")
    obj
}

# =============================================================================
# 1.  .build_stability_df — univariate
# =============================================================================

test_that(".build_stability_df UV: returns data.frame with required columns", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    sdf <- .build_stability_df(obj)
    expect_s3_class(sdf, "data.frame")
    expect_true(all(c("qtl_label", "qtl_idx", "iter", "effect", "se",
                       "lo", "hi", "trial") %in% names(sdf)))
})

test_that(".build_stability_df UV: has one row per (QTL, iter) combination", {
    obj  <- make_mock_qtlAim(n_qtl = 2)
    sdf  <- .build_stability_df(obj)
    # QTL 1 appears in iters 1,2; QTL 2 only in iter 2 → 3 rows
    expect_equal(nrow(sdf), 3L)
})

test_that(".build_stability_df UV: se values are non-negative", {
    obj <- make_mock_qtlAim(n_qtl = 3)
    sdf <- .build_stability_df(obj)
    expect_true(all(sdf$se >= 0))
})

test_that(".build_stability_df UV: lo = effect - se, hi = effect + se", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    sdf <- .build_stability_df(obj)
    expect_equal(sdf$lo, sdf$effect - sdf$se)
    expect_equal(sdf$hi, sdf$effect + sdf$se)
})

# =============================================================================
# 2.  .build_stability_df — multivariate
# =============================================================================

test_that(".build_stability_df MV: returns data.frame with trial column", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    sdf <- .build_stability_df(obj)
    expect_s3_class(sdf, "data.frame")
    expect_true("trial" %in% names(sdf))
})

test_that(".build_stability_df MV: has ntrait rows per (QTL, iter) combination", {
    trials <- c("S1", "S2")
    obj    <- .make_mv_aimTrace_qtlAim(n_qtl = 2, trials = trials)
    sdf    <- .build_stability_df(obj)
    nt     <- length(trials)
    # QTL 1: iters 1,2 → 2*nt rows; QTL 2: iter 2 → 1*nt rows
    expect_equal(nrow(sdf), (2L + 1L) * nt)
})

test_that(".build_stability_df MV: trial values match trait.levels", {
    trials <- c("S1", "S2")
    obj    <- .make_mv_aimTrace_qtlAim(n_qtl = 1, trials = trials)
    sdf    <- .build_stability_df(obj)
    expect_true(all(sdf$trial %in% trials))
})

# =============================================================================
# 3.  .plot_lrt — univariate and multivariate threshold
# =============================================================================

test_that(".plot_lrt UV: returns a ggplot", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(.plot_lrt(obj))
    expect_s3_class(gp, "ggplot")
})

test_that(".plot_lrt MV: returns a ggplot with mixture threshold", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(.plot_lrt(obj))
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 4.  .plot_stability — univariate and multivariate
# =============================================================================

test_that(".plot_stability UV: returns a ggplot", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(.plot_stability(obj))
    expect_s3_class(gp, "ggplot")
})

test_that(".plot_stability MV: returns a ggplot", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(.plot_stability(obj))
    expect_s3_class(gp, "ggplot")
})

test_that(".plot_stability MV: uses per-trial colour scheme", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(.plot_stability(obj))
    # MV plot has a colour scale (legend present for trial)
    expect_false(is.null(gp$scales$scales))
})

# =============================================================================
# 5.  aimTrace.qtlAim — MV object
# =============================================================================

test_that("aimTrace.qtlAim MV: plot=FALSE runs without error", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    expect_no_error(capture.output(aimTrace(obj, plot = FALSE, lik.out = FALSE)))
})

test_that("aimTrace.qtlAim MV: console output contains '(main effect per QTL)'", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    out <- capture.output(aimTrace(obj, plot = FALSE, lik.out = FALSE))
    expect_true(any(grepl("main effect", out, ignore.case = TRUE)))
})

test_that("aimTrace.qtlAim MV: plot='lrt' returns ggplot", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    result <- suppressWarnings(
        capture.output(aimTrace(obj, plot = "lrt", lik.out = FALSE))
    )
    # Re-run capturing the return value
    gp <- suppressWarnings(
        suppressMessages(aimTrace(.make_mv_aimTrace_qtlAim(), plot = "lrt", lik.out = FALSE))
    )
    expect_s3_class(gp, "ggplot")
})

test_that("aimTrace.qtlAim MV: plot='stability' returns ggplot", {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(
        suppressMessages(aimTrace(obj, plot = "stability", lik.out = FALSE))
    )
    expect_s3_class(gp, "ggplot")
})

test_that("aimTrace.qtlAim MV: plot='both' returns named list", {
    obj    <- .make_mv_aimTrace_qtlAim(n_qtl = 2)
    result <- suppressWarnings(
        suppressMessages(aimTrace(obj, plot = "both", lik.out = FALSE))
    )
    expect_type(result, "list")
    expect_named(result, c("lrt", "stability"))
})

# =============================================================================
# 6.  print.out argument — suppresses all console output
# =============================================================================

test_that("aimTrace.qtlAim: print.out=FALSE produces no console output", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    out <- capture.output(aimTrace(obj, plot = FALSE, print.out = FALSE))
    expect_length(out, 0L)
})

test_that("aimTrace.qtlAim: print.out=FALSE still returns plot when requested", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    gp  <- suppressWarnings(aimTrace(obj, plot = "lrt", print.out = FALSE))
    expect_s3_class(gp, "ggplot")
})

test_that("aimTrace.qtlAim: print.out=FALSE suppresses output even with lik.out=TRUE", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    out <- capture.output(aimTrace(obj, plot = FALSE,
                                   print.out = FALSE, lik.out = TRUE))
    expect_length(out, 0L)
})

test_that("aimTrace.qtlAim: print.out=TRUE (default) still prints", {
    obj <- make_mock_qtlAim(n_qtl = 2)
    out <- capture.output(aimTrace(obj, plot = FALSE, lik.out = FALSE))
    expect_true(any(nchar(out) > 0L))
})

# =============================================================================
# 7.  aimTrace.gwasAim — MV object (using gwasAim fixture adapted for MV)
# =============================================================================

.make_mv_aimTrace_gwasAim <- function(n_qtl = 2, trials = c("S1", "S2")) {
    obj <- .make_mv_aimTrace_qtlAim(n_qtl, trials)
    class(obj) <- c("gwasAim", "asreml")
    obj$QTL$n.markers <- 20L
    obj$QTL$type      <- "marker"
    obj
}

test_that("aimTrace.gwasAim MV: plot=FALSE runs without error", {
    obj <- .make_mv_aimTrace_gwasAim(n_qtl = 2)
    expect_no_error(capture.output(aimTrace(obj, plot = FALSE, lik.out = FALSE)))
})

test_that("aimTrace.gwasAim MV: plot='stability' returns ggplot", {
    obj <- .make_mv_aimTrace_gwasAim(n_qtl = 2)
    gp  <- suppressWarnings(
        suppressMessages(aimTrace(obj, plot = "stability", lik.out = FALSE))
    )
    expect_s3_class(gp, "ggplot")
})

test_that("aimTrace.gwasAim: print.out=FALSE produces no console output", {
    obj <- .make_mv_aimTrace_gwasAim(n_qtl = 2)
    out <- capture.output(aimTrace(obj, plot = FALSE, print.out = FALSE))
    expect_length(out, 0L)
})

test_that("aimTrace.gwasAim: print.out=FALSE still returns plot when requested", {
    obj <- .make_mv_aimTrace_gwasAim(n_qtl = 2)
    gp  <- suppressWarnings(aimTrace(obj, plot = "stability", print.out = FALSE))
    expect_s3_class(gp, "ggplot")
})
