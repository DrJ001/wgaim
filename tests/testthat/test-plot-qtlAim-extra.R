# =============================================================================
# test-plot-qtlAim-extra.R
# Additional tests for plot.qtlAim.R targeting the MV helper branches and
# any uncovered univariate helpers not hit in test-plot-qtlAim.R.
#
# Targets:
#   - .build_mv_blups_df      (MV blups data construction)
#   - .add_mv_blups_flags     (MV blups vline+text overlay)
#   - .build_mv_effects_df    (MV effects lollipop data)
#   - .plot_mv_effects        (MV effects plot)
#   - .match_perc_var         (effects helper; interval + marker)
#   - .build_effects_df       (univariate effects data)
#   - plot.qtlAim type='effects' with marker-type genObj (UV)
#   - plot.qtlAim MV object: type='blups', type='effects'
#   - plot.gwasAim MV object: type='blups', type='effects'
# =============================================================================

skip_if_not_installed("ggplot2")

# ---- helpers ----------------------------------------------------------------

# Build a mock summary table in the exact column layout that
# summary.qtlAim returns for interval gen.type with a Trait column:
# col1=Trait col2=Chr col3=LMark col4=dist col5=IMark col6=dist(cM) col7=RMark col8=dist col9=Size col10=Prob col11=Perc.Var
.make_mv_interval_summ <- function(qtl_fit, genObj) {
    n_qtl   <- length(qtl_fit$QTL$qtl)
    chrs    <- names(genObj$geno)
    chr_idx <- chrs[seq_len(n_qtl) %% length(chrs) + 1]
    int_idx <- seq_len(n_qtl) + 1L
    dist_val <- vapply(seq_len(n_qtl), function(k)
        as.numeric(genObj$geno[[chr_idx[k]]]$inferred.map[int_idx[k]]),
        numeric(1L))
    data.frame(
        "Env"         = rep("T1", n_qtl),
        Chromosome    = chr_idx,
        "Left Marker" = paste0("mk_L", seq_len(n_qtl)),
        "dist(cM).lh" = dist_val - 5,
        "Infer. Marker" = paste0("int", seq_len(n_qtl)),
        "dist(cM)"    = dist_val,
        "Right Marker" = paste0("mk_R", seq_len(n_qtl)),
        "dist(cM).rh" = dist_val + 5,
        Size          = round(qtl_fit$QTL$effects[seq_len(n_qtl)], 4),
        Prob          = 0.01,
        "Perc.Var"    = c(5.2, 4.1, 3.0)[seq_len(n_qtl)],
        check.names   = FALSE,
        stringsAsFactors = FALSE
    )
}

# Build a mock summary table for GWAS / marker gen.type (with Trait column):
# col1=Trait col2=Chr col3=Marker col4=dist(cM) col5=Size col6=Prob col7=Perc.Var
.make_mv_marker_summ <- function(n_qtl = 2, chrs = c("P1", "P2")) {
    chr_idx <- chrs[seq_len(n_qtl) %% length(chrs) + 1]
    data.frame(
        "Env"       = rep("T1", n_qtl),
        Chromosome  = chr_idx,
        Marker      = paste0("snp", seq_len(n_qtl)),
        "dist(cM)"  = seq(10, 10 * n_qtl, by = 10)[seq_len(n_qtl)],
        Size        = rnorm(n_qtl, 0.3),
        Prob        = 0.02,
        "Perc.Var"  = c(4.0, 3.5)[seq_len(n_qtl)],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

.build_mv_blups_df   <- wgAim:::.build_mv_blups_df
.build_mv_effects_df <- wgAim:::.build_mv_effects_df
.plot_mv_effects     <- wgAim:::.plot_mv_effects
.match_perc_var      <- wgAim:::.match_perc_var
.build_effects_df    <- wgAim:::.build_effects_df
.build_cumpos        <- wgAim:::.build_cumpos

# Build a minimal MV qtlAim mock with $QTL$Trait set and appropriate blups
.make_mv_qtlAim <- function(n_qtl = 2, trials = c("T1", "T2"), n_lines = 30,
                              n_chr = 2, n_mar = 5) {
    set.seed(88)
    genObj <- make_wgCross_interval(n_lines, n_chr, n_mar)
    chrs   <- names(genObj$geno)
    ids    <- paste0("L", seq_len(n_lines))
    nt     <- length(trials)

    chr_idx  <- chrs[seq_len(n_qtl) %% n_chr + 1]
    int_idx  <- seq_len(n_qtl) + 1L
    qtl_keys <- paste("Chr", chr_idx, int_idx, sep = ".")
    x_keys   <- paste("X",   chr_idx, int_idx, sep = ".")

    n_total  <- n_chr * n_mar
    all_keys <- paste("Chr", rep(chrs, each = n_mar),
                      rep(seq_len(n_mar), times = n_chr), sep = ".")

    # blups slot is ntrait-column matrix for MV
    blup_list <- lapply(seq_len(n_qtl), function(k) {
        m <- matrix(rnorm(n_total * nt), n_total, nt,
                    dimnames = list(all_keys, trials))
        m
    })
    oint_list <- lapply(seq_len(n_qtl), function(k) {
        v <- runif(n_total, 0, 2); names(v) <- all_keys; v
    })
    lik_list <- lapply(seq_len(n_qtl), function(k)
        list(baseLogL = -50 + k, stat = 4 + k, pvalue = 0.02, pass = TRUE))
    lik.mat <- matrix(
        c(sapply(lik_list, function(l) c(l$baseLogL, l$baseLogL + l$stat/2, l$stat, l$pvalue))),
        ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue")))

    coef.list  <- lapply(seq_len(n_qtl), function(k) {
        mains <- setNames(rnorm(k, 0.4), x_keys[seq_len(k)])
        ints  <- setNames(rnorm(k, 0.1), paste0(trials[2], ":", x_keys[seq_len(k)]))
        c(mains, ints)
    })
    vcoef.list <- lapply(coef.list, function(ef)
        setNames(runif(length(ef), 0.01, 0.03), names(ef)))

    state <- rep(1, n_total); names(state) <- all_keys

    # Final model coefficients (pruned model) — needed by .build_mv_effects_df
    final_coefs <- coef.list[[n_qtl]]
    coef_mat    <- matrix(final_coefs, ncol = 1,
                          dimnames = list(names(final_coefs), "effect"))
    vcoef_vec   <- vcoef.list[[n_qtl]]

    QTL <- list(
        qtl            = qtl_keys,
        effects        = setNames(rnorm(n_qtl, 0.4), x_keys),
        veffects       = setNames(runif(n_qtl, 0.02, 0.04), x_keys),
        method         = "fixed",
        type           = "interval",
        selection      = "interval",
        TypeI          = 0.05,
        iterations     = n_qtl + 1L,
        breakout       = FALSE,
        Trait          = "Trial",
        trait.levels   = trials,
        is.interaction = rep(FALSE, n_qtl),
        wald.test      = data.frame("Wald Statistic" = rep(0.5, n_qtl),
                                    "P-Value"         = rep(0.48, n_qtl),
                                    check.names = FALSE,
                                    row.names = paste0("q", seq_len(n_qtl))),
        final.terms    = x_keys,
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
        coefficients    = list(fixed  = coef_mat,
                               random = matrix(rnorm(n_lines), n_lines, 1,
                                               dimnames = list(paste0("id_", ids), "effect"))),
        vcoeff          = list(fixed  = vcoef_vec),
        formulae        = list(fixed = yld ~ 1, random = ~ Trial:id),
        mf              = data.frame(id = factor(ids), yld = rnorm(n_lines)),
        call            = call("asreml"),
        QTL             = QTL
    )
    class(obj) <- c("qtlAim", "asreml")
    attr(obj, "genObj") <- genObj
    obj
}

# =============================================================================
# 1.  .build_mv_blups_df
# =============================================================================

test_that(".build_mv_blups_df: returns data.frame with required columns", {
    obj    <- .make_mv_qtlAim()
    genObj <- attr(obj, "genObj")
    chrs   <- names(genObj$geno)
    cp     <- .build_cumpos(genObj, "interval", chrs)

    df <- .build_mv_blups_df(obj, iter = 1L, chr = chrs, cp = cp)
    expect_s3_class(df, "data.frame")
    expect_true(all(c("values", "chr", "dist", "trial") %in% names(df)))
})

test_that(".build_mv_blups_df: has one row per (marker, trial)", {
    obj    <- .make_mv_qtlAim(n_qtl = 1, trials = c("T1", "T2"),
                               n_chr = 2, n_mar = 4)
    genObj <- attr(obj, "genObj")
    chrs   <- names(genObj$geno)
    cp     <- .build_cumpos(genObj, "interval", chrs)
    n_total <- 2 * 4
    nt      <- 2

    df <- .build_mv_blups_df(obj, iter = 1L, chr = chrs, cp = cp)
    expect_equal(nrow(df), n_total * nt)
})

test_that(".build_mv_blups_df: trial column contains trait levels", {
    trials <- c("E1", "E2")
    obj    <- .make_mv_qtlAim(n_qtl = 1, trials = trials, n_chr = 1, n_mar = 4)
    genObj <- attr(obj, "genObj")
    chrs   <- names(genObj$geno)
    cp     <- .build_cumpos(genObj, "interval", chrs)

    df <- .build_mv_blups_df(obj, iter = 1L, chr = chrs, cp = cp)
    expect_true(all(df$trial %in% trials))
})

# =============================================================================
# 2.  plot.qtlAim MV — type='blups'
# =============================================================================

test_that("plot.qtlAim MV type='blups': returns ggplot", {
    obj    <- .make_mv_qtlAim()
    genObj <- attr(obj, "genObj")
    gp     <- plot.qtlAim(obj, genObj, type = "blups")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim MV type='blups': iter filter works", {
    obj    <- .make_mv_qtlAim(n_qtl = 2)
    genObj <- attr(obj, "genObj")
    gp     <- plot.qtlAim(obj, genObj, type = "blups", iter = 1L)
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 3.  .build_mv_effects_df
# =============================================================================

test_that(".build_mv_effects_df: returns data.frame with required columns", {
    obj    <- .make_mv_qtlAim()
    genObj <- attr(obj, "genObj")
    # MV summary has Trait prepended → interval pos at col 6
    # cols: Trait | Chromosome | LMark | dist | IMark | dist(cM) | RMark | dist | Size | Prob | Perc.Var
    mock_summ <- .make_mv_interval_summ(obj, genObj)
    edf <- with_mocked_bindings(
        `summary.qtlAim` = function(object, genObj, ...) mock_summ,
        .package = "wgAim",
        { .build_mv_effects_df(obj, genObj) }
    )
    expect_s3_class(edf, "data.frame")
    expect_true(all(c("effect", "se", "y_row", "qtl_label") %in% names(edf)))
})

# =============================================================================
# 4.  plot.qtlAim MV — type='effects'
# =============================================================================

test_that("plot.qtlAim MV type='effects': returns ggplot", {
    obj    <- .make_mv_qtlAim()
    genObj <- attr(obj, "genObj")
    mock_summ <- .make_mv_interval_summ(obj, genObj)
    gp <- with_mocked_bindings(
        `summary.qtlAim` = function(object, genObj, ...) mock_summ,
        .package = "wgAim",
        { plot.qtlAim(obj, genObj, type = "effects") }
    )
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 5.  plot.gwasAim MV — type='blups' and type='effects'
# =============================================================================

.make_mv_gwasAim <- function(n_qtl = 2, trials = c("T1", "T2")) {
    obj <- .make_mv_qtlAim(n_qtl, trials)
    class(obj) <- c("gwasAim", "asreml")
    obj$QTL$type      <- "marker"
    obj$QTL$n.markers <- 20L
    # For GWAS genObj should be wgPanel but plot.gwasAim uses the same helpers
    panel <- make_wgPanel(n_lines = 30, n_chr = 2, n_mar = 5)
    attr(obj, "genObj") <- panel

    # Adjust blups to use panel marker keys
    chrs     <- names(panel$geno)
    n_total  <- 2 * 5
    all_keys <- paste("Chr", rep(chrs, each = 5),
                      rep(seq_len(5), times = 2), sep = ".")
    nt <- length(trials)
    obj$QTL$diag$blups <- lapply(seq_len(n_qtl), function(k) {
        m <- matrix(rnorm(n_total * nt), n_total, nt,
                    dimnames = list(all_keys, trials))
        m
    })
    obj$QTL$diag$oint <- lapply(seq_len(n_qtl), function(k) {
        v <- runif(n_total, 0, 2); names(v) <- all_keys; v
    })
    obj$QTL$diag$state <- setNames(rep(1, n_total), all_keys)
    obj
}

test_that("plot.gwasAim MV type='blups': returns ggplot", {
    obj    <- .make_mv_gwasAim()
    genObj <- attr(obj, "genObj")
    gp     <- plot.gwasAim(obj, genObj, type = "blups")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim MV type='effects': returns ggplot", {
    obj    <- .make_mv_gwasAim()
    genObj <- attr(obj, "genObj")
    chrs   <- names(genObj$geno)
    n_qtl  <- length(obj$QTL$qtl)
    mock_summ <- .make_mv_marker_summ(n_qtl = n_qtl, chrs = chrs)

    # getQTL() tries to look up map positions from genObj; mock it so it
    # returns a minimal data.frame matching the qtl_keys in the object
    mock_qtlm <- data.frame(
        Chromosome   = chrs[seq_len(n_qtl) %% length(chrs) + 1],
        Marker       = paste0("snp", seq_len(n_qtl)),
        "dist(cM)"   = seq(10, 10 * n_qtl, by = 10),
        "dist(cM).q" = seq(10, 10 * n_qtl, by = 10),
        check.names  = FALSE,
        stringsAsFactors = FALSE
    )
    gp <- with_mocked_bindings(
        `summary.gwasAim` = function(object, genObj, ...) mock_summ,
        `getQTL`          = function(object, genObj, ...) mock_qtlm,
        .package = "wgAim",
        { plot.gwasAim(obj, genObj, type = "effects") }
    )
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 6.  UV helper: .build_effects_df and .match_perc_var
# =============================================================================

test_that(".build_effects_df: returns data.frame with effect and se columns", {
    obj    <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
    genObj <- attr(obj, "genObj")
    # UV interval summary: chr | LMark | dist | IMark | dist(cM) | RMark | dist | Size | Prob | Perc.Var
    chrs    <- names(genObj$geno)
    n_qtl   <- length(obj$QTL$qtl)
    chr_idx <- chrs[seq_len(n_qtl) %% length(chrs) + 1]
    int_idx <- seq_len(n_qtl) + 1L
    dist_val <- vapply(seq_len(n_qtl), function(k)
        as.numeric(genObj$geno[[chr_idx[k]]]$inferred.map[int_idx[k]]), numeric(1L))
    mock_summ <- data.frame(
        Chromosome     = chr_idx,
        "Left Marker"  = paste0("mk_L", seq_len(n_qtl)),
        "dist(cM).lh"  = dist_val - 5,
        "Infer. Marker"= paste0("int", seq_len(n_qtl)),
        "dist(cM)"     = dist_val,
        "Right Marker" = paste0("mk_R", seq_len(n_qtl)),
        "dist(cM).rh"  = dist_val + 5,
        Size           = round(obj$QTL$effects, 4),
        Prob           = 0.01,
        "Perc.Var"     = c(5.2, 4.1)[seq_len(n_qtl)],
        check.names    = FALSE,
        stringsAsFactors = FALSE
    )
    edf <- with_mocked_bindings(
        `summary.qtlAim` = function(object, genObj, ...) mock_summ,
        .package = "wgAim",
        { .build_effects_df(obj, genObj) }
    )
    expect_s3_class(edf, "data.frame")
    expect_true(all(c("effect", "se") %in% names(edf)))
})

test_that(".match_perc_var: returns numeric vector of length n_qtl", {
    obj    <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
    genObj <- attr(obj, "genObj")
    chrs   <- names(genObj$geno)
    n_qtl  <- length(obj$QTL$qtl)
    chr_idx <- chrs[seq_len(n_qtl) %% length(chrs) + 1]
    int_idx <- seq_len(n_qtl) + 1L
    dist_val <- vapply(seq_len(n_qtl), function(k)
        as.numeric(genObj$geno[[chr_idx[k]]]$inferred.map[int_idx[k]]), numeric(1L))
    # UV interval: cols chr | LMark | dist | IMark | dist(cM)=col5 | RMark | dist | Size | Prob | Perc.Var
    mock_summ <- data.frame(
        Chromosome     = chr_idx,
        "Left Marker"  = paste0("mk_L", seq_len(n_qtl)),
        "dist(cM).lh"  = dist_val - 5,
        "Infer. Marker"= paste0("int", seq_len(n_qtl)),
        "dist(cM)"     = dist_val,
        "Right Marker" = paste0("mk_R", seq_len(n_qtl)),
        "dist(cM).rh"  = dist_val + 5,
        Size           = round(obj$QTL$effects, 4),
        Prob           = 0.01,
        "Perc.Var"     = c(5.2, 4.1)[seq_len(n_qtl)],
        check.names    = FALSE
    )
    effects <- obj$QTL$effects
    pv <- suppressWarnings(
        .match_perc_var(mock_summ, obj, genObj, effects, "interval")
    )
    expect_length(pv, length(effects))
})
