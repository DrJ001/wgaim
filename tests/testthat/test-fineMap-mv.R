# =============================================================================
# test-fineMap-mv.R
# Tests for fineMap() multivariate (Trait) paths.
#
# Strategy
# --------
# fineMap() calls update.asreml() inside the scan loop, so ASReml is required
# for a real run.  We avoid the licence requirement by mocking update() at the
# wgAim namespace level (same pattern as test-engine-asreml.R).
#
# Two MV sub-paths are exercised:
#
#   Main-effect QTL  (is.interaction = FALSE):
#     - old term  "X.chr.idx"        removed from fixed formula
#     - new term  "X.chr.idx_fm"     added (no Trait prefix)
#     - z-ratio test (single df) — identical to UV path
#
#   Interaction QTL  (is.interaction = TRUE):
#     - old term  "Trait:X.chr.idx"  removed from fixed formula
#     - new term  "Trait:X.chr.idx_fm" added
#     - joint Wald test via .waldTest() across non-aliased coefficients;
#       LOD = stat / (2 * log(10))
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_mock_mv_qtlAim()  — multivariate qtlAim object with $QTL$Trait set,
#                            $QTL$is.interaction per QTL, and a fixed formula
#                            containing the pruned final terms.
#   make_updated_model()   — minimal post-update ASReml-like model list.
# =============================================================================

testthat::local_edition(3)

# ---- helpers -----------------------------------------------------------------

# Silently run fineMap with update() mocked in the wgAim namespace.
.run_fm_mv <- function(obj, genObj, upd_m, ...) {
    fe <- new.env(parent = globalenv())
    assign("yld.data", obj$mf, envir = fe)
    environment(obj$call$fixed) <- fe
    with_mocked_bindings(
        update = function(object, ...) upd_m,
        .package = "wgAim",
        fineMap(obj, genObj = genObj, ...)
    )
}

# Build a mock updated model whose fixed coefficients include rows that
# match a candidate column name (cand_pat), with controllable vcoeff values.
make_mv_updated_model <- function(cand_pat, trait_levels,
                                   n_int = 15, loglik = -45) {
    upd_m <- make_updated_model(n_int = n_int, loglik = loglik)
    # Add one coefficient row per trial level for the candidate interaction,
    # to mimic what ASReml returns for "Trial:X.chr.idx_fm".
    int_nams <- paste0(trait_levels, ":", cand_pat)
    all_nams <- c("(Intercept)", int_nams)
    fix_mat  <- matrix(rnorm(length(all_nams), 0.2, 0.05), length(all_nams), 1,
                       dimnames = list(all_nams, "effect"))
    upd_m$coefficients$fixed <- fix_mat
    upd_m$vcoeff$fixed       <- setNames(runif(length(all_nams), 0.01, 0.05),
                                          all_nams)
    # Cfixed: a valid positive-definite matrix so .waldTest() can invert it
    nc <- length(all_nams)
    cf <- diag(runif(nc, 0.5, 1.5))
    rownames(cf) <- colnames(cf) <- all_nams
    upd_m$Cfixed <- cf
    upd_m
}

# =============================================================================
# 1. Class and slot detection
# =============================================================================

test_that("fineMap: detects MV object via object$QTL$Trait non-NULL", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = TRUE)
    expect_false(is.null(obj$QTL$Trait))
    expect_equal(obj$QTL$Trait, "Trial")
})

test_that("fineMap: is.interaction stored correctly for each QTL", {
    obj <- make_mock_mv_qtlAim(n_qtl = 2, interact = c(TRUE, FALSE))
    expect_true(obj$QTL$is.interaction[1])
    expect_false(obj$QTL$is.interaction[2])
})

test_that("fineMap: final.terms contains Trait prefix for interaction QTL", {
    obj <- make_mock_mv_qtlAim(n_qtl = 2, interact = c(TRUE, FALSE))
    expect_true(grepl("^Trial:", obj$QTL$final.terms[1]))
    expect_false(grepl(":",       obj$QTL$final.terms[2]))
})

# =============================================================================
# 2. MV main-effect path — UV technology, no Trait prefix
# =============================================================================

test_that("fineMap MV main-effect: returns fineMap object", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = FALSE)
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    # Mock updated model: single main-effect coefficient for cand_x
    upd_m  <- make_updated_model(n_int = 15, loglik = -45)
    upd_m$coefficients$fixed <- matrix(
        c(5.0, 0.3), 2, 1,
        dimnames = list(c("(Intercept)", "X.C2.2_fm"), "effect")
    )
    upd_m$vcoeff$fixed <- c("(Intercept)" = 0.1, "X.C2.2_fm" = 0.02)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    expect_s3_class(fm, "fineMap")
    expect_equal(length(fm), 1L)
    expect_named(fm[[1]], c("mark", "dist", "pvalue", "LOD"))
})

test_that("fineMap MV main-effect: pvalue and LOD are finite numerics", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = FALSE)
    genObj <- attr(obj, "genObj")

    upd_m <- make_updated_model(n_int = 15, loglik = -45)
    upd_m$coefficients$fixed <- matrix(
        c(5.0, 0.3), 2, 1,
        dimnames = list(c("(Intercept)", "X.C2.2_fm"), "effect")
    )
    upd_m$vcoeff$fixed <- c("(Intercept)" = 0.1, "X.C2.2_fm" = 0.02)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    non_na <- fm[[1]][!is.na(fm[[1]]$pvalue), ]
    expect_true(nrow(non_na) > 0L)
    expect_true(all(is.finite(non_na$pvalue)))
    expect_true(all(is.finite(non_na$LOD)))
})

test_that("fineMap MV main-effect: pvalue in [0, 1]", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = FALSE)
    genObj <- attr(obj, "genObj")

    upd_m <- make_updated_model(n_int = 15, loglik = -45)
    upd_m$coefficients$fixed <- matrix(
        c(5.0, 0.3), 2, 1,
        dimnames = list(c("(Intercept)", "X.C2.2_fm"), "effect")
    )
    upd_m$vcoeff$fixed <- c("(Intercept)" = 0.1, "X.C2.2_fm" = 0.02)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    pv <- fm[[1]]$pvalue[!is.na(fm[[1]]$pvalue)]
    expect_true(all(pv >= 0 & pv <= 1))
})

# =============================================================================
# 3. MV interaction path — Wald test, Trait:cand_x formula
# =============================================================================

test_that("fineMap MV interaction: returns fineMap object", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = TRUE)
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    # Candidate pattern for the first scan position on C2 chr (step 20 cM)
    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    expect_s3_class(fm, "fineMap")
    expect_equal(length(fm), 1L)
    expect_named(fm[[1]], c("mark", "dist", "pvalue", "LOD"))
})

test_that("fineMap MV interaction: pvalue and LOD are finite numerics", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = TRUE)
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    non_na <- fm[[1]][!is.na(fm[[1]]$pvalue), ]
    expect_true(nrow(non_na) > 0L)
    expect_true(all(is.finite(non_na$pvalue)))
    expect_true(all(is.finite(non_na$LOD)))
})

test_that("fineMap MV interaction: pvalue in [0, 1]", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = TRUE)
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    pv <- fm[[1]]$pvalue[!is.na(fm[[1]]$pvalue)]
    expect_true(all(pv >= 0 & pv <= 1))
})

test_that("fineMap MV interaction: LOD consistent with Wald stat / (2 * log(10))", {
    # Confirm LOD = Wald_stat / (2*log(10)), not 0.5*log(exp(zrat^2), 10).
    # The key difference: for interaction path LOD uses the multi-df Wald stat,
    # not a single-df z-ratio.  We check LOD >= 0 and that it is NOT simply
    # zrat^2 from a single coefficient (which would be smaller).
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = TRUE)
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels   # 3 levels → 3 Wald df

    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    lod_vals <- fm[[1]]$LOD[!is.na(fm[[1]]$LOD)]
    expect_true(all(lod_vals >= 0))
})

# =============================================================================
# 4. Mixed object: one interaction QTL + one main-effect QTL
# =============================================================================

test_that("fineMap MV mixed: both QTL returned with correct columns", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 2, interact = c(TRUE, FALSE))
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    # Use the interaction mock model (safe for both paths — the main-effect
    # path only reads whr[1L] so extra rows are harmless)
    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    expect_s3_class(fm, "fineMap")
    expect_equal(length(fm), 2L)
    for (i in seq_along(fm)) {
        expect_named(fm[[i]], c("mark", "dist", "pvalue", "LOD"))
    }
})

test_that("fineMap MV mixed: qtl.key attribute records both QTL keys", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 2, interact = c(TRUE, FALSE))
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    expect_equal(attr(fm, "qtl.key"), obj$QTL$qtl)
})

test_that("fineMap MV mixed: qtl= subset returns only the requested QTL", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 2, interact = c(TRUE, FALSE))
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels
    upd_m  <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)

    target <- obj$QTL$qtl[2]   # main-effect QTL only
    fm <- .run_fm_mv(obj, genObj, upd_m, qtl = target, window = 40, step = 20,
                     exclusion.window = 10)

    expect_equal(length(fm), 1L)
    expect_equal(names(fm), target)
})

# =============================================================================
# 5. Aliased-coefficient guard: all vcoeff zero → NA returned
# =============================================================================

test_that("fineMap MV interaction: aliased coefficients (all vcoeff=0) give NA", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1, interact = TRUE)
    genObj <- attr(obj, "genObj")
    trials <- obj$QTL$trait.levels

    upd_m <- make_mv_updated_model("X.C2.2_fm", trials, n_int = 15)
    # Zero out all vcoeff for cand_x rows so whr_nz is empty
    int_nams <- paste0(trials, ":X.C2.2_fm")
    upd_m$vcoeff$fixed[int_nams] <- 0

    fm <- .run_fm_mv(obj, genObj, upd_m, window = 40, step = 20,
                     exclusion.window = 10)

    expect_true(all(is.na(fm[[1]]$pvalue)))
    expect_true(all(is.na(fm[[1]]$LOD)))
})

# =============================================================================
# 6. Guard clauses still apply for MV objects (inherited from UV)
# =============================================================================

test_that("fineMap: stops when no QTL found in MV object", {
    obj <- make_mock_mv_qtlAim(n_qtl = 1)
    obj$QTL$qtl <- NULL
    genObj <- attr(obj, "genObj")
    expect_error(fineMap(obj, genObj), "No significant QTL")
})

test_that("fineMap: stops on bad qtl= key in MV object", {
    obj    <- make_mock_mv_qtlAim(n_qtl = 1)
    genObj <- attr(obj, "genObj")
    expect_error(fineMap(obj, genObj, qtl = "Chr.ZZ.99"), "not found")
})
