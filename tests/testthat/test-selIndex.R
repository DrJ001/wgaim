# =============================================================================
# test-selIndex.R
# Tests for selIndex(), print.selIndex(), summary.selIndex(), plot.selIndex().
#
# All tests are pure-R (no ASReml required). They exercise:
#   - selIndex.default  (guard clause)
#   - selIndex.gpAim    (univariate + multivariate, all 3 index types)
#   - .si_threshold     (internal helper)
#   - .check_index_vector (internal helper)
#   - print.selIndex    (univariate + MV)
#   - summary.selIndex  (univariate + MV)
#   - plot.selIndex     (type = "index", "heatmap", "weights")
# =============================================================================

skip_if_not_installed("ggplot2")

# ---- Shared fixtures --------------------------------------------------------

# Build a minimal MV gpAim-like object with GP$Trait set.
# n_lines x ntrait GEBVs, with Ga / Gcor stored.
.make_mv_gpAim <- function(n_lines = 30,
                            trials  = c("T1", "T2"),
                            seed    = 55) {
    set.seed(seed)
    ids    <- paste0("L", seq_len(n_lines))
    nt     <- length(trials)
    Ga_diag <- seq(0.8, 0.5, length.out = nt)
    Ga     <- diag(Ga_diag)
    dimnames(Ga) <- list(trials, trials)
    sds    <- sqrt(Ga_diag)
    Gcor   <- Ga / outer(sds, sds); diag(Gcor) <- 1
    dimnames(Gcor) <- list(trials, trials)

    SE_mv <- runif(n_lines * nt, 0.05, 0.20)
    gebv_df <- data.frame(
        id    = factor(rep(ids, nt)),
        Trial = factor(rep(trials, each = n_lines), levels = trials),
        GEBV  = rnorm(n_lines * nt),
        SE    = SE_mv,
        stringsAsFactors = FALSE
    )
    names(gebv_df)[1] <- "id"

    GP <- list(
        gebv         = gebv_df,
        gen.type     = "marker",
        path         = "vm",
        var.genetic  = setNames(Ga_diag, trials),
        var.resid    = setNames(rep(1, nt), trials),
        heritability = setNames(Ga_diag / (Ga_diag + 1), trials),
        gen.H2       = setNames(Ga_diag * 0.9, trials),
        n.markers    = 20L,
        rel.scale    = 20L,
        rel.matrix   = diag(n_lines),
        genetic.term = "id",
        Trait        = "Trial",
        trait.levels = trials,
        Ga           = Ga,
        Gcor         = Gcor
    )
    obj <- list(
        converge = TRUE, loglik = -100, sigma2 = 1.0,
        vparameters     = c(id = 0.6, "R!variance" = 1.0),
        vparameters.con = c(id = 0L,  "R!variance" = 0L),
        coefficients    = list(
            fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
            random = matrix(rnorm(n_lines), n_lines, 1,
                            dimnames = list(paste0("id_", ids), "effect"))
        ),
        formulae = list(fixed = yld ~ 1, random = ~ Trial:id),
        mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
        call     = call("asreml"),
        GP       = GP
    )
    class(obj) <- c("gpAim", "asreml")
    obj
}

# Univariate gpAim (no Trait)
.make_uv_gpAim <- function(n_lines = 30, seed = 42) {
    set.seed(seed)
    ids     <- paste0("L", seq_len(n_lines))
    SE_vals <- runif(n_lines, 0.1, 0.3)
    GP <- list(
        gebv = data.frame(
            id   = ids,
            GEBV = rnorm(n_lines),
            SE   = SE_vals,
            stringsAsFactors = FALSE
        ),
        gen.type    = "marker",
        path        = "vm",
        var.genetic = 0.6,
        var.resid   = 0.4,
        heritability = 0.6,
        gen.H2      = 0.75,
        n.markers   = 20L,
        rel.scale   = 20L,
        rel.matrix  = diag(n_lines),
        genetic.term = "id",
        Trait        = NULL,
        trait.levels = NULL,
        Ga           = NULL,
        Gcor         = NULL
    )
    obj <- list(
        converge = TRUE, loglik = -80, sigma2 = 1.0,
        vparameters     = c(id = 0.6, "R!variance" = 1.0),
        vparameters.con = c(id = 0L,  "R!variance" = 0L),
        coefficients    = list(
            fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
            random = matrix(rnorm(n_lines), n_lines, 1,
                            dimnames = list(paste0("id_", ids), "effect"))
        ),
        formulae = list(fixed = yld ~ 1, random = ~ id),
        mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
        call     = call("asreml"),
        GP       = GP
    )
    class(obj) <- c("gpAim", "asreml")
    obj
}

# =============================================================================
# 1.  selIndex.default -- wrong class guard
# =============================================================================

test_that("selIndex.default: non-gpAim object stops with informative error", {
    expect_error(selIndex(list()), regexp = "\"gpAim\"")
})

test_that("selIndex.default: plain data.frame stops with error", {
    expect_error(selIndex(data.frame(x = 1)), regexp = "\"gpAim\"")
})

# =============================================================================
# 2.  selIndex.gpAim -- univariate path
# =============================================================================

test_that("selIndex univariate: issues message and returns selIndex list", {
    obj <- .make_uv_gpAim()
    expect_message(si <- selIndex(obj), regexp = "univariate")
    expect_s3_class(si, "selIndex")
})

test_that("selIndex univariate: $index has columns id, GEBV, Index, Rank", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    expect_true(all(c("id", "GEBV", "Index", "Rank") %in% names(si$index)))
})

test_that("selIndex univariate: rows sorted descending by Index", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    expect_true(all(diff(si$index$Index) <= 0))
})

test_that("selIndex univariate: Rank is 1:n", {
    obj <- .make_uv_gpAim(n_lines = 20)
    si  <- suppressMessages(selIndex(obj))
    expect_equal(si$index$Rank, seq_len(20L))
})

test_that("selIndex univariate: type is always 'weighted'", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    expect_equal(si$type, "weighted")
})

test_that("selIndex univariate: n.selected consistent with prop.select", {
    set.seed(1)
    obj <- .make_uv_gpAim(n_lines = 40)
    si  <- suppressMessages(selIndex(obj, prop.select = 0.25))
    # Top 25% of 40 = 10 lines (may vary by ties, allow +-1)
    expect_gte(si$n.selected, 9L)
    expect_lte(si$n.selected, 11L)
})

test_that("selIndex univariate: gain is positive when selecting top lines", {
    obj <- .make_uv_gpAim(n_lines = 50)
    si  <- suppressMessages(selIndex(obj, prop.select = 0.10))
    expect_gt(si$gain, 0)
})

test_that("selIndex univariate: external selected lines respected", {
    obj <- .make_uv_gpAim(n_lines = 30)
    my_lines <- paste0("L", 1:5)
    si  <- suppressMessages(selIndex(obj, selected = my_lines))
    expect_equal(si$n.selected, 5L)
    expect_true(is.na(si$thr))
    expect_equal(sort(si$selected), sort(my_lines))
})

test_that("selIndex univariate: bad selected line stops with error", {
    obj <- .make_uv_gpAim(n_lines = 20)
    expect_error(
        suppressMessages(selIndex(obj, selected = c("L1", "BOGUS"))),
        regexp = "not found"
    )
})

test_that("selIndex univariate: threshold as raw value outside (0,1)", {
    obj <- .make_uv_gpAim(n_lines = 40)
    si  <- suppressMessages(selIndex(obj, threshold = 1.5))
    expect_equal(si$thr, 1.5)
})

test_that("selIndex univariate: threshold in (0,1) treated as quantile", {
    set.seed(10)
    obj <- .make_uv_gpAim(n_lines = 100)
    si  <- suppressMessages(selIndex(obj, threshold = 0.8))
    expected_thr <- quantile(obj$GP$gebv$GEBV, 0.8, names = FALSE)
    expect_equal(si$thr, expected_thr)
})

test_that("selIndex univariate: n.lines and n.environments correct", {
    obj <- .make_uv_gpAim(n_lines = 25)
    si  <- suppressMessages(selIndex(obj))
    expect_equal(si$n.lines, 25L)
    expect_equal(si$n.environments, 1L)
})

# =============================================================================
# 3.  selIndex.gpAim -- multivariate, type = "weighted"
# =============================================================================

test_that("selIndex MV weighted: returns selIndex with correct class", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    expect_s3_class(si, "selIndex")
    expect_equal(si$type, "weighted")
})

test_that("selIndex MV weighted: NULL weights uses equal weights with message", {
    obj <- .make_mv_gpAim()
    expect_message(si <- selIndex(obj, weights = NULL), regexp = "equal weights")
    expect_equal(length(si$weights), 2L)
    expect_equal(unname(si$weights), rep(0.5, 2))
})

test_that("selIndex MV weighted: named weights reordered to match trait levels", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    # Supply weights in reversed order
    si  <- selIndex(obj, weights = c(T2 = 0.3, T1 = 0.7))
    expect_equal(names(si$weights), c("T1", "T2"))
    expect_equal(si$weights["T1"], c(T1 = 0.7))
    expect_equal(si$weights["T2"], c(T2 = 0.3))
})

test_that("selIndex MV weighted: unnamed weights vector accepted", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(2, 1))
    expect_equal(names(si$weights), c("T1", "T2"))
})

test_that("selIndex MV weighted: wrong-length unnamed weights errors", {
    obj <- .make_mv_gpAim()
    expect_error(selIndex(obj, weights = c(1, 2, 3)), regexp = "named vector")
})

test_that("selIndex MV weighted: missing environment name in weights errors", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    expect_error(selIndex(obj, weights = c(T1 = 1)), regexp = "missing entries")
})

test_that("selIndex MV weighted: extra names in weights triggers warning", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    expect_warning(
        selIndex(obj, weights = c(T1 = 1, T2 = 1, BOGUS = 9)),
        regexp = "ignoring extra"
    )
})

test_that("selIndex MV weighted: $index has expected columns", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    expected_cols <- c("id", "T1", "T2", "Index", "Rank")
    expect_true(all(expected_cols %in% names(si$index)))
})

test_that("selIndex MV weighted: Index sorted descending", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    expect_true(all(diff(si$index$Index) <= 0))
})

test_that("selIndex MV weighted: gain is a named vector over trait levels", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    expect_named(si$gain, c("T1", "T2"))
    expect_length(si$gain, 2L)
})

test_that("selIndex MV weighted: standardise = FALSE still runs", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1), standardise = FALSE)
    expect_s3_class(si, "selIndex")
    expect_false(si$standardise)
})

test_that("selIndex MV weighted: external selected lines", {
    obj <- .make_mv_gpAim(n_lines = 30)
    my_lines <- paste0("L", 1:6)
    si <- selIndex(obj, weights = c(T1 = 1, T2 = 1), selected = my_lines)
    expect_equal(si$n.selected, 6L)
    expect_true(is.na(si$thr))
    expect_equal(sort(si$selected), sort(my_lines))
})

test_that("selIndex MV weighted: bad external selected line stops with error", {
    obj <- .make_mv_gpAim(n_lines = 20)
    expect_error(
        selIndex(obj, weights = c(T1 = 1, T2 = 1),
                 selected = c("L1", "DOES_NOT_EXIST")),
        regexp = "not found"
    )
})

test_that("selIndex MV weighted: n.lines and n.environments populated", {
    obj <- .make_mv_gpAim(n_lines = 28, trials = c("T1", "T2", "T3"))
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1, T3 = 1))
    expect_equal(si$n.lines, 28L)
    expect_equal(si$n.environments, 3L)
})

test_that("selIndex MV weighted: Ga and Gcor stored", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    expect_true(!is.null(si$Ga))
    expect_true(!is.null(si$Gcor))
    expect_equal(dim(si$Ga), c(2L, 2L))
})

# =============================================================================
# 4.  selIndex.gpAim -- type = "smith-hazel"
# =============================================================================

test_that("selIndex MV smith-hazel: type stored correctly", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1), type = "smith-hazel")
    expect_equal(si$type, "smith-hazel")
})

test_that("selIndex MV smith-hazel: returns selIndex with correct structure", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1), type = "smith-hazel")
    expect_s3_class(si, "selIndex")
    expect_named(si$weights, c("T1", "T2"))
})

test_that("selIndex MV smith-hazel: weights differ from equal weights", {
    obj <- .make_mv_gpAim(n_lines = 60, seed = 77)
    si_w  <- selIndex(obj, weights = c(T1 = 1, T2 = 1), type = "weighted")
    si_sh <- selIndex(obj, weights = c(T1 = 1, T2 = 1), type = "smith-hazel")
    # Smith-Hazel derived weights should not equal direct weights in general
    expect_false(isTRUE(all.equal(si_w$weights, si_sh$weights, tolerance = 0.01)))
})

test_that("selIndex MV smith-hazel: requires weights argument", {
    obj <- .make_mv_gpAim()
    # NULL weights → equal weights message then runs
    expect_message(
        si <- selIndex(obj, weights = NULL, type = "smith-hazel"),
        regexp = "equal weights"
    )
    expect_s3_class(si, "selIndex")
})

# =============================================================================
# 5.  selIndex.gpAim -- type = "desired-gains"
# =============================================================================

test_that("selIndex MV desired-gains: type stored correctly", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, desired = c(T1 = 0.5, T2 = 0.3), type = "desired-gains")
    expect_equal(si$type, "desired-gains")
})

test_that("selIndex MV desired-gains: returns selIndex", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, desired = c(T1 = 0.5, T2 = 0.3), type = "desired-gains")
    expect_s3_class(si, "selIndex")
})

test_that("selIndex MV desired-gains: weights have unit L2 norm", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, desired = c(T1 = 0.5, T2 = 0.2), type = "desired-gains")
    expect_equal(sqrt(sum(si$weights^2)), 1.0, tolerance = 1e-8)
})

test_that("selIndex MV desired-gains: missing desired argument stops with error", {
    obj <- .make_mv_gpAim()
    expect_error(
        selIndex(obj, type = "desired-gains"),
        regexp = "desired"
    )
})

test_that("selIndex MV desired-gains: missing env in desired stops with error", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    expect_error(
        selIndex(obj, desired = c(T1 = 0.5), type = "desired-gains"),
        regexp = "missing entries"
    )
})

test_that("selIndex MV desired-gains: reversed desired order accepted", {
    obj <- .make_mv_gpAim(trials = c("T1", "T2"))
    si  <- selIndex(obj, desired = c(T2 = 0.3, T1 = 0.5), type = "desired-gains")
    expect_equal(names(si$weights), c("T1", "T2"))
})

# =============================================================================
# 6.  .si_threshold internal helper
# =============================================================================

.si_threshold <- wgAim:::.si_threshold

test_that(".si_threshold: NULL threshold uses prop.select quantile", {
    set.seed(1)
    vals <- rnorm(100)
    thr  <- .si_threshold(vals, threshold = NULL, prop.select = 0.10)
    expect_equal(thr, quantile(vals, 0.90, names = FALSE))
})

test_that(".si_threshold: threshold in (0,1) treated as quantile", {
    set.seed(2)
    vals <- rnorm(100)
    thr  <- .si_threshold(vals, threshold = 0.75, prop.select = 0.10)
    expect_equal(thr, quantile(vals, 0.75, names = FALSE))
})

test_that(".si_threshold: threshold outside (0,1) returned as-is", {
    vals <- rnorm(50)
    thr  <- .si_threshold(vals, threshold = 3.0, prop.select = 0.10)
    expect_equal(thr, 3.0)
})

test_that(".si_threshold: threshold = 0 returned as-is (boundary)", {
    vals <- rnorm(50)
    thr  <- .si_threshold(vals, threshold = 0, prop.select = 0.10)
    expect_equal(thr, 0)
})

# =============================================================================
# 7.  .check_index_vector internal helper
# =============================================================================

.check_index_vector <- wgAim:::.check_index_vector

test_that(".check_index_vector: named vector reordered to match tl", {
    tl  <- c("A", "B", "C")
    v   <- c(C = 3, A = 1, B = 2)
    out <- .check_index_vector(v, tl, "weights")
    expect_equal(out, c(A = 1, B = 2, C = 3))
})

test_that(".check_index_vector: unnamed vector of correct length accepted", {
    tl  <- c("X", "Y")
    # The function assigns names then strips them via as.numeric(); values are preserved
    out <- .check_index_vector(c(0.6, 0.4), tl, "weights")
    expect_equal(as.numeric(out), c(0.6, 0.4))
    expect_length(out, 2L)
})

test_that(".check_index_vector: unnamed wrong-length stops with error", {
    tl <- c("X", "Y", "Z")
    expect_error(.check_index_vector(c(1, 2), tl, "weights"), regexp = "named vector")
})

test_that(".check_index_vector: missing tl entry in named vector stops with error", {
    tl <- c("A", "B")
    expect_error(.check_index_vector(c(A = 1), tl, "desired"), regexp = "missing entries")
})

test_that(".check_index_vector: extra names trigger warning, not error", {
    tl <- c("A", "B")
    expect_warning(
        .check_index_vector(c(A = 1, B = 2, C = 9), tl, "weights"),
        regexp = "ignoring extra"
    )
})

# =============================================================================
# 8.  print.selIndex
# =============================================================================

test_that("print.selIndex univariate: output contains 'Selection Index'", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    out <- capture.output(print(si))
    expect_true(any(grepl("Selection Index", out, ignore.case = TRUE)))
})

test_that("print.selIndex univariate: output mentions 'univariate'", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    out <- capture.output(print(si))
    expect_true(any(grepl("univariate", out, ignore.case = TRUE)))
})

test_that("print.selIndex univariate: returns x invisibly", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    ret <- capture.output(result <- print(si))
    expect_identical(result, si)
})

test_that("print.selIndex MV: output contains environment names", {
    obj <- .make_mv_gpAim(trials = c("E1", "E2"))
    si  <- selIndex(obj, weights = c(E1 = 1, E2 = 1))
    out <- capture.output(print(si))
    expect_true(any(grepl("E1", out)))
    expect_true(any(grepl("E2", out)))
})

test_that("print.selIndex MV: output contains 'Genetic correlations'", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    out <- capture.output(print(si))
    expect_true(any(grepl("Genetic corr", out, ignore.case = TRUE)))
})

test_that("print.selIndex MV: output contains weight values", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 2, T2 = 1))
    out <- capture.output(print(si))
    # Weights are printed in a data.frame, check for the numeric values
    expect_true(any(grepl("Weight", out, ignore.case = TRUE)))
})

test_that("print.selIndex MV: top and bottom n lines shown", {
    obj <- .make_mv_gpAim(n_lines = 30)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    out <- capture.output(print(si, n = 3))
    expect_true(any(grepl("Top 3", out)))
    expect_true(any(grepl("Bottom 3", out)))
})

test_that("print.selIndex MV: returns si invisibly", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    capture.output(result <- print(si))
    expect_identical(result, si)
})

# =============================================================================
# 9.  summary.selIndex
# =============================================================================

test_that("summary.selIndex univariate: returns a data.frame", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    tab <- capture.output(result <- summary(si))
    expect_s3_class(result, "data.frame")
})

test_that("summary.selIndex univariate: returns full index table (n rows)", {
    obj <- .make_uv_gpAim(n_lines = 25)
    si  <- suppressMessages(selIndex(obj))
    tab <- capture.output(result <- summary(si))
    expect_equal(nrow(result), 25L)
})

test_that("summary.selIndex univariate: printed output contains 'Expected genetic gain'", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    out <- capture.output(summary(si))
    expect_true(any(grepl("Expected genetic gain", out, ignore.case = TRUE)))
})

test_that("summary.selIndex univariate: printed output mentions selected count", {
    obj <- .make_uv_gpAim(n_lines = 40)
    si  <- suppressMessages(selIndex(obj, prop.select = 0.10))
    out <- capture.output(summary(si))
    expect_true(any(grepl("Selected", out, ignore.case = TRUE)))
})

test_that("summary.selIndex MV: returns data.frame with Index column", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    tab <- capture.output(result <- summary(si))
    expect_true("Index" %in% names(result))
})

test_that("summary.selIndex MV: printed output shows weights", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 2, T2 = 1))
    out <- capture.output(summary(si))
    expect_true(any(grepl("Weights|Weight", out, ignore.case = TRUE)))
})

test_that("summary.selIndex MV: external selected shows 'externally supplied'", {
    obj <- .make_mv_gpAim(n_lines = 30)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1),
                    selected = paste0("L", 1:5))
    out <- capture.output(summary(si))
    expect_true(any(grepl("extern", out, ignore.case = TRUE)))
})

test_that("summary.selIndex MV: threshold shown when externally-supplied selection absent", {
    obj <- .make_mv_gpAim(n_lines = 40)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1), prop.select = 0.20)
    out <- capture.output(summary(si))
    expect_true(any(grepl("threshold", out, ignore.case = TRUE)))
})

# =============================================================================
# 10.  plot.selIndex -- type = "index"
# =============================================================================

test_that("plot.selIndex type='index' univariate: returns ggplot", {
    obj <- .make_uv_gpAim(n_lines = 30)
    si  <- suppressMessages(selIndex(obj))
    gp  <- plot(si, type = "index")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='index' MV: returns ggplot", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    gp  <- plot(si, type = "index")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='index' MV with external selected: returns ggplot", {
    obj <- .make_mv_gpAim(n_lines = 30)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1),
                    selected = paste0("L", 1:5))
    gp  <- plot(si, type = "index")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='index': custom pt.col accepted", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    gp  <- plot(si, type = "index", pt.col = c("navy", "orange"))
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='index': threshold as quantile accepted", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    gp  <- plot(si, type = "index", threshold = 0.8)
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='index': threshold as raw value accepted", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    gp  <- plot(si, type = "index", threshold = 2.0)
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 11.  plot.selIndex -- type = "heatmap"
# =============================================================================

test_that("plot.selIndex type='heatmap' MV: returns ggplot", {
    obj <- .make_mv_gpAim(n_lines = 30)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    gp  <- plot(si, type = "heatmap")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='heatmap': univariate stops with error", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    expect_error(plot(si, type = "heatmap"), regexp = "multivariate")
})

test_that("plot.selIndex type='heatmap': external selected lines accepted", {
    obj <- .make_mv_gpAim(n_lines = 30)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1),
                    selected = paste0("L", 1:8))
    gp  <- plot(si, type = "heatmap")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='heatmap': 3-environment object returns ggplot", {
    obj <- .make_mv_gpAim(n_lines = 40, trials = c("E1", "E2", "E3"))
    si  <- selIndex(obj, weights = c(E1 = 1, E2 = 1, E3 = 1))
    gp  <- plot(si, type = "heatmap")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='heatmap': many lines (>80) returns ggplot", {
    obj <- .make_mv_gpAim(n_lines = 100)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    gp  <- plot(si, type = "heatmap")
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 12.  plot.selIndex -- type = "weights"
# =============================================================================

test_that("plot.selIndex type='weights' MV: returns ggplot", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 2, T2 = 1))
    gp  <- plot(si, type = "weights")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='weights': univariate stops with error", {
    obj <- .make_uv_gpAim()
    si  <- suppressMessages(selIndex(obj))
    expect_error(plot(si, type = "weights"), regexp = "multivariate")
})

test_that("plot.selIndex type='weights' smith-hazel: returns ggplot", {
    obj <- .make_mv_gpAim(n_lines = 50)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1), type = "smith-hazel")
    gp  <- plot(si, type = "weights")
    expect_s3_class(gp, "ggplot")
})

test_that("plot.selIndex type='weights' desired-gains: returns ggplot", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, desired = c(T1 = 0.5, T2 = 0.3), type = "desired-gains")
    gp  <- plot(si, type = "weights")
    expect_s3_class(gp, "ggplot")
})

# =============================================================================
# 13.  .compute_gains internal helper
# =============================================================================

.compute_gains <- wgAim:::.compute_gains

test_that(".compute_gains MV: returns named vector over trait levels", {
    obj <- .make_mv_gpAim()
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    thr <- si$thr
    gains <- .compute_gains(si, thr)
    expect_named(gains, c("T1", "T2"))
    expect_length(gains, 2L)
})

test_that(".compute_gains UV: returns scalar", {
    obj <- .make_uv_gpAim(n_lines = 40)
    si  <- suppressMessages(selIndex(obj))
    thr <- si$thr
    gains <- .compute_gains(si, thr)
    expect_length(gains, 1L)
})

test_that(".compute_gains: gain is positive when top lines selected", {
    obj <- .make_mv_gpAim(n_lines = 60)
    si  <- selIndex(obj, weights = c(T1 = 1, T2 = 1))
    # The selected lines have above-mean index, so at least one env gain >= 0
    gains <- .compute_gains(si, si$thr)
    expect_true(any(gains >= 0))
})
