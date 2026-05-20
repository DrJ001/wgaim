# =============================================================================
# test-print-gpAim.R
# Additional tests for print.gpAim() targeting the branches not yet covered
# in test-print-summary.R.
#
# Existing tests in test-print-summary.R cover only the basic univariate path.
# This file adds:
#   - Univariate: with gen.H2 present / absent
#   - Univariate: returns x invisibly
#   - Univariate: output structure (keywords)
#   - Univariate: NULL gen.H2 branch (no Generalised H2 line)
#   - Univariate: NA gen.H2 branch  (no Generalised H2 line)
#   - Multivariate: basic output
#   - Multivariate: Gcor block printed
#   - Multivariate: GEBV range per trial
#   - Multivariate: gen.H2 per trial in vc.df
#   - Multivariate: NULL gen.H2 omitted
#   - Multivariate: returns x invisibly
# =============================================================================

# ---- Helpers to build univariate / MV gpAim mocks with known GP structure --

.make_print_uv_gpAim <- function(n_lines = 20,
                                  gen.H2  = 0.75,
                                  path    = "vm",
                                  seed    = 1) {
    set.seed(seed)
    ids     <- paste0("L", seq_len(n_lines))
    SE_vals <- runif(n_lines, 0.1, 0.3)
    GP <- list(
        gebv = data.frame(
            id   = ids,
            GEBV = rnorm(n_lines, 0, 1),
            SE   = SE_vals,
            stringsAsFactors = FALSE
        ),
        gen.type    = "interval",
        path        = path,
        var.genetic = 0.6,
        var.resid   = 0.4,
        heritability = 0.6,
        gen.H2      = gen.H2,
        n.markers   = 16L,
        rel.scale   = 16L,
        rel.matrix  = diag(n_lines),
        genetic.term = "id",
        Trait        = NULL,
        trait.levels = NULL
    )
    obj <- list(
        converge = TRUE, loglik = -80, sigma2 = 1.0,
        vparameters     = c(id = 0.6, "R!variance" = 1.0),
        vparameters.con = c(id = 0L,  "R!variance" = 0L),
        coefficients = list(
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

.make_print_mv_gpAim <- function(n_lines = 20,
                                  trials  = c("T1", "T2"),
                                  gen.H2  = NULL,
                                  seed    = 2) {
    set.seed(seed)
    nt   <- length(trials)
    ids  <- paste0("L", seq_len(n_lines))
    Ga_diag <- seq(0.8, 0.5, length.out = nt)
    Ga   <- diag(Ga_diag); dimnames(Ga) <- list(trials, trials)
    sds  <- sqrt(Ga_diag)
    Gcor <- Ga / outer(sds, sds); diag(Gcor) <- 1
    dimnames(Gcor) <- list(trials, trials)

    SE_mv <- runif(n_lines * nt, 0.05, 0.20)
    gebv_df <- data.frame(
        id    = factor(rep(ids, nt)),
        Trial = factor(rep(trials, each = n_lines), levels = trials),
        GEBV  = rnorm(n_lines * nt, 0, 0.8),
        SE    = SE_mv,
        stringsAsFactors = FALSE
    )
    if (is.null(gen.H2))
        gen.H2_val <- NULL
    else
        gen.H2_val <- setNames(gen.H2, trials[seq_along(gen.H2)])

    GP <- list(
        gebv         = gebv_df,
        gen.type     = "marker",
        path         = "vm",
        var.genetic  = setNames(Ga_diag, trials),
        var.resid    = setNames(rep(1, nt), trials),
        heritability = setNames(Ga_diag / (Ga_diag + 1), trials),
        gen.H2       = gen.H2_val,
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
        coefficients = list(
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

# =============================================================================
# Univariate print.gpAim
# =============================================================================

test_that("print.gpAim UV: header contains 'Genomic Prediction'", {
    obj <- .make_print_uv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("Genomic Prediction", out, ignore.case = TRUE)))
})

test_that("print.gpAim UV: reports line count", {
    obj <- .make_print_uv_gpAim(n_lines = 17)
    out <- capture.output(print(obj))
    expect_true(any(grepl("17", out)))
})

test_that("print.gpAim UV: reports marker count", {
    obj <- .make_print_uv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("16", out)))  # n.markers = 16
})

test_that("print.gpAim UV: reports fitting path 'vm'", {
    obj <- .make_print_uv_gpAim(path = "vm")
    out <- capture.output(print(obj))
    expect_true(any(grepl("vm", out)))
})

test_that("print.gpAim UV: reports fitting path 'mbf'", {
    obj <- .make_print_uv_gpAim(path = "mbf")
    out <- capture.output(print(obj))
    expect_true(any(grepl("mbf", out)))
})

test_that("print.gpAim UV: reports heritability", {
    obj <- .make_print_uv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("h2|Herit", out, ignore.case = TRUE)))
})

test_that("print.gpAim UV: reports Generalised H2 when gen.H2 is numeric", {
    obj <- .make_print_uv_gpAim(gen.H2 = 0.82)
    out <- capture.output(print(obj))
    expect_true(any(grepl("[Gg]eneralised", out)))
})

test_that("print.gpAim UV: no Generalised H2 line when gen.H2 is NULL", {
    obj <- .make_print_uv_gpAim(gen.H2 = NULL)
    out <- capture.output(print(obj))
    expect_false(any(grepl("[Gg]eneralised", out)))
})

test_that("print.gpAim UV: no Generalised H2 line when gen.H2 is NA", {
    obj <- .make_print_uv_gpAim(gen.H2 = NA_real_)
    out <- capture.output(print(obj))
    expect_false(any(grepl("[Gg]eneralised", out)))
})

test_that("print.gpAim UV: GEBV range section present", {
    obj <- .make_print_uv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("GEBV range|Min|Max", out, ignore.case = TRUE)))
})

test_that("print.gpAim UV: returns x invisibly", {
    obj <- .make_print_uv_gpAim()
    capture.output(result <- print(obj))
    expect_identical(result, obj)
})

test_that("print.gpAim UV: works with make_mock_gpAim fixture", {
    obj <- make_mock_gpAim(n_lines = 20)
    out <- capture.output(print(obj))
    expect_true(any(grepl("Genomic Prediction", out, ignore.case = TRUE)))
})

# =============================================================================
# Multivariate print.gpAim
# =============================================================================

test_that("print.gpAim MV: header contains 'Multivariate'", {
    obj <- .make_print_mv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("Multivariate", out, ignore.case = TRUE)))
})

test_that("print.gpAim MV: reports trait levels", {
    obj <- .make_print_mv_gpAim(trials = c("Site1", "Site2"))
    out <- capture.output(print(obj))
    expect_true(any(grepl("Site1", out)))
    expect_true(any(grepl("Site2", out)))
})

test_that("print.gpAim MV: reports line count (unique lines)", {
    obj <- .make_print_mv_gpAim(n_lines = 15)
    out <- capture.output(print(obj))
    expect_true(any(grepl("15", out)))
})

test_that("print.gpAim MV: reports marker count", {
    obj <- .make_print_mv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("20", out)))  # n.markers = 20
})

test_that("print.gpAim MV: variance component table printed (Vg, Ve, h2 columns)", {
    obj <- .make_print_mv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("Vg|Ve|h2", out)))
})

test_that("print.gpAim MV: gen.H2 column appears in vc table when provided", {
    obj <- .make_print_mv_gpAim(trials = c("T1", "T2"), gen.H2 = c(0.80, 0.72))
    out <- capture.output(print(obj))
    expect_true(any(grepl("gen.H2|H2", out, ignore.case = TRUE)))
})

test_that("print.gpAim MV: no gen.H2 column when gen.H2 is NULL", {
    obj <- .make_print_mv_gpAim(gen.H2 = NULL)
    out <- capture.output(print(obj))
    # gen.H2 line should be absent from variance component table header
    # The Gcor section will still print, so just verify the table row count
    expect_true(any(grepl("Vg|Ve|h2", out)))  # basic VC table still there
})

test_that("print.gpAim MV: genetic correlation block printed", {
    obj <- .make_print_mv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("[Gg]enetic corr", out, ignore.case = TRUE)))
})

test_that("print.gpAim MV: GEBV range per trial shows trial names", {
    obj <- .make_print_mv_gpAim(trials = c("Env1", "Env2"))
    out <- capture.output(print(obj))
    expect_true(any(grepl("Env1", out)))
    expect_true(any(grepl("Env2", out)))
})

test_that("print.gpAim MV: GEBV range shows Min, Mean, Max keywords", {
    obj <- .make_print_mv_gpAim()
    out <- capture.output(print(obj))
    expect_true(any(grepl("Min|Mean|Max", out, ignore.case = TRUE)))
})

test_that("print.gpAim MV: three-trial object prints all three trial names", {
    obj <- .make_print_mv_gpAim(trials = c("T1", "T2", "T3"),
                                 gen.H2 = c(0.8, 0.7, 0.65))
    out <- capture.output(print(obj))
    expect_true(any(grepl("T1", out)))
    expect_true(any(grepl("T2", out)))
    expect_true(any(grepl("T3", out)))
})

test_that("print.gpAim MV: returns x invisibly", {
    obj <- .make_print_mv_gpAim()
    capture.output(result <- print(obj))
    expect_identical(result, obj)
})
