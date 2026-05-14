# =============================================================================
# test-gwasAim.R
# Tests for gwasAim() generic, gwasAim.default(), and the guard clauses /
# argument validation in gwasAim.asreml().
#
# Same strategy as test-qtlAim.R: mock .validateModel for Group B guards.
# gwasAim differs from qtlAim in:
#   - genObj must be "wgPanel" (not "wgCross")
#   - method / selection / gen.type are hard-coded ("fixed"/"interval"/"marker")
#   - $QTL$n.markers is set in results
# =============================================================================

# Reuse the same inline-model helper defined in test-qtlAim.R via helper-fixtures.R
# (testthat loads all helper-*.R before any test-*.R)

make_inline_model_gwas <- function(ids = paste0("S", 1:30),
                                    merge.by = "id",
                                    converge = TRUE) {
    df <- data.frame(
        id  = factor(ids),
        yld = rnorm(length(ids)),
        stringsAsFactors = FALSE
    )
    names(df)[1] <- merge.by
    m <- list(
        converge        = converge,
        loglik          = -60,
        sigma2          = 1,
        vparameters     = c(id = 0.4, "R!variance" = 1),
        vparameters.con = c(0L, 0L),
        coefficients    = list(fixed  = matrix(10, 1, 1,
                                               dimnames = list("(Intercept)", "effect")),
                               random = matrix(rnorm(length(ids)), length(ids), 1,
                                               dimnames = list(paste0(merge.by, "_", ids),
                                                               "effect"))),
        formulae        = list(fixed = yld ~ 1, random = as.formula(paste("~", merge.by))),
        mf              = df,
        call            = bquote(asreml(fixed  = yld ~ 1,
                                        random = .(as.formula(paste("~", merge.by))),
                                        data   = .(df)))
    )
    class(m) <- c("asreml", "list")
    m
}

make_validate_return_gwas <- function(model, phenoData) {
    list(baseModel = model,
         asremlEnv = list(),
         phenoData = phenoData)
}

# =============================================================================
# 1. Generic dispatch
# =============================================================================

test_that("gwasAim() is a function whose body calls UseMethod", {
    expect_true(is.function(wgAim:::gwasAim))
    expect_true(grepl("UseMethod", paste(deparse(body(wgAim:::gwasAim)), collapse = " ")))
})

test_that("gwasAim.default() stops with 'asreml' message", {
    expect_error(gwasAim(list()),        regexp = "asreml")
    expect_error(gwasAim("not_a_model"), regexp = "asreml")
    expect_error(gwasAim(42),            regexp = "asreml")
})

# =============================================================================
# 2. .validateModel guard clauses (shared with qtlAim)
# =============================================================================

test_that("gwasAim: merge.by=NULL stops in .validateModel", {
    m <- make_inline_model_gwas()
    expect_error(
        gwasAim(m, genObj = make_wgPanel(), merge.by = NULL),
        regexp = "merge.by"
    )
})

test_that("gwasAim: bad breakout stops in .validateModel", {
    m <- make_inline_model_gwas()
    expect_error(
        gwasAim(m, genObj = make_wgPanel(), merge.by = "id", breakout = 0),
        regexp = "breakout"
    )
})

# =============================================================================
# 3. gwasAim-specific guard clauses (Group B — after .validateModel)
# =============================================================================

test_that("gwasAim: missing genObj stops", {
    m  <- make_inline_model_gwas()
    vd <- make_validate_return_gwas(m, m$mf)
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(gwasAim(m, merge.by = "id"), regexp = "genObj")
        }
    )
})

test_that("gwasAim: wrong genObj class (wgCross, not wgPanel) stops", {
    m      <- make_inline_model_gwas()
    vd     <- make_validate_return_gwas(m, m$mf)
    # wgCross does NOT inherit wgPanel — this should trigger the class check
    bad_go <- make_wgCross_interval()
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                gwasAim(m, genObj = bad_go, merge.by = "id"),
                regexp = "wgPanel"
            )
        }
    )
})

test_that("gwasAim: plain list genObj stops", {
    m  <- make_inline_model_gwas()
    vd <- make_validate_return_gwas(m, m$mf)
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                gwasAim(m, genObj = list(), merge.by = "id"),
                regexp = "wgPanel"
            )
        }
    )
})

test_that("gwasAim: no overlapping lines between panel and phenoData stops", {
    m     <- make_inline_model_gwas(ids = paste0("S", 1:30))
    vd    <- make_validate_return_gwas(m, m$mf)
    panel <- make_wgPanel()
    # Rename panel pheno ids so nothing matches phenoData
    panel$pheno$id <- paste0("NOMATCH", seq_len(nrow(panel$pheno)))
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                gwasAim(m, genObj = panel, merge.by = "id"),
                regexp = "do not match"
            )
        }
    )
})

test_that("gwasAim: merge.by column missing from panel$pheno stops", {
    m     <- make_inline_model_gwas()
    vd    <- make_validate_return_gwas(m, m$mf)
    panel <- make_wgPanel()
    # Drop the id column from panel$pheno
    panel$pheno <- panel$pheno[, names(panel$pheno) != "id", drop = FALSE]
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                gwasAim(m, genObj = panel, merge.by = "id"),
                regexp = "column"
            )
        }
    )
})

# =============================================================================
# 4. Hard-coded engine constants
# =============================================================================

test_that("gwasAim: method/selection/type are hard-coded (not user args)", {
    # gwasAim.asreml does not expose method= or selection= in its signature
    formals_list <- formals(wgAim:::gwasAim.asreml)
    expect_false("method"    %in% names(formals_list))
    expect_false("selection" %in% names(formals_list))
    expect_false("gen.type"  %in% names(formals_list))
})

test_that("gwasAim: TypeI and exclusion.window defaults are correct", {
    formals_list <- formals(wgAim:::gwasAim.asreml)
    expect_equal(formals_list$TypeI,            0.05)
    expect_equal(formals_list$exclusion.window, 20)
    expect_equal(eval(formals_list$breakout),   -1)    # stored as call `-1`
    expect_equal(formals_list$verboseLev,       0)
})

# =============================================================================
# 5. n.markers captured correctly
# =============================================================================

test_that("gwasAim: n.markers captured and passed into .buildGenoData", {
    m     <- make_inline_model_gwas(ids = paste0("S", 1:30))
    vd    <- make_validate_return_gwas(m, m$mf)
    panel <- make_wgPanel(n_lines = 30, n_chr = 2, n_mar = 10)
    panel$pheno$id <- factor(paste0("S", seq_len(nrow(panel$pheno))))

    # The total marker count is n_chr * n_mar = 20
    # Capture it via mocking .buildGenoData
    captured_ncol <- NULL
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .buildGenoData = function(intervalObj, gen.type, glines, plines) {
            nd <- sum(sapply(intervalObj$geno, function(ch) ncol(ch$imputed.data)))
            captured_ncol <<- nd
            list(
                genoData = matrix(0, nrow = length(intersect(glines, plines)),
                                  ncol = nd,
                                  dimnames = list(intersect(as.character(glines),
                                                            as.character(plines)),
                                                  paste0("Chr.P1.", seq_len(nd)))),
                mnams    = paste0("Chr.P1.", seq_len(nd)),
                state    = rep(1, nd)
            )
        },
        .fixLines = function(...) stop("stop after capture"),
        .package  = "wgAim",
        {
            tryCatch(
                gwasAim(m, genObj = panel, merge.by = "id"),
                error = function(e) NULL
            )
        }
    )
    expect_equal(captured_ncol, 20L)  # 2 chr x 10 markers
})

# =============================================================================
# 6. Trace file creation
# =============================================================================

test_that("gwasAim: trace=character creates file and sink cleaned up", {
    skip_if_not_installed("asreml")
    tf <- tempfile(fileext = ".txt")
    on.exit(unlink(tf), add = TRUE)
    m <- make_inline_model_gwas()
    tryCatch(
        gwasAim(m, genObj = make_wgPanel(),
                merge.by = "id", trace = tf),
        error = function(e) NULL
    )
    expect_true(file.exists(tf))
    expect_equal(sink.number(), 0L)
})
