# =============================================================================
# test-qtlAim.R
# Tests for qtlAim() generic, qtlAim.default(), and the guard clauses /
# argument validation in qtlAim.asreml().
#
# Strategy:
#   - Generic dispatch and .default() need no ASReml at all.
#   - Guard clauses in qtlAim.asreml() fire in two groups:
#       Group A: inside .validateModel() — merge.by, method, selection, breakout
#       Group B: after .validateModel() — genObj class, gen.type, column matching
#     For Group B we mock .validateModel() so it returns successfully, letting
#     us test the qtlAim-specific guards without a real ASReml model.
#   - Trace-file behaviour needs only a converged mock model.
#   - The full forward-selection loop requires ASReml — tested only in
#     test-integration.R (skipped if ASReml is absent).
# =============================================================================

# ---- helper: build a mock model whose call$data eval() works ----------------
# .validateModel does: phenoData <- eval(baseModel$call$data)
# We embed the data frame directly in the call via bquote(local(.)) so that
# eval() always finds it regardless of frame.
make_inline_model <- function(ids = paste0("L", 1:20),
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
        loglik          = -50,
        sigma2          = 1,
        vparameters     = c(id = 0.5, "R!variance" = 1),
        vparameters.con = c(0L, 0L),
        coefficients    = list(fixed  = matrix(5, 1, 1,
                                               dimnames = list("(Intercept)", "effect")),
                               random = matrix(rnorm(length(ids)), length(ids), 1,
                                               dimnames = list(paste0(merge.by, "_", ids), "effect"))),
        formulae        = list(fixed = yld ~ 1, random = as.formula(paste("~", merge.by))),
        mf              = df,
        call            = bquote(asreml(fixed  = yld ~ 1,
                                        random = .(as.formula(paste("~", merge.by))),
                                        data   = .(df)))
    )
    class(m) <- c("asreml", "list")
    m
}

# ---- mock .validateModel return value ----------------------------------------
# Returns the minimum list that qtlAim.asreml expects from .validateModel(),
# so we can test the post-validate guard clauses in isolation.
make_validate_return <- function(model, phenoData) {
    list(baseModel = model,
         asremlEnv = list(),
         phenoData = phenoData)
}

# =============================================================================
# 1. Generic dispatch
# =============================================================================

test_that("qtlAim() dispatches via UseMethod", {
    expect_true(isS3stdGeneric(wgAim:::qtlAim))
})

test_that("qtlAim.default() stops with 'asreml' message", {
    expect_error(qtlAim(list()), regexp = "asreml")
    expect_error(qtlAim("not_a_model"), regexp = "asreml")
    expect_error(qtlAim(42), regexp = "asreml")
})

test_that("gwasAim.default() stops with 'asreml' message", {
    expect_error(gwasAim(list()), regexp = "asreml")
})

test_that("gpAim.default() stops with 'asreml' message", {
    expect_error(gpAim(list()), regexp = "asreml")
})

# =============================================================================
# 2. .validateModel guard clauses (tested via qtlAim.asreml)
# =============================================================================

test_that("qtlAim: merge.by=NULL stops before genObj check", {
    m <- make_inline_model()
    # merge.by NULL triggers .validateModel stop — no ASReml needed because
    # the stop fires before update() is called
    expect_error(
        qtlAim(m, genObj = make_wgCross_interval(), merge.by = NULL),
        regexp = "merge.by"
    )
})

test_that("qtlAim: bad method stops", {
    m <- make_inline_model()
    expect_error(
        qtlAim(m, genObj = make_wgCross_interval(),
                merge.by = "id", method = "bayes"),
        regexp = "method"
    )
})

test_that("qtlAim: bad selection stops", {
    m <- make_inline_model()
    expect_error(
        qtlAim(m, genObj = make_wgCross_interval(),
                merge.by = "id", selection = "random"),
        regexp = "selection"
    )
})

test_that("qtlAim: breakout=0 stops", {
    m <- make_inline_model()
    expect_error(
        qtlAim(m, genObj = make_wgCross_interval(),
                merge.by = "id", breakout = 0),
        regexp = "breakout"
    )
})

test_that("qtlAim: breakout=-2 stops", {
    m <- make_inline_model()
    expect_error(
        qtlAim(m, genObj = make_wgCross_interval(),
                merge.by = "id", breakout = -2),
        regexp = "breakout"
    )
})

# =============================================================================
# 3. qtlAim-specific guard clauses (Group B — after .validateModel)
#    We mock .validateModel to return successfully, then test the remaining
#    guards which are unique to qtlAim.asreml.
# =============================================================================

test_that("qtlAim: missing genObj stops", {
    m   <- make_inline_model()
    vd  <- make_validate_return(m, m$mf)
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                qtlAim(m, merge.by = "id"),
                regexp = "genObj"
            )
        }
    )
})

test_that("qtlAim: wrong genObj class (not wgCross) stops", {
    m      <- make_inline_model()
    vd     <- make_validate_return(m, m$mf)
    bad_go <- list()  # plain list, not wgCross
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                qtlAim(m, genObj = bad_go, merge.by = "id"),
                regexp = "wgCross"
            )
        }
    )
})

test_that("qtlAim: gen.type='interval' with marker-only wgCross stops", {
    m      <- make_inline_model()
    vd     <- make_validate_return(m, m$mf)
    gen_m  <- make_wgCross_marker()  # no interval.data
    # give it a matching id column
    gen_m$pheno$id <- paste0("L", seq_len(nrow(gen_m$pheno)))
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                qtlAim(m, genObj = gen_m, merge.by = "id",
                        gen.type = "interval"),
                regexp = "interval"
            )
        }
    )
})

test_that("qtlAim: gen.type='interval' with wgPanel stops", {
    m     <- make_inline_model()
    vd    <- make_validate_return(m, m$mf)
    panel <- make_wgPanel()
    panel$pheno$id <- paste0("L", seq_len(nrow(panel$pheno)))
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            # wgPanel inherits wgCross so passes the class check, then hits
            # the gen.type="interval" + wgPanel check... actually qtlAim
            # checks inherits(genObj, "wgCross") which wgPanel passes.
            # But gen.type="interval" with wgPanel: wgPanel$type = "marker"
            # so the stop("gen.type='interval' but genObj contains no interval") fires.
            expect_error(
                qtlAim(m, genObj = panel, merge.by = "id",
                        gen.type = "interval"),
                regexp = "interval"
            )
        }
    )
})

test_that("qtlAim: no overlapping lines between genObj and phenoData stops", {
    m      <- make_inline_model(ids = paste0("L", 1:20))
    vd     <- make_validate_return(m, m$mf)
    gen_i  <- make_wgCross_interval()
    # Rename pheno ids so nothing matches
    gen_i$pheno$id <- paste0("X", seq_len(nrow(gen_i$pheno)))
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                qtlAim(m, genObj = gen_i, merge.by = "id"),
                regexp = "do not match"
            )
        }
    )
})

test_that("qtlAim: merge.by column missing from genObj stops", {
    m     <- make_inline_model()
    vd    <- make_validate_return(m, m$mf)
    gen_i <- make_wgCross_interval()
    # Remove the id column from genObj$pheno
    gen_i$pheno <- gen_i$pheno[, names(gen_i$pheno) != "id", drop = FALSE]
    with_mocked_bindings(
        .validateModel = function(...) vd,
        .package = "wgAim",
        {
            expect_error(
                qtlAim(m, genObj = gen_i, merge.by = "id"),
                regexp = "column"
            )
        }
    )
})

# =============================================================================
# 4. Trace file creation
# =============================================================================

test_that("qtlAim: trace=character creates a file (sink cleaned up on exit)", {
    skip_if_not_installed("asreml")  # needs working model for phase 3+
    # But we can test the file creation path by checking that the file is
    # created immediately on entry, even if the analysis later fails.
    # Use tryCatch so a downstream error doesn't mask the assertion.
    tf <- tempfile(fileext = ".txt")
    on.exit(unlink(tf), add = TRUE)
    m <- make_inline_model()
    tryCatch(
        qtlAim(m, genObj = make_wgCross_interval(),
                merge.by = "id", trace = tf),
        error = function(e) NULL
    )
    # File must exist (created by sink) even if the analysis errored
    expect_true(file.exists(tf))
    # Sink must have been cleaned up (no active sinks remaining)
    expect_equal(sink.number(), 0L)
})

# =============================================================================
# 5. Argument pass-through defaults (visible via mocking)
# =============================================================================

test_that("qtlAim: gen.type defaults to genObj$type when not supplied", {
    m    <- make_inline_model()
    vd   <- make_validate_return(m, m$mf)
    gen_i <- make_wgCross_interval()  # $type = "interval"
    gen_i$pheno$id <- paste0("L", seq_len(nrow(gen_i$pheno)))
    m$mf$id <- factor(paste0("L", seq_len(nrow(m$mf))))

    # Capture the gen.type that reaches .buildGenoData
    captured_type <- NULL
    with_mocked_bindings(
        .validateModel  = function(...) vd,
        .buildGenoData  = function(intervalObj, gen.type, glines, plines) {
            captured_type <<- gen.type
            # Return minimal structure to avoid further errors
            nd <- ncol(intervalObj$geno[[1]]$imputed.data) * length(intervalObj$geno)
            list(genoData = matrix(0, nrow = 1, ncol = nd),
                 mnams    = paste0("Chr.C1.", seq_len(nd)),
                 state    = rep(1, nd))
        },
        .fixLines       = function(...) stop("stop here"),
        .package = "wgAim",
        {
            tryCatch(
                qtlAim(m, genObj = gen_i, merge.by = "id"),
                error = function(e) NULL
            )
        }
    )
    expect_equal(captured_type, "interval")
})

test_that("qtlAim: TypeI, exclusion.window, verboseLev defaults are numeric", {
    # Just confirm that qtlAim.asreml() signature defaults are as documented
    formals_list <- formals(wgAim:::qtlAim.asreml)
    expect_equal(formals_list$TypeI,            0.05)
    expect_equal(formals_list$exclusion.window, 20)
    expect_equal(formals_list$verboseLev,       0)
    expect_equal(eval(formals_list$breakout),   -1)    # stored as call `-1`
    expect_equal(as.character(formals_list$method),    "fixed")
    expect_equal(as.character(formals_list$selection), "interval")
})
