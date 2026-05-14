# =============================================================================
# test-gpAim.R
# Tests for gpAim() generic, gpAim.default(), and the guard clauses /
# argument validation in gpAim.asreml().
#
# gpAim does validation INLINE (not via .validateModel), so all guard clauses
# up to .buildGenoData can be reached with the mock model from fixtures alone,
# without any mocking of internal functions.
#
# The vm-vs-mbf path selection is also entirely testable because it is based
# purely on ncol(genoData) vs nrow(genoData), which .buildGenoData controls.
# We mock .buildGenoData to test this branching.
# =============================================================================

# ---- helper: converged inline model for gpAim (no .validateModel call) ------
# gpAim does its own inline validation, so the mock just needs:
#   $converge, $formulae (for asremlEnv), $call$data (for phenoData eval)

make_gp_model <- function(ids       = paste0("L", 1:40),
                            merge.by  = "id",
                            converge  = TRUE) {
    df <- data.frame(
        id  = factor(ids),
        yld = rnorm(length(ids)),
        stringsAsFactors = FALSE
    )
    names(df)[1] <- merge.by
    m <- list(
        converge        = converge,
        loglik          = -80,
        sigma2          = 1,
        vparameters     = c(id = 0.6, "R!variance" = 1),
        vparameters.con = c(0L, 0L),
        coefficients    = list(fixed  = matrix(5, 1, 1,
                                               dimnames = list("(Intercept)", "effect")),
                               random = matrix(rnorm(length(ids)), length(ids), 1,
                                               dimnames = list(paste0(merge.by, "_", ids),
                                                               "effect"))),
        formulae = list(fixed  = yld ~ 1,
                        random = as.formula(paste("~", merge.by))),
        mf   = df,
        call = bquote(asreml(fixed  = yld ~ 1,
                             random = .(as.formula(paste("~", merge.by))),
                             data   = .(df)))
    )
    class(m) <- c("asreml", "list")
    m
}

# =============================================================================
# 1. Generic dispatch
# =============================================================================

test_that("gpAim() is a function whose body calls UseMethod", {
    expect_true(is.function(wgAim:::gpAim))
    expect_true(grepl("UseMethod", paste(deparse(body(wgAim:::gpAim)), collapse = " ")))
})

test_that("gpAim.default() stops with 'asreml' message", {
    expect_error(gpAim(list()),        regexp = "asreml")
    expect_error(gpAim("not_a_model"), regexp = "asreml")
    expect_error(gpAim(42),            regexp = "asreml")
})

# =============================================================================
# 2. Inline guard clauses in gpAim.asreml (no mocking needed)
# =============================================================================

test_that("gpAim: merge.by=NULL stops immediately", {
    m <- make_gp_model()
    expect_error(
        gpAim(m, genObj = make_wgCross_interval(), merge.by = NULL),
        regexp = "merge.by"
    )
})

test_that("gpAim: bad gen.type stops", {
    m <- make_gp_model()
    expect_error(
        gpAim(m, genObj = make_wgCross_interval(),
               merge.by = "id", gen.type = "chromosome"),
        regexp = "gen.type"
    )
})

test_that("gpAim: gen.type='interval' with wgPanel stops", {
    m     <- make_gp_model()
    panel <- make_wgPanel()
    panel$pheno$id <- factor(paste0("L", seq_len(nrow(panel$pheno))))
    expect_error(
        gpAim(m, genObj = panel, merge.by = "id", gen.type = "interval"),
        regexp = "wgCross|interval"
    )
})

test_that("gpAim: gen.type='interval' with marker-only wgCross stops", {
    m     <- make_gp_model()
    gen_m <- make_wgCross_marker()
    gen_m$pheno$id <- factor(paste0("L", seq_len(nrow(gen_m$pheno))))
    expect_error(
        gpAim(m, genObj = gen_m, merge.by = "id", gen.type = "interval"),
        regexp = "interval"
    )
})

test_that("gpAim: wrong genObj class (plain list) stops", {
    m <- make_gp_model()
    expect_error(
        gpAim(m, genObj = list(), merge.by = "id"),
        regexp = "wgCross|wgPanel"
    )
})

test_that("gpAim: non-matching merge.by stops", {
    m     <- make_gp_model(ids = paste0("L", 1:40))
    gen_i <- make_wgCross_interval()
    # Override pheno ids so nothing matches the model's 'L1..L40' ids
    gen_i$pheno$id <- factor(paste0("NOMATCH", seq_len(nrow(gen_i$pheno))))
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id"),
        regexp = "do not match"
    )
})

test_that("gpAim: merge.by column missing from genObj$pheno stops", {
    m     <- make_gp_model()
    gen_i <- make_wgCross_interval()
    gen_i$pheno <- gen_i$pheno[, names(gen_i$pheno) != "id", drop = FALSE]
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id"),
        regexp = "column"
    )
})

test_that("gpAim: non-converged model stops if re-update also fails", {
    m            <- make_gp_model(converge = FALSE)
    m_bad        <- m
    # update() needs to return something with $converge=FALSE to trigger stop
    with_mocked_bindings(
        update = function(object, ...) { object$converge <- FALSE; object },
        .package = "wgAim",
        {
            expect_error(
                gpAim(m_bad, genObj = make_wgCross_interval(), merge.by = "id"),
                regexp = "not converged"
            )
        }
    )
})

# =============================================================================
# 3. vm vs mbf path selection
# =============================================================================

test_that("gpAim: vm path selected when markers > lines", {
    # markers > lines → use.vm = TRUE
    m     <- make_gp_model(ids = paste0("L", 1:10))   # 10 lines
    gen_i <- make_wgCross_interval(n_lines = 10, n_chr = 3, n_mar = 5)  # 15 markers
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    m$mf$id        <- factor(paste0("L", 1:10))
    m$call         <- bquote(asreml(fixed  = yld ~ 1,
                                    random = ~ id,
                                    data   = .(m$mf)))

    # Capture the path chosen by intercepting .constructCM (vm path calls this)
    vm_called <- FALSE
    with_mocked_bindings(
        .constructCM = function(...) {
            vm_called <<- TRUE
            # Return a minimal environment so assign works downstream
            e <- new.env()
            nd <- 10  # n_lines
            e$relm  <- diag(nd)
            dimnames(e$relm) <- list(paste0("L", 1:nd), paste0("L", 1:nd))
            e$scale <- 1
            mat <- matrix(rnorm(15 * nd), 15, nd)
            e$trans <- mat
            e
        },
        # Stop immediately after covObj assignment so we don't need update()
        update = function(object, ...) stop("stop after vm setup"),
        .package = "wgAim",
        {
            tryCatch(
                gpAim(m, genObj = gen_i, merge.by = "id"),
                error = function(e) NULL
            )
        }
    )
    expect_true(vm_called)
})

test_that("gpAim: mbf path selected when lines >= markers", {
    # lines >= markers → use.vm = FALSE
    m     <- make_gp_model(ids = paste0("L", 1:40))  # 40 lines
    gen_i <- make_wgCross_interval(n_lines = 40, n_chr = 2, n_mar = 5)  # 10 markers
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    m$mf$id        <- factor(paste0("L", 1:40))
    m$call         <- bquote(asreml(fixed  = yld ~ 1,
                                    random = ~ id,
                                    data   = .(m$mf)))

    # If mbf path: .constructCM is NOT called; instead assign(covObj) is a df
    vm_called <- FALSE
    with_mocked_bindings(
        .constructCM = function(...) { vm_called <<- TRUE; list() },
        update       = function(object, ...) stop("stop after mbf setup"),
        .package     = "wgAim",
        {
            tryCatch(
                gpAim(m, genObj = gen_i, merge.by = "id"),
                error = function(e) NULL
            )
        }
    )
    expect_false(vm_called)
})

test_that("gpAim: force=TRUE forces mbf path even when markers > lines", {
    m     <- make_gp_model(ids = paste0("L", 1:10))   # 10 lines
    gen_i <- make_wgCross_interval(n_lines = 10, n_chr = 3, n_mar = 5)  # 15 markers > 10 lines
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    m$mf$id        <- factor(paste0("L", 1:10))
    m$call         <- bquote(asreml(fixed  = yld ~ 1,
                                    random = ~ id,
                                    data   = .(m$mf)))

    vm_called <- FALSE
    with_mocked_bindings(
        .constructCM = function(...) { vm_called <<- TRUE; list() },
        update       = function(object, ...) stop("stop after mbf setup"),
        .package     = "wgAim",
        {
            tryCatch(
                gpAim(m, genObj = gen_i, merge.by = "id", force = TRUE),
                error = function(e) NULL
            )
        }
    )
    # force=TRUE → mbf path → .constructCM never called
    expect_false(vm_called)
})

# =============================================================================
# 4. Argument defaults
# =============================================================================

test_that("gpAim: default argument values match documentation", {
    formals_list <- formals(wgAim:::gpAim.asreml)
    expect_equal(as.character(formals_list$gen.type), "marker")
    expect_equal(formals_list$fix.lines, TRUE)
    expect_equal(formals_list$force,     FALSE)
    expect_equal(formals_list$trace,     TRUE)
})

# =============================================================================
# 5. Trace file creation
# =============================================================================

test_that("gpAim: trace=character creates file and sink cleaned up", {
    skip_if_not_installed("asreml")
    tf <- tempfile(fileext = ".txt")
    on.exit(unlink(tf), add = TRUE)
    m <- make_gp_model()
    tryCatch(
        gpAim(m, genObj = make_wgCross_interval(),
               merge.by = "id", trace = tf),
        error = function(e) NULL
    )
    expect_true(file.exists(tf))
    expect_equal(sink.number(), 0L)
})

# =============================================================================
# 6. gen.type="marker" works with both wgCross and wgPanel
# =============================================================================

test_that("gpAim: gen.type='marker' accepted for wgCross input (no stop)", {
    m     <- make_gp_model()
    gen_m <- make_wgCross_marker()
    gen_m$pheno$id <- factor(paste0("L", seq_len(nrow(gen_m$pheno))))
    m$mf$id        <- factor(paste0("L", seq_len(nrow(m$mf))))
    m$call         <- bquote(asreml(fixed  = yld ~ 1,
                                    random = ~ id,
                                    data   = .(m$mf)))

    # Should pass all guards; first ASReml-dependent call is update()
    with_mocked_bindings(
        update   = function(object, ...) stop("stop at update"),
        .package = "wgAim",
        {
            err <- tryCatch(
                gpAim(m, genObj = gen_m, merge.by = "id", gen.type = "marker"),
                error = function(e) conditionMessage(e)
            )
        }
    )
    # Error must be the mocked "stop at update", NOT a guard clause error
    expect_match(err, "stop at update")
})

test_that("gpAim: wgPanel accepted with gen.type='marker' (no stop)", {
    m     <- make_gp_model(ids = paste0("S", 1:50))
    panel <- make_wgPanel(n_lines = 50)
    panel$pheno$id <- factor(paste0("S", seq_len(nrow(panel$pheno))))
    m$mf$id        <- factor(paste0("S", 1:50))
    m$call         <- bquote(asreml(fixed  = yld ~ 1,
                                    random = ~ id,
                                    data   = .(m$mf)))

    with_mocked_bindings(
        update   = function(object, ...) stop("stop at update"),
        .package = "wgAim",
        {
            err <- tryCatch(
                gpAim(m, genObj = panel, merge.by = "id", gen.type = "marker"),
                error = function(e) conditionMessage(e)
            )
        }
    )
    expect_match(err, "stop at update")
})
