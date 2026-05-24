# =============================================================================
# test-engine-asreml.R
# Coverage tests for the ASReml-dependent engine functions.
#
# Every uncovered line in engine_select.R, engine_model.R, engine_effects.R,
# qtlAim.R, gwasAim.R, and fineMap.R is inside code that calls ASReml's
# update() or predict(). We cover those paths by mocking update() and
# predict() at the wgAim namespace level so no licence is needed.
#
# Mocking strategy:
#   update()  -> returns a minimal model-like list with all slots the
#                engine reads ($loglik, $vparameters, $vparameters.con,
#                $sigma2, $coefficients, $vcoeff, $call, $converge, $G.param)
#   predict() -> returns a minimal list with $pvals and $vcov
# =============================================================================

# ---- shared mock factories ---------------------------------------------------
# Note: make_updated_model() is defined in helper-fixtures.R so it is
# available to all test files without duplication.

# Build a minimal predict() return value for the vm path in .qtlSelect
make_predict_return <- function(merge.by = "id", ids = paste0("L", 1:10)) {
    pvals <- data.frame(
        id              = factor(ids),
        predicted.value = rnorm(length(ids), 0, 0.2),
        std.error       = runif(length(ids), 0.1, 0.3),
        stringsAsFactors = FALSE
    )
    names(pvals)[1] <- merge.by
    list(
        pvals = pvals,
        vcov  = diag(0.01, length(ids))
    )
}

# =============================================================================
# 1.  .qtlSelect — mbf path (attr(intervalObj, "env") is NULL)
# =============================================================================

test_that(".qtlSelect mbf path: computes oint, blups, state correctly", {
    set.seed(1)
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 5)
    n_int  <- 10  # 2 chr x 5 markers
    all_keys <- paste("Chr", rep(c("C1","C2"), each = 5),
                      rep(1:5, times = 2), sep = ".")
    state <- rep(1, n_int); names(state) <- all_keys

    asm <- make_updated_model(n_int = n_int)

    res <- wgAim:::.qtlSelect(
        asm        = asm,
        phenoData  = data.frame(id = factor(paste0("L", 1:20)), yld = rnorm(20)),
        intervalObj = genObj,
        gen.type   = "marker",
        selection  = "interval",
        exclusion.window = 20,
        state      = state,
        verboseLev = 0
    )

    expect_named(res, c("state", "qtl", "ochr", "oint", "blups"))
    expect_equal(length(res$oint), n_int)
    expect_equal(length(res$blups), n_int)
    expect_true(is.character(res$qtl))
    expect_true(all(res$state %in% c(0, 1)))
})

test_that(".qtlSelect mbf path with verboseLev=1 prints output", {
    set.seed(2)
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 5)
    all_keys <- paste("Chr", rep(c("C1","C2"), each = 5),
                      rep(1:5, times = 2), sep = ".")
    state <- rep(1, 10); names(state) <- all_keys
    asm <- make_updated_model(n_int = 10)

    expect_output(
        wgAim:::.qtlSelect(asm, phenoData = data.frame(id = factor(paste0("L", 1:20)),
                                                        yld = rnorm(20)),
                            genObj, "marker", "interval", 20, state, verboseLev = 1),
        "Outlier"
    )
})

test_that(".qtlSelect chromosome selection path works", {
    set.seed(3)
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 5)
    all_keys <- paste("Chr", rep(c("C1","C2"), each = 5),
                      rep(1:5, times = 2), sep = ".")
    state <- rep(1, 10); names(state) <- all_keys
    asm <- make_updated_model(n_int = 10)

    res <- wgAim:::.qtlSelect(asm,
                               phenoData  = data.frame(id = factor(paste0("L", 1:20)),
                                                       yld = rnorm(20)),
                               intervalObj = genObj,
                               gen.type   = "marker",
                               selection  = "chromosome",
                               exclusion.window = 20,
                               state      = state,
                               verboseLev = 0)

    expect_false(is.null(res$ochr))
    expect_equal(length(res$ochr), 2L)
})

test_that(".qtlSelect chromosome selection with verboseLev=1 prints output", {
    set.seed(4)
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 5)
    all_keys <- paste("Chr", rep(c("C1","C2"), each = 5),
                      rep(1:5, times = 2), sep = ".")
    state <- rep(1, 10); names(state) <- all_keys
    asm <- make_updated_model(n_int = 10)

    expect_output(
        wgAim:::.qtlSelect(asm, phenoData = data.frame(id = factor(paste0("L", 1:20)),
                                                        yld = rnorm(20)),
                            genObj, "marker", "chromosome", 20, state, verboseLev = 1),
        "chromosome"
    )
})

test_that(".qtlSelect sigma2=1 branch fires when vpar.con==4", {
    set.seed(5)
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 5)
    all_keys <- paste("Chr", rep(c("C1","C2"), each = 5),
                      rep(1:5, times = 2), sep = ".")
    state <- rep(1, 10); names(state) <- all_keys
    asm <- make_updated_model(n_int = 10, sigma2 = 99)
    # Force vpar.con last element to 4 (boundary constraint code)
    asm$vparameters.con[length(asm$vparameters.con)] <- 4L

    # Should not error — sigma2 is silently set to 1 internally
    res <- wgAim:::.qtlSelect(asm, data.frame(id = factor(paste0("L", 1:20)),
                                               yld = rnorm(20)),
                               genObj, "marker", "interval", 20, state, 0)
    expect_true(is.list(res))
})

test_that(".qtlSelect interval gen.type uses inferred.map for exclusion", {
    set.seed(6)
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 5)
    all_keys <- paste("Chr", rep(c("C1","C2"), each = 5),
                      rep(1:5, times = 2), sep = ".")
    state <- rep(1, 10); names(state) <- all_keys
    asm <- make_updated_model(n_int = 10)

    res <- wgAim:::.qtlSelect(asm, data.frame(id = factor(paste0("L", 1:20)),
                                               yld = rnorm(20)),
                               genObj, "interval", "interval", 20, state, 0)
    expect_true(is.list(res))
    expect_true(all(res$state %in% c(0, 1)))
})

# =============================================================================
# 2.  .qtlSelect — vm path (attr(intervalObj, "env") is set)
# =============================================================================

test_that(".qtlSelect vm path: uses predict() and trans matrix", {
    set.seed(7)
    n_lines <- 10; n_int <- 6
    genObj  <- make_wgCross_interval(n_lines = n_lines, n_chr = 2, n_mar = 3)
    ids     <- paste0("L", 1:n_lines)

    # Build a cov.env with correct dimensions
    genoData <- do.call(cbind, lapply(genObj$geno, `[[`, "imputed.data"))
    cov.env  <- wgAim:::.constructCM(genoData)

    # Attach env attribute to genObj (vm path signal)
    attr(genObj, "env") <- cov.env

    all_keys <- paste("Chr", rep(c("C1","C2"), each = 3),
                      rep(1:3, times = 2), sep = ".")
    state <- rep(1, n_int); names(state) <- all_keys

    asm <- make_updated_model(n_int = n_int)
    # Give asm a vm-style vparameter name
    vm_name <- paste0("vm(id, covObj)!vm!var")
    names(asm$vparameters)[1] <- vm_name
    names(asm$vparameters.con)[1] <- vm_name

    # Give asm a random formula with vm term so grep finds it
    asm$call$random <- ~ vm(id, covObj) + id

    pv_ret <- make_predict_return(ids = ids)

    with_mocked_bindings(
        predict = function(object, classify, only, vcov, data, ...) pv_ret,
        .package = "wgAim",
        {
            res <- wgAim:::.qtlSelect(
                asm        = asm,
                phenoData  = data.frame(id = factor(ids), yld = rnorm(n_lines)),
                intervalObj = genObj,
                gen.type   = "marker",
                selection  = "interval",
                exclusion.window = 10,
                state      = state,
                verboseLev = 0
            )
        }
    )
    expect_named(res, c("state", "qtl", "ochr", "oint", "blups"))
    expect_equal(length(res$oint), n_int)
})

# =============================================================================
# 3.  .buildGenomeModel — mbf path (lines >= markers)
# =============================================================================

test_that(".buildGenomeModel mbf path: assigns covObj, returns qtlModel list", {
    set.seed(8)
    n_lines <- 20; n_int <- 6
    genObj  <- make_wgCross_interval(n_lines = n_lines, n_chr = 2, n_mar = 3)
    ids     <- paste0("L", 1:n_lines)
    genoData <- do.call(cbind, lapply(genObj$geno, `[[`, "imputed.data"))
    rownames(genoData) <- ids

    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines),
                             stringsAsFactors = FALSE)
    base_m <- make_updated_model(n_int = n_int)
    base_m$call <- bquote(asreml(fixed = yld ~ 1, random = ~ id, data = .(phenoData)))

    updated_m <- make_updated_model(n_int = n_int, loglik = -48)
    caller_e  <- new.env(parent = emptyenv())

    with_mocked_bindings(
        update = function(object, ...) updated_m,
        .package = "wgAim",
        {
            res <- wgAim:::.buildGenomeModel(
                baseModel  = base_m,
                genoData   = genoData,   # 20 lines x 6 markers: lines > markers -> mbf
                phenoData  = phenoData,
                merge.by   = "id",
                intervalObj = genObj,
                force      = FALSE,
                rterms     = character(0),
                caller.env = caller_e
            )
        }
    )
    expect_named(res, c("qtlModel", "intervalObj", "cov.env", "vm", "vmterms",
                        "n.fa", "vm.struct", "upgrade.to.fa", "gterm", "resid.term"),
                 ignore.order = TRUE)
    expect_false(res$vm)
    expect_null(res$cov.env)
    expect_true(exists("covObj", envir = caller_e))
})

test_that(".buildGenomeModel vm path: builds GRM, sets vm=TRUE", {
    set.seed(9)
    # markers (10) > lines (6) -> vm path
    n_lines <- 6; n_mar <- 5; n_chr <- 2; n_int <- n_mar * n_chr
    genObj  <- make_wgCross_interval(n_lines = n_lines, n_chr = n_chr, n_mar = n_mar)
    ids     <- paste0("L", 1:n_lines)
    genoData <- do.call(cbind, lapply(genObj$geno, `[[`, "imputed.data"))
    rownames(genoData) <- ids

    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines),
                             stringsAsFactors = FALSE)
    base_m <- make_updated_model(n_int = n_int)
    base_m$call <- bquote(asreml(fixed = yld ~ 1, random = ~ id, data = .(phenoData)))

    updated_m <- make_updated_model(n_int = n_int, loglik = -45)
    caller_e  <- new.env(parent = emptyenv())

    with_mocked_bindings(
        update = function(object, ...) updated_m,
        .package = "wgAim",
        {
            res <- wgAim:::.buildGenomeModel(
                baseModel  = base_m,
                genoData   = genoData,   # 10 markers > 6 lines -> vm
                phenoData  = phenoData,
                merge.by   = "id",
                intervalObj = genObj,
                force      = FALSE,
                rterms     = character(0),
                caller_e
            )
        }
    )
    expect_true(res$vm)
    expect_false(is.null(res$cov.env))
})

# =============================================================================
# 4.  .rebuildCovObj — vm->mbf switch path
# =============================================================================

test_that(".rebuildCovObj mbf path: returns correct structure", {
    set.seed(10)
    n_lines <- 20; n_int <- 6
    genObj  <- make_wgCross_interval(n_lines = n_lines, n_chr = 2, n_mar = 3)
    ids     <- paste0("L", 1:n_lines)
    genoData <- do.call(cbind, lapply(genObj$geno, `[[`, "imputed.data"))
    rownames(genoData) <- ids
    all_keys <- colnames(genoData)
    state <- rep(1, n_int); names(state) <- all_keys

    qtlModel <- make_updated_model(n_int = n_int)
    caller_e <- new.env(parent = emptyenv())

    res <- wgAim:::.rebuildCovObj(
        genoData   = genoData,
        state      = state,
        merge.by   = "id",
        intervalObj = genObj,
        force      = FALSE,
        vm         = FALSE,
        vmterms    = NULL,
        qtlModel   = qtlModel,
        caller.env = caller_e
    )

    expect_named(res, c("cov.env", "intervalObj", "qtlModel"))
    expect_true(exists("covObj", envir = caller_e))
})

test_that(".rebuildCovObj vm->mbf switch: fires when ncol<=nrow after exclusion", {
    set.seed(11)
    # Start with markers > lines so vm=TRUE; after excluding some, lines >= markers
    n_lines <- 6; n_int <- 8
    genObj  <- make_wgCross_interval(n_lines = n_lines, n_chr = 2, n_mar = 4)
    ids     <- paste0("L", 1:n_lines)
    genoData <- do.call(cbind, lapply(genObj$geno, `[[`, "imputed.data"))
    rownames(genoData) <- ids
    all_keys <- colnames(genoData)

    # Exclude 5 columns so only 3 remain — 3 < 6 lines -> mbf path
    state <- rep(0, n_int); names(state) <- all_keys
    state[1:3] <- 1   # keep only 3 columns

    qtlModel <- make_updated_model(n_int = n_int)
    qtlModel$call$random <- ~ vm(id, covObj) + id  # was vm path
    caller_e <- new.env(parent = emptyenv())

    res <- wgAim:::.rebuildCovObj(
        genoData   = genoData,
        state      = state,
        merge.by   = "id",
        intervalObj = genObj,
        force      = FALSE,
        vm         = TRUE,      # was vm, now should switch
        vmterms    = c("vm(id, covObj)", "id"),
        qtlModel   = qtlModel,
        caller.env = caller_e
    )

    expect_null(res$cov.env)  # switched to mbf so cov.env is NULL
})

# =============================================================================
# 5.  .vModify — resets G.param constraints
# =============================================================================

test_that(".vModify: modifies G.param con from B to P for matching terms", {
    model <- make_updated_model()
    # Add a G.param entry that matches the mbf pattern
    model$G.param <- list(
        `mbf(ints)_id!mbf(ints)!var` = list(
            list(con = c("B", "U"), initial = c(0.0, 0.5))
        )
    )
    result <- wgAim:::.vModify(model, "id")
    expect_equal(result$G.param[[1]][[1]]$con[1], "P")
    expect_equal(result$G.param[[1]][[1]]$initial[1], 0.1)
})

test_that(".vModify: leaves non-matching G.param terms unchanged", {
    model <- make_updated_model()
    model$G.param <- list(
        `someOtherTerm` = list(list(con = c("U"), initial = c(0.5)))
    )
    result <- wgAim:::.vModify(model, "id")
    expect_equal(result$G.param[[1]][[1]]$con[1], "U")
})

# =============================================================================
# 6.  .envFix — cleans up formula environments
# =============================================================================

test_that(".envFix: returns model with cleaned environments", {
    model <- make_updated_model()
    # Add the minimum structure .envFix expects
    e <- new.env()
    attr(model$formulae$fixed, ".Environment")  <- e
    attr(model$formulae$random, ".Environment") <- e
    model$call$fixed  <- quote(yld ~ 1)
    model$call$random <- quote(~ id)
    environment(model$call$fixed)  <- e
    environment(model$call$random) <- e
    attr(model$mf, "model.terms") <- list()

    asremlEnv <- list(fixed = e, random = e)
    result <- wgAim:::.envFix(model, asremlEnv)
    expect_true(is.list(result))
})

# =============================================================================
# 7.  .addEffect — fixed and random paths
# =============================================================================

test_that(".addEffect fixed path: updates both models and returns coefs", {
    set.seed(12)
    n_lines <- 20
    ids     <- paste0("L", 1:n_lines)
    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines),
                             `X.C1.2` = rnorm(n_lines),
                             stringsAsFactors = FALSE, check.names = FALSE)

    base_m  <- make_updated_model(n_int = 6)
    qtl_m   <- make_updated_model(n_int = 6)

    # Give fixed coefficients an X. entry so grep finds it
    fix_coefs <- matrix(c(5.0, 0.4), 2, 1,
                        dimnames = list(c("(Intercept)", "X.C1.2"), "effect"))
    qtl_m$coefficients$fixed  <- fix_coefs
    qtl_m$vcoeff$fixed        <- c("(Intercept)" = 0.1, "X.C1.2" = 0.02)

    updated <- qtl_m  # update() returns same structure

    with_mocked_bindings(
        update = function(object, ...) updated,
        .package = "wgAim",
        {
            res <- wgAim:::.addEffect(
                baseModel = base_m,
                qtlModel  = qtl_m,
                phenoData = phenoData,
                merge.by  = "id",
                qtl.x     = "X.C1.2",
                method    = "fixed",
                iter      = 1L
            )
        }
    )
    expect_named(res, c("baseModel", "qtlModel", "coefs", "vcoefs"))
    expect_true(is.numeric(res$coefs))
    expect_true("X.C1.2" %in% names(res$coefs))
})

test_that(".addEffect random path: updates both models and returns coefs", {
    set.seed(13)
    n_lines <- 20
    ids     <- paste0("L", 1:n_lines)
    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines),
                             `X.C1.2` = rnorm(n_lines),
                             stringsAsFactors = FALSE, check.names = FALSE)

    base_m  <- make_updated_model(n_int = 6)
    qtl_m   <- make_updated_model(n_int = 6)

    # Give random coefficients an X. entry
    all_rnam <- c(paste0("mbf_", 1:6), "X.C1.2_L1", "X.C1.2_L2")
    rand_mat <- matrix(rnorm(length(all_rnam)), length(all_rnam), 1,
                       dimnames = list(all_rnam, "effect"))
    qtl_m$coefficients$random <- rand_mat
    qtl_m$vcoeff$random       <- setNames(runif(length(all_rnam), 0.01, 0.1), all_rnam)

    base_m$G.param <- qtl_m$G.param <- list(
        `mbf(ints)_id!mbf(ints)!var` = list(list(con = "B", initial = 0))
    )
    updated <- qtl_m

    with_mocked_bindings(
        update = function(object, ...) updated,
        .package = "wgAim",
        {
            res <- wgAim:::.addEffect(
                baseModel = base_m,
                qtlModel  = qtl_m,
                phenoData = phenoData,
                merge.by  = "id",
                qtl.x     = "X.C1.2",
                method    = "random",
                iter      = 1L
            )
        }
    )
    expect_named(res, c("baseModel", "qtlModel", "coefs", "vcoefs"))
})

# =============================================================================
# 8.  qtlAim / gwasAim full loop — mocked to run one iteration then break
# =============================================================================

# Build a complete mock environment to drive one qtlAim iteration.
# We need: .validateModel, .buildGenoData, .fixLines, .buildGenomeModel,
#          .qtlSelect, .lrtTest, .packResults to all return plausible values,
# and update() to return a minimal model.

make_full_mock_env <- function(n_lines = 20, n_chr = 2, n_mar = 3,
                                lrt_pass = FALSE) {
    set.seed(42)
    ids      <- paste0("L", 1:n_lines)
    genObj   <- make_wgCross_interval(n_lines, n_chr, n_mar)
    genObj$pheno$id <- factor(ids)
    n_int    <- n_chr * n_mar
    all_keys <- paste("Chr", rep(paste0("C", 1:n_chr), each = n_mar),
                      rep(1:n_mar, times = n_chr), sep = ".")
    genoData <- do.call(cbind, lapply(genObj$geno, `[[`, "imputed.data"))
    rownames(genoData) <- ids
    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines),
                             stringsAsFactors = FALSE)

    base_m <- make_updated_model(n_int = n_int)
    base_m$call <- bquote(asreml(fixed = yld ~ 1, random = ~ id, data = .(phenoData)))

    upd_m <- make_updated_model(n_int = n_int,
                                 loglik = if (lrt_pass) -44 else -52)

    state <- rep(1, n_int); names(state) <- all_keys

    list(
        genObj    = genObj,
        phenoData = phenoData,
        base_m    = base_m,
        upd_m     = upd_m,
        genoData  = genoData,
        all_keys  = all_keys,
        state     = state,
        ids       = ids,
        n_int     = n_int
    )
}

test_that("qtlAim full loop: LRT fails -> breakout -> packResults called", {
    e <- make_full_mock_env(lrt_pass = FALSE)

    caller_e <- new.env(parent = globalenv())
    caller_e$phenoData <- e$phenoData
    caller_e$genObj    <- e$genObj
    bm <- e$base_m; class(bm) <- c("asreml", "list")
    caller_e$base_m    <- bm    # needs "asreml" class for S3 dispatch

    with_mocked_bindings(
        .validateModel = function(baseModel, ...) {
            list(baseModel = baseModel,
                 asremlEnv = list(fixed = new.env(), random = new.env()),
                 phenoData = e$phenoData)
        },
        .buildGenoData = function(intervalObj, gen.type, glines, plines) {
            gd <- e$genoData[as.character(glines) %in% as.character(plines), , drop = FALSE]
            list(genoData = gd, mnams = colnames(gd),
                 state = setNames(rep(1, ncol(gd)), colnames(gd)))
        },
        .fixLines = function(baseModel, phenoData, genoData, merge.by,
                              plines, fix.lines, ...) {
            list(baseModel = baseModel, phenoData = phenoData,
                 merge.by = merge.by, rterms = character(0),
                 genetic.term = merge.by)
        },
        .buildGenomeModel = function(...) {
            list(qtlModel = e$upd_m, intervalObj = e$genObj,
                 cov.env = NULL, vm = FALSE, vmterms = NULL)
        },
        .qtlSelect = function(...) {
            oint <- setNames(runif(e$n_int), e$all_keys)
            list(state = e$state, qtl = e$all_keys[1],
                 ochr = NULL, oint = oint,
                 blups = setNames(rnorm(e$n_int), e$all_keys))
        },
        .lrtTest = function(qtlModel, baseModel, TypeI, ...) {
            list(baseLogL = -52, stat = 1.0, pvalue = 0.3, pass = FALSE)
        },
        .envFix = function(model, asremlEnv) model,
        update  = function(object, ...) e$upd_m,
        .package = "wgAim",
        {
            result <- eval(quote(
                qtlAim(base_m, genObj = genObj,
                        merge.by = "id", trace = FALSE)
            ), envir = caller_e)
        }
    )

    expect_s3_class(result, "qtlAim")
    expect_equal(result$QTL$iterations, 1L)
    expect_null(result$QTL$effects)  # no QTL found (LRT failed)
})

test_that("gwasAim full loop: LRT fails -> packResults with n.markers set", {
    e <- make_full_mock_env(lrt_pass = FALSE)
    panel <- make_wgPanel(n_lines = 20, n_chr = 2, n_mar = 3)
    panel$pheno$id <- factor(e$ids)

    caller_e <- new.env(parent = globalenv())
    caller_e$phenoData <- e$phenoData
    caller_e$panel     <- panel
    bm <- e$base_m; class(bm) <- c("asreml", "list")
    caller_e$base_m    <- bm    # needs "asreml" class for S3 dispatch

    with_mocked_bindings(
        .validateModel = function(baseModel, ...) {
            list(baseModel = baseModel,
                 asremlEnv = list(fixed = new.env(), random = new.env()),
                 phenoData = e$phenoData)
        },
        .buildGenoData = function(intervalObj, gen.type, glines, plines) {
            n <- e$n_int
            gd <- matrix(sample(c(-1,0,1), 20*n, replace=TRUE), 20, n,
                         dimnames = list(e$ids, e$all_keys))
            list(genoData = gd, mnams = e$all_keys,
                 state = setNames(rep(1, n), e$all_keys))
        },
        .fixLines = function(baseModel, phenoData, genoData, merge.by,
                              plines, fix.lines, ...) {
            list(baseModel = baseModel, phenoData = phenoData,
                 merge.by = merge.by, rterms = character(0),
                 genetic.term = merge.by)
        },
        .buildGenomeModel = function(...) {
            list(qtlModel = e$upd_m, intervalObj = panel,
                 cov.env = NULL, vm = FALSE, vmterms = NULL)
        },
        .qtlSelect = function(...) {
            oint <- setNames(runif(e$n_int), e$all_keys)
            list(state = e$state, qtl = e$all_keys[1],
                 ochr = NULL, oint = oint,
                 blups = setNames(rnorm(e$n_int), e$all_keys))
        },
        .lrtTest = function(...) list(baseLogL=-60, stat=0.5, pvalue=0.5, pass=FALSE),
        .envFix  = function(model, asremlEnv) model,
        update   = function(object, ...) e$upd_m,
        .package  = "wgAim",
        {
            result <- eval(quote(
                gwasAim(base_m, genObj = panel,
                         merge.by = "id", trace = FALSE)
            ), envir = caller_e)
        }
    )

    expect_s3_class(result, "gwasAim")
    expect_equal(result$QTL$n.markers, e$n_int)
})

# =============================================================================
# 9.  fineMap — qtlAim path with mocked update()
# =============================================================================

test_that("fineMap qtlAim interval path: runs scan loop with mocked update()", {
    set.seed(20)
    qtl_obj <- make_mock_qtlAim(n_qtl = 1, n_lines = 20, n_chr = 2, n_mar = 4)
    genObj  <- attr(qtl_obj, "genObj")
    ids     <- paste0("L", 1:20)

    # fineMap retrieves phenoData from formula environment.
    # $call$fixed must be a proper formula (not a call) so terms() works on it.
    phenoData <- qtl_obj$mf
    fe <- new.env(parent = globalenv())
    assign("yld.data", phenoData, envir = fe)
    qtl_obj$call$fixed <- yld ~ 1
    environment(qtl_obj$call$fixed) <- fe

    n_int <- 8  # 2 chr x 4 mar
    upd_m <- make_updated_model(n_int = n_int, loglik = -45)
    upd_m$coefficients$fixed <- matrix(
        c(5.0, 0.3), 2, 1,
        dimnames = list(c("(Intercept)", "X.C1.2_fm"), "effect")
    )
    upd_m$vcoeff$fixed <- c("(Intercept)" = 0.1, "X.C1.2_fm" = 0.02)

    with_mocked_bindings(
        update = function(object, ...) upd_m,
        .package = "wgAim",
        {
            fm <- fineMap(qtl_obj, genObj = genObj,
                          window = 40, step = 20, exclusion.window = 10)
        }
    )

    expect_s3_class(fm, "fineMap")
    expect_true(length(fm) >= 1L)
    expect_true(is.data.frame(fm[[1]]))
    expect_named(fm[[1]], c("mark", "dist", "pvalue", "LOD"))
})

test_that("fineMap gwasAim path: runs marker scan with mocked update()", {
    set.seed(21)
    gwas_obj <- make_mock_gwasAim(n_qtl = 1, n_lines = 50, n_chr = 2, n_mar = 10)
    genObj   <- attr(gwas_obj, "genObj")
    ids      <- paste0("S", 1:50)

    phenoData <- gwas_obj$mf
    fe <- new.env(parent = globalenv())
    assign("yld.data", phenoData, envir = fe)
    gwas_obj$call$fixed <- yld ~ 1
    environment(gwas_obj$call$fixed) <- fe

    n_mar_total <- 20
    upd_m <- make_updated_model(n_int = n_mar_total, loglik = -55)
    upd_m$coefficients$fixed <- matrix(
        c(10.0, 0.2), 2, 1,
        dimnames = list(c("(Intercept)", "X.P1.3_fm"), "effect")
    )
    upd_m$vcoeff$fixed <- c("(Intercept)" = 0.1, "X.P1.3_fm" = 0.015)

    with_mocked_bindings(
        update = function(object, ...) upd_m,
        .package = "wgAim",
        {
            fm <- fineMap(gwas_obj, genObj = genObj,
                          window = 50, exclusion.window = 10)
        }
    )

    expect_s3_class(fm, "fineMap")
    expect_true(is.data.frame(fm[[1]]))
})

test_that("fineMap: qtl= subset argument returns only requested QTL", {
    set.seed(22)
    qtl_obj <- make_mock_qtlAim(n_qtl = 2, n_lines = 20, n_chr = 2, n_mar = 4)
    genObj  <- attr(qtl_obj, "genObj")

    phenoData <- qtl_obj$mf
    fe <- new.env(parent = globalenv())
    assign("yld.data", phenoData, envir = fe)
    qtl_obj$call$fixed <- yld ~ 1
    environment(qtl_obj$call$fixed) <- fe

    n_int <- 8
    upd_m <- make_updated_model(n_int = n_int, loglik = -45)
    upd_m$coefficients$fixed <- matrix(c(5.0, 0.3), 2, 1,
        dimnames = list(c("(Intercept)", "X.C1.2_fm"), "effect"))
    upd_m$vcoeff$fixed <- c("(Intercept)" = 0.1, "X.C1.2_fm" = 0.02)

    # Request only the first QTL by key
    first_key <- qtl_obj$QTL$qtl[1]

    with_mocked_bindings(
        update = function(object, ...) upd_m,
        .package = "wgAim",
        {
            fm <- fineMap(qtl_obj, genObj = genObj,
                          qtl = first_key, window = 40, step = 20)
        }
    )

    expect_s3_class(fm, "fineMap")
    expect_equal(length(fm), 1L)
    expect_equal(names(fm), first_key)
})
