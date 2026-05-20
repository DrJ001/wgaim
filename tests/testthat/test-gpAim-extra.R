# =============================================================================
# test-gpAim-extra.R
# Additional tests for gpAim.R, summary.gpAim.R, and engine_genodata.R
# targeting the branches not covered in test-gpAim.R:
#
#  gpAim.R:
#   - Trait guards (not a factor, <2 levels, Trait not found)
#   - str guards (fa no number, fa too large)
#   - str accepted with Trait (reaches update)
#   - .cullis_H2 helper
#   - .extract_Ga_gp helper
#
#  summary.gpAim.R:
#   - MV summary path (Trait non-NULL)
#   - gen.H2 present vs absent in MV
#   - Accuracy column rounding
#
#  engine_genodata.R:
#   - .buildGenoData: interval vs marker mode, column naming
#   - .fixLines: fix.lines=FALSE path (no update called)
#   - .fixLines: fix.lines=TRUE path basics
# =============================================================================

.cullis_H2    <- wgAim:::.cullis_H2
.extract_Ga_gp <- wgAim:::.extract_Ga_gp

# ---- shared helper: MV gpAim mock with GP$Trait populated -------------------
.make_mv_gpAim_obj <- function(n_lines = 20, trials = c("T1", "T2"),
                                gen.H2 = setNames(c(0.80, 0.72), c("T1", "T2")),
                                seed = 5) {
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
        Accuracy = sqrt(pmax(0, 1 - SE_mv^2 / rep(Ga_diag, each = n_lines))),
        gen.H2 = if (!is.null(gen.H2))
            rep(gen.H2[trials], each = n_lines) else NA_real_,
        stringsAsFactors = FALSE
    )

    GP <- list(
        gebv         = gebv_df,
        gen.type     = "marker",
        path         = "vm",
        var.genetic  = setNames(Ga_diag, trials),
        var.resid    = setNames(rep(1, nt), trials),
        heritability = setNames(Ga_diag / (Ga_diag + 1), trials),
        gen.H2       = gen.H2,
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
# 1.  gpAim Trait guards
# =============================================================================

make_gp_model_mv <- function(ids = paste0("L", 1:20),
                               trials = c("T1", "T2")) {
    df <- data.frame(
        id    = factor(rep(ids, length(trials))),
        Trial = factor(rep(trials, each = length(ids))),
        yld   = rnorm(length(ids) * length(trials)),
        stringsAsFactors = FALSE
    )
    m <- list(
        converge        = TRUE,
        loglik          = -80,
        sigma2          = 1,
        vparameters     = c(Trial.id = 0.4, "R!variance" = 1),
        vparameters.con = c(0L, 0L),
        coefficients    = list(
            fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
            random = matrix(rnorm(length(ids)), length(ids), 1,
                            dimnames = list(paste0("id_", ids), "effect"))),
        formulae = list(fixed  = yld ~ Trial,
                        random = as.formula("~ corh(Trial):id")),
        mf   = df,
        call = bquote(asreml(fixed   = yld ~ Trial,
                             random  = ~ corh(Trial):id,
                             data    = .(df)))
    )
    class(m) <- c("asreml", "list")
    m
}

test_that("gpAim: Trait column not found in phenoData stops with error", {
    m     <- make_gp_model_mv()
    gen_i <- make_wgCross_interval(n_lines = 20)
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id", Trait = "NOTEXIST"),
        regexp = "not found"
    )
})

test_that("gpAim: Trait with only 1 level stops with error", {
    df <- data.frame(id = factor(paste0("L", 1:20)),
                     Trial = factor(rep("T1", 20)),
                     yld = rnorm(20))
    m <- list(
        converge = TRUE, loglik = -80, sigma2 = 1,
        vparameters = c(id = 0.4, "R!variance" = 1),
        vparameters.con = c(0L, 0L),
        coefficients = list(
            fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
            random = matrix(rnorm(20), 20, 1, dimnames = list(paste0("id_L", 1:20), "e"))),
        formulae = list(fixed = yld ~ 1, random = ~ id),
        mf = df,
        call = bquote(asreml(fixed = yld ~ 1, random = ~ id, data = .(df)))
    )
    class(m) <- c("asreml", "list")
    gen_i <- make_wgCross_interval(n_lines = 20)
    gen_i$pheno$id <- factor(paste0("L", seq_len(20)))
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id", Trait = "Trial"),
        regexp = "at least 2 levels"
    )
})

test_that("gpAim: str='fa' (no number) stops with error", {
    m     <- make_gp_model_mv()
    gen_i <- make_wgCross_interval(n_lines = 20)
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id", Trait = "Trial", str = "fa"),
        regexp = "number of factors"
    )
})

test_that("gpAim: str='fa3' too large for 2 traits stops with error", {
    m     <- make_gp_model_mv()
    gen_i <- make_wgCross_interval(n_lines = 20)
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id", Trait = "Trial", str = "fa3"),
        regexp = "too large|exceeds unstructured"
    )
})

test_that("gpAim: str='fa1' too large for 2 traits stops with error", {
    # ntrait=2: n.par.us = 3; fa1: (1+1)*2 - 1*0/2 = 4 > 3 → error
    m     <- make_gp_model_mv()
    gen_i <- make_wgCross_interval(n_lines = 20)
    gen_i$pheno$id <- factor(paste0("L", seq_len(nrow(gen_i$pheno))))
    expect_error(
        gpAim(m, genObj = gen_i, merge.by = "id",
              Trait = "Trial", str = "fa1"),
        regexp = "too large|exceeds unstructured"
    )
})

# =============================================================================
# 2.  .cullis_H2 internal helper
# =============================================================================

test_that(".cullis_H2: returns NA for NULL sed_sub", {
    expect_equal(.cullis_H2(NULL, 1.0), NA_real_)
})

test_that(".cullis_H2: returns NA for NA G_jj", {
    expect_equal(.cullis_H2(diag(3), NA_real_), NA_real_)
})

test_that(".cullis_H2: returns NA for G_jj <= 0", {
    expect_equal(.cullis_H2(diag(3), 0), NA_real_)
    expect_equal(.cullis_H2(diag(3), -1), NA_real_)
})

test_that(".cullis_H2: returns value in [0,1] for valid SED matrix", {
    set.seed(1)
    n   <- 10
    sed <- matrix(0.1, n, n); diag(sed) <- 0
    val <- .cullis_H2(sed, 0.8)
    expect_gte(val, 0)
    expect_lte(val, 1)
})

test_that(".cullis_H2: returns NA when no finite positive upper-triangle values", {
    sed <- matrix(c(0, -1, -1, 0), 2, 2)
    expect_equal(.cullis_H2(sed, 1.0), NA_real_)
})

# =============================================================================
# 3.  .extract_Ga_gp internal helper
# =============================================================================

.make_gpModel_diag <- function(ntrait = 2, sigma2 = 1.0) {
    tl   <- paste0("T", seq_len(ntrait))
    vpar <- setNames(c(rep(0.5, ntrait), 1.0),
                     c(paste0("Trial:vm(id, covObj)!Trial_", seq_len(ntrait)), "R!variance"))
    g <- list(
        vparameters = vpar,
        sigma2      = sigma2,
        G.param     = setNames(list(list()), "Trial:vm(id, covObj)")
    )
    class(g) <- c("asreml", "list")
    g
}

.make_gpModel_corh <- function(ntrait = 2, sigma2 = 1.0) {
    tl   <- paste0("T", seq_len(ntrait))
    # corh stores: correlations first (ending .cor), then variances
    cor.nms <- "corh(Trial):vm(id, covObj)!Trial.cor"
    var.nms <- paste0("corh(Trial):vm(id, covObj)!Trial_", seq_len(ntrait))
    vpar    <- setNames(c(0.6, rep(0.5, ntrait), 1.0),
                        c(cor.nms, var.nms, "R!variance"))
    g <- list(
        vparameters = vpar,
        sigma2      = sigma2,
        G.param     = setNames(list(list()), "Trial:vm(id, covObj)")
    )
    class(g) <- c("asreml", "list")
    g
}

test_that(".extract_Ga_gp: diag structure gives diagonal matrix", {
    gm <- .make_gpModel_diag(ntrait = 2)
    Ga <- .extract_Ga_gp(gm, 1.0, 2L, c("T1", "T2"))
    expect_equal(dim(Ga), c(2L, 2L))
    expect_equal(Ga[1, 2], 0)
    expect_equal(Ga[2, 1], 0)
})

test_that(".extract_Ga_gp: diag variances scaled by sigma2", {
    gm <- .make_gpModel_diag(ntrait = 2, sigma2 = 2.0)
    Ga <- .extract_Ga_gp(gm, 2.0, 2L, c("T1", "T2"))
    # vpar = 0.5, sigma2 = 2.0 → Ga diag = 1.0
    expect_equal(diag(Ga), c(T1 = 1.0, T2 = 1.0))
})

test_that(".extract_Ga_gp: corh structure returns symmetric matrix", {
    gm <- .make_gpModel_corh(ntrait = 2)
    Ga <- .extract_Ga_gp(gm, 1.0, 2L, c("T1", "T2"))
    expect_equal(dim(Ga), c(2L, 2L))
    expect_equal(Ga, t(Ga))
})

test_that(".extract_Ga_gp: dimnames match trait levels", {
    gm <- .make_gpModel_diag(ntrait = 3)
    tl <- c("E1", "E2", "E3")
    Ga <- .extract_Ga_gp(gm, 1.0, 3L, tl)
    expect_equal(rownames(Ga), tl)
    expect_equal(colnames(Ga), tl)
})

# =============================================================================
# 4.  summary.gpAim — MV path
# =============================================================================

test_that("summary.gpAim MV: returns data.frame", {
    obj <- .make_mv_gpAim_obj()
    tab <- capture.output(result <- summary(obj))
    expect_s3_class(result, "data.frame")
})

test_that("summary.gpAim MV: output contains 'Multivariate'", {
    obj <- .make_mv_gpAim_obj()
    out <- capture.output(summary(obj))
    expect_true(any(grepl("Multivariate", out, ignore.case = TRUE)))
})

test_that("summary.gpAim MV: Accuracy column rounded", {
    obj    <- .make_mv_gpAim_obj()
    capture.output(result <- summary(obj))
    expect_true("Accuracy" %in% names(result))
    # Rounded to 4 decimal places
    expect_equal(result$Accuracy, round(result$Accuracy, 4))
})

test_that("summary.gpAim MV: gen.H2 shown when present", {
    obj <- .make_mv_gpAim_obj(gen.H2 = setNames(c(0.80, 0.72), c("T1", "T2")))
    out <- capture.output(summary(obj))
    expect_true(any(grepl("H2|h2", out, ignore.case = TRUE)))
})

test_that("summary.gpAim MV: no gen.H2 block when gen.H2 is NULL (all NA)", {
    obj <- .make_mv_gpAim_obj(gen.H2 = NULL)
    # GP$gen.H2 = NULL → summary skips the generalised H2 block
    obj$GP$gen.H2 <- NULL
    out <- capture.output(summary(obj))
    expect_false(any(grepl("Generalised H2", out)))
})

test_that("summary.gpAim MV: trial names appear in output", {
    obj <- .make_mv_gpAim_obj(trials = c("Site1", "Site2"))
    out <- capture.output(summary(obj))
    expect_true(any(grepl("Site1", out)))
    expect_true(any(grepl("Site2", out)))
})

test_that("summary.gpAim MV: rows sorted by trial then descending GEBV", {
    obj <- .make_mv_gpAim_obj()
    capture.output(result <- summary(obj))
    for (t in c("T1", "T2")) {
        sub_t <- result[as.character(result$Trial) == t, ]
        expect_true(all(diff(sub_t$GEBV) <= 0))
    }
})

# =============================================================================
# 5.  engine_genodata.R — .buildGenoData
# =============================================================================

.buildGenoData <- wgAim:::.buildGenoData

test_that(".buildGenoData: interval mode uses interval.data", {
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 4)
    glines <- paste0("L", 1:20)
    plines <- glines
    out    <- .buildGenoData(genObj, "interval", glines, plines)
    # Column names should encode interval indices
    expect_true(all(grepl("^Chr\\.", colnames(out$genoData))))
})

test_that(".buildGenoData: marker mode uses imputed.data", {
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 4)
    glines <- paste0("L", 1:20)
    plines <- glines
    out    <- .buildGenoData(genObj, "marker", glines, plines)
    expect_true(all(grepl("^Chr\\.", colnames(out$genoData))))
})

test_that(".buildGenoData: only lines in plines returned", {
    genObj <- make_wgCross_interval(n_lines = 30, n_chr = 2, n_mar = 3)
    glines <- paste0("L", 1:30)
    plines <- paste0("L", 1:15)   # subset
    out    <- .buildGenoData(genObj, "marker", glines, plines)
    expect_equal(nrow(out$genoData), 15L)
})

test_that(".buildGenoData: state vector has length = ncol(genoData)", {
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 3, n_mar = 4)
    glines <- paste0("L", 1:20)
    plines <- glines
    out    <- .buildGenoData(genObj, "interval", glines, plines)
    expect_equal(length(out$state), ncol(out$genoData))
})

test_that(".buildGenoData: column names match mnams", {
    genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 3)
    glines <- paste0("L", 1:20)
    plines <- glines
    out    <- .buildGenoData(genObj, "marker", glines, plines)
    expect_equal(colnames(out$genoData), out$mnams)
})

# =============================================================================
# 6.  engine_genodata.R — .fixLines
# =============================================================================

.fixLines <- wgAim:::.fixLines

test_that(".fixLines: fix.lines=FALSE returns baseModel unchanged", {
    m     <- make_mock_asreml()
    pheno <- data.frame(id = factor(paste0("L", 1:10)),
                        yld = rnorm(10))
    geno  <- matrix(sample(c(-1, 1), 10 * 5, TRUE), 10, 5,
                    dimnames = list(paste0("L", 1:10), paste0("m", 1:5)))
    fl <- .fixLines(m, pheno, geno, "id", pheno$id, fix.lines = FALSE)
    # Should not have added Gomit column
    expect_false("Gomit" %in% names(fl$phenoData))
    expect_equal(fl$merge.by, "id")
})

test_that(".fixLines: fix.lines=FALSE preserves original rterms", {
    m     <- make_mock_asreml()
    pheno <- data.frame(id = factor(paste0("L", 1:10)), yld = rnorm(10))
    geno  <- matrix(sample(c(-1, 1), 10 * 5, TRUE), 10, 5,
                    dimnames = list(paste0("L", 1:10), paste0("m", 1:5)))
    fl <- .fixLines(m, pheno, geno, "id", pheno$id, fix.lines = FALSE)
    expect_equal(fl$genetic.term, "id")
})

test_that(".fixLines: fix.lines=TRUE with all lines matched skips update", {
    # When all pheno lines ARE in geno, no Gomit/Gsave is needed
    m     <- make_mock_asreml()
    ids   <- paste0("L", 1:10)
    pheno <- data.frame(id = factor(ids), yld = rnorm(10))
    geno  <- matrix(sample(c(-1, 1), 10 * 5, TRUE), 10, 5,
                    dimnames = list(ids, paste0("m", 1:5)))
    fl <- .fixLines(m, pheno, geno, "id", pheno$id, fix.lines = TRUE)
    expect_false("Gomit" %in% names(fl$phenoData))
})
