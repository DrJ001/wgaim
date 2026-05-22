# =============================================================================
# helper-fixtures.R
# Shared test fixtures for wgAim tests.
#
# All objects here are built without requiring ASReml-R, using synthetic data
# that mirrors the internal structures the package actually works with.
# ASReml-dependent tests use skip_if_not_installed("asreml").
#
# Loaded automatically by testthat before any test file runs.
# =============================================================================

# ---- 1. Minimal synthetic qtl cross objects ----------------------------------

#' Build a minimal 'bc' (backcross) qtl cross object for testing.
make_bc_cross <- function(n_lines = 40, n_chr = 3, n_mar = 5) {
    set.seed(42)
    chrs <- paste0("C", 1:n_chr)
    pos_list  <- lapply(seq_len(n_chr), function(i) seq(0, 80, length.out = n_mar))
    geno_list <- vector("list", n_chr)
    for (i in seq_len(n_chr)) {
        g <- matrix(sample(c(1L, 2L, NA_integer_),
                           n_lines * n_mar, replace = TRUE,
                           prob = c(0.45, 0.45, 0.1)),
                    nrow = n_lines)
        colnames(g) <- paste0(chrs[i], "_m", seq_len(n_mar))
        rownames(g) <- paste0("L", seq_len(n_lines))
        mp          <- pos_list[[i]]
        names(mp)   <- colnames(g)
        class(g)    <- "A"
        geno_list[[i]] <- list(data = g, map = mp)
    }
    names(geno_list) <- chrs
    pheno <- data.frame(
        id   = paste0("L", seq_len(n_lines)),
        yld  = rnorm(n_lines, 5, 1),
        stringsAsFactors = FALSE
    )
    obj <- list(geno = geno_list, pheno = pheno)
    class(obj) <- c("bc", "cross")
    obj
}

#' Build a minimal 'dh' cross (same structure as bc, just different class).
make_dh_cross <- function(n_lines = 30, n_chr = 2, n_mar = 4) {
    cr <- make_bc_cross(n_lines, n_chr, n_mar)
    class(cr) <- c("dh", "cross")
    cr
}

#' Build a minimal 'f2' cross (genotype codes 1/2/3 for AA/AB/BB).
make_f2_cross <- function(n_lines = 30, n_chr = 2, n_mar = 4) {
    cr <- make_bc_cross(n_lines, n_chr, n_mar)
    for (i in seq_along(cr$geno)) {
        cr$geno[[i]]$data[] <- sample(c(1L, 2L, 3L, NA_integer_),
                                       length(cr$geno[[i]]$data),
                                       replace = TRUE, prob = c(0.25, 0.5, 0.25, 0))
    }
    class(cr) <- c("f2", "cross")
    cr
}

# ---- 2. Minimal wgCross objects (post-primeCross) ----------------------------

#' Build a minimal wgCross (interval type) without calling primeCross().
#' Mirrors the exact structure that primeCross() produces.
make_wgCross_interval <- function(n_lines = 40, n_chr = 3, n_mar = 5,
                                   n_int = NULL) {
    set.seed(42)
    if (is.null(n_int)) n_int <- n_mar  # one interval per marker gap + edges
    chrs   <- paste0("C", 1:n_chr)
    ids    <- paste0("L", seq_len(n_lines))
    geno   <- vector("list", n_chr)

    for (i in seq_len(n_chr)) {
        # Map
        mp <- seq(0, 80, length.out = n_mar)
        names(mp) <- paste0(chrs[i], "_m", seq_len(n_mar))

        # Imputed marker data: lines x markers, coded -1/+1
        imp <- matrix(sample(c(-1, 1), n_lines * n_mar, replace = TRUE),
                      nrow = n_lines,
                      dimnames = list(ids, names(mp)))

        # Interval data: lines x intervals, coded -1/+1
        n_intervals <- n_mar  # n_mar intervals for n_mar markers
        int_nams <- paste0(chrs[i], "_int", seq_len(n_intervals))
        int_pos  <- seq(0, 80, length.out = n_intervals)
        names(int_pos) <- int_nams
        int_dat  <- matrix(sample(c(-1, 1), n_lines * n_intervals, replace = TRUE),
                           nrow = n_lines,
                           dimnames = list(ids, int_nams))

        geno[[i]] <- list(
            data         = imp,
            map          = mp,
            imputed.data = imp,
            interval.data = int_dat,
            inferred.map  = int_pos
        )
    }
    names(geno) <- chrs

    pheno <- data.frame(id = ids, yld = rnorm(n_lines, 5, 1),
                        stringsAsFactors = FALSE)
    obj <- list(geno = geno, pheno = pheno, type = "interval")
    class(obj) <- c("wgCross", "bc", "cross")
    obj
}

#' Build a minimal wgCross (marker type).
make_wgCross_marker <- function(n_lines = 40, n_chr = 3, n_mar = 5) {
    obj       <- make_wgCross_interval(n_lines, n_chr, n_mar)
    obj$type  <- "marker"
    # Remove interval-specific slots
    for (i in seq_along(obj$geno)) {
        obj$geno[[i]]$interval.data <- NULL
        obj$geno[[i]]$inferred.map  <- NULL
    }
    obj
}

# ---- 3. Minimal wgPanel object (post-primePanel) ----------------------------

#' Build a minimal wgPanel without calling primePanel().
make_wgPanel <- function(n_lines = 50, n_chr = 2, n_mar = 10) {
    set.seed(99)
    chrs <- paste0("P", 1:n_chr)
    ids  <- paste0("S", seq_len(n_lines))
    geno <- vector("list", n_chr)

    for (i in seq_len(n_chr)) {
        mp <- seq(0, 100, length.out = n_mar)
        names(mp) <- paste0(chrs[i], "_snp", seq_len(n_mar))
        imp <- matrix(sample(c(-1, 0, 1), n_lines * n_mar,
                             replace = TRUE, prob = c(0.25, 0.5, 0.25)),
                      nrow = n_lines,
                      dimnames = list(ids, names(mp)))
        geno[[i]] <- list(data = imp, map = mp, imputed.data = imp)
    }
    names(geno) <- chrs
    pheno <- data.frame(id = ids, yld = rnorm(n_lines, 10, 2),
                        stringsAsFactors = FALSE)
    obj <- list(geno = geno, pheno = pheno, type = "marker")
    class(obj) <- c("wgPanel", "wgCross", "cross")
    obj
}

# ---- 4. Mock ASReml-like model objects ---------------------------------------
# These are minimal lists that satisfy the checks wgAim functions do on an
# ASReml model object, without requiring ASReml to be installed.

#' Build a minimal mock ASReml model object.
#' @param converge Logical — whether $converge is TRUE.
#' @param loglik   Numeric — the log-likelihood.
#' @param data_name Character — name of the data object in a temp env.
#' @param data_val  Data frame — the phenoData to embed in the model call.
make_mock_asreml <- function(converge = TRUE,
                              loglik   = -50,
                              data_name = "phenoData",
                              data_val  = NULL,
                              merge.by  = "id") {
    if (is.null(data_val))
        data_val <- data.frame(id = paste0("L", 1:10), yld = rnorm(10),
                               stringsAsFactors = FALSE)
    data_val[[merge.by]] <- factor(data_val[[merge.by]])

    # Embed data in the function's local environment so eval() finds it
    e <- new.env(parent = emptyenv())
    assign(data_name, data_val, envir = e)

    call_obj         <- call("asreml",
                             fixed    = quote(yld ~ 1),
                             random   = as.formula(paste("~", merge.by)),
                             data     = as.name(data_name))

    # vparameters: one variance component + residual (constraint 0 = U, 4 = B)
    vpar     <- c(0.5, 1.0); names(vpar) <- c(merge.by, "R!variance")
    vpar.con <- c(0L, 0L);  names(vpar.con) <- names(vpar)
    sigma2   <- 1.0

    # Minimal coefficient list
    coef_rand <- matrix(rnorm(nrow(data_val)), ncol = 1,
                        dimnames = list(paste0(merge.by, "_", data_val[[merge.by]]),
                                        "effect"))
    coef_fix  <- matrix(c(5.0), ncol = 1, dimnames = list("(Intercept)", "effect"))

    formulae <- list(
        fixed  = as.formula(paste("yld ~ 1")),
        random = as.formula(paste("~", merge.by))
    )
    for (f in formulae) attr(f, ".Environment") <- e

    list(
        converge          = converge,
        loglik            = loglik,
        sigma2            = sigma2,
        vparameters       = vpar,
        vparameters.con   = vpar.con,
        vcoeff            = list(random = setNames(rep(0.1, nrow(data_val)),
                                                   rownames(coef_rand))),
        coefficients      = list(fixed = coef_fix, random = coef_rand),
        call              = call_obj,
        formulae          = formulae,
        mf                = data_val,
        converge          = converge
    )
}

# ---- 5. Fully-populated mock qtlAim object -----------------------------------
# Avoids needing ASReml to run the full pipeline for display/summary tests.

#' Build a mock qtlAim result object.
#' n_qtl: how many QTL to simulate in the $QTL slot.
make_mock_qtlAim <- function(n_qtl = 2,
                              gen.type  = "interval",
                              method    = "fixed",
                              selection = "interval",
                              n_lines   = 40,
                              n_chr     = 3,
                              n_mar     = 5) {
    set.seed(7)
    genObj  <- make_wgCross_interval(n_lines, n_chr, n_mar)
    chrs    <- names(genObj$geno)
    n_int   <- n_mar  # intervals per chr

    # Pick n_qtl distinct (chr, int-index) positions
    chr_idx <- chrs[seq_len(n_qtl) %% n_chr + 1]
    int_idx <- seq_len(n_qtl) + 1L

    qtl_keys  <- paste("Chr", chr_idx, int_idx, sep = ".")
    eff_names <- paste("X",   chr_idx, int_idx, sep = ".")
    effects   <- setNames(rnorm(n_qtl, 0.4, 0.1), eff_names)
    veffects  <- setNames(runif(n_qtl, 0.02, 0.05), eff_names)

    # Build per-iteration effect/veffect lists
    coef.list  <- lapply(seq_len(n_qtl), function(k) effects[seq_len(k)])
    vcoef.list <- lapply(seq_len(n_qtl), function(k) veffects[seq_len(k)])

    # Minimal diag
    n_ints_total <- n_chr * n_mar
    state <- rep(1, n_ints_total)
    all_keys <- paste("Chr", rep(chrs, each = n_mar),
                      rep(seq_len(n_mar), times = n_chr), sep = ".")
    names(state) <- all_keys

    oint_list  <- lapply(seq_len(n_qtl), function(k) {
        v <- runif(n_ints_total, 0, 2); names(v) <- all_keys; v
    })
    blup_list  <- lapply(seq_len(n_qtl), function(k) {
        v <- rnorm(n_ints_total);       names(v) <- all_keys; v
    })
    lik_list   <- lapply(seq_len(n_qtl), function(k)
        list(baseLogL = -50 + k, stat = 4 + k, pvalue = 0.02, pass = TRUE))
    ochr_list  <- lapply(seq_len(n_qtl), function(k) {
        v <- runif(n_chr, 0, 5); names(v) <- chrs; v
    })

    lik.mat <- matrix(
        c(sapply(lik_list, function(l) c(l$baseLogL, l$baseLogL + l$stat/2, l$stat, l$pvalue))),
        ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue"))
    )

    QTL <- list(
        qtl        = qtl_keys,
        effects    = effects,
        veffects   = veffects,
        method     = method,
        type       = gen.type,
        selection  = selection,
        TypeI      = 0.05,
        iterations = n_qtl + 1L,
        breakout   = FALSE,
        diag       = list(
            oint       = oint_list,
            blups      = blup_list,
            lik        = lik_list,
            ochr       = ochr_list,
            lik.mat    = lik.mat,
            state      = state,
            genetic.term = "id",
            rel.scale  = 1,
            coef.list  = coef.list,
            vcoef.list = vcoef.list
        )
    )

    # Attach QTL covariate columns to a phenoData (needed for contrast plots)
    ids      <- paste0("L", seq_len(n_lines))
    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines, 5, 1),
                             stringsAsFactors = FALSE)
    for (k in seq_len(n_qtl)) {
        xcol <- eff_names[k]
        chr_k <- chr_idx[k]; int_k <- int_idx[k]
        phenoData[[xcol]] <- genObj$geno[[chr_k]]$interval.data[ids, int_k]
    }

    # Build minimal mock ASReml body
    vpar     <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
    vpar.con <- c(0L, 0L); names(vpar.con) <- names(vpar)
    formulae <- list(fixed  = yld ~ 1, random = ~ id)
    mf       <- phenoData

    obj <- list(
        converge        = TRUE,
        loglik          = -45,
        sigma2          = 1.0,
        vparameters     = vpar,
        vparameters.con = vpar.con,
        coefficients    = list(fixed  = matrix(5, 1, 1,
                                               dimnames = list("(Intercept)", "effect")),
                               random = matrix(rnorm(n_lines), n_lines, 1,
                                               dimnames = list(paste0("id_", ids), "effect"))),
        formulae        = formulae,
        mf              = mf,
        call            = call("asreml", fixed = quote(yld ~ 1),
                               random = quote(~ id), data = quote(phenoData)),
        QTL             = QTL
    )
    class(obj) <- c("qtlAim", "asreml")
    attr(obj, "genObj") <- genObj   # stash for convenience
    obj
}

#' Build a mock gwasAim result object.
make_mock_gwasAim <- function(n_qtl = 2,
                               n_lines = 50,
                               n_chr   = 2,
                               n_mar   = 10) {
    set.seed(13)
    genObj <- make_wgPanel(n_lines, n_chr, n_mar)
    chrs   <- names(genObj$geno)

    chr_idx  <- chrs[seq_len(n_qtl) %% n_chr + 1]
    int_idx  <- seq_len(n_qtl) + 2L
    qtl_keys <- paste("Chr", chr_idx, int_idx, sep = ".")
    eff_names <- paste("X", chr_idx, int_idx, sep = ".")
    effects  <- setNames(rnorm(n_qtl, 0.3, 0.08), eff_names)
    veffects <- setNames(runif(n_qtl, 0.01, 0.03), eff_names)

    coef.list  <- lapply(seq_len(n_qtl), function(k) effects[seq_len(k)])
    vcoef.list <- lapply(seq_len(n_qtl), function(k) veffects[seq_len(k)])

    n_total <- n_chr * n_mar
    all_keys <- paste("Chr", rep(chrs, each = n_mar),
                      rep(seq_len(n_mar), times = n_chr), sep = ".")
    state <- rep(1, n_total); names(state) <- all_keys

    oint_list <- lapply(seq_len(n_qtl), function(k) {
        v <- runif(n_total, 0, 2); names(v) <- all_keys; v
    })
    blup_list <- lapply(seq_len(n_qtl), function(k) {
        v <- rnorm(n_total); names(v) <- all_keys; v
    })
    lik_list  <- lapply(seq_len(n_qtl), function(k)
        list(baseLogL = -60 + k, stat = 5 + k, pvalue = 0.01, pass = TRUE))
    lik.mat <- matrix(
        c(sapply(lik_list, function(l) c(l$baseLogL, l$baseLogL + l$stat/2, l$stat, l$pvalue))),
        ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue"))
    )

    QTL <- list(
        qtl        = qtl_keys,
        effects    = effects,
        veffects   = veffects,
        method     = "fixed",
        type       = "marker",
        selection  = "interval",
        TypeI      = 0.05,
        n.markers  = n_total,
        iterations = n_qtl + 1L,
        breakout   = FALSE,
        diag       = list(
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

    ids       <- paste0("S", seq_len(n_lines))
    phenoData <- data.frame(id = factor(ids), yld = rnorm(n_lines, 10, 2),
                             stringsAsFactors = FALSE)
    for (k in seq_len(n_qtl)) {
        xcol  <- eff_names[k]
        chr_k <- chr_idx[k]; int_k <- int_idx[k]
        phenoData[[xcol]] <- genObj$geno[[chr_k]]$imputed.data[ids, int_k]
    }

    vpar     <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
    vpar.con <- c(0L, 0L); names(vpar.con) <- names(vpar)

    obj <- list(
        converge        = TRUE,
        loglik          = -55,
        sigma2          = 1.0,
        vparameters     = vpar,
        vparameters.con = vpar.con,
        coefficients    = list(fixed  = matrix(10, 1, 1,
                                               dimnames = list("(Intercept)", "effect")),
                               random = matrix(rnorm(n_lines), n_lines, 1,
                                               dimnames = list(paste0("id_", ids), "effect"))),
        formulae        = list(fixed = yld ~ 1, random = ~ id),
        mf              = phenoData,
        call            = call("asreml", fixed = quote(yld ~ 1),
                               random = quote(~ id), data = quote(phenoData)),
        QTL             = QTL
    )
    class(obj) <- c("gwasAim", "asreml")
    attr(obj, "genObj") <- genObj
    obj
}

#' Build a mock gpAim result object.
make_mock_gpAim <- function(n_lines = 40, n_chr = 2, n_mar = 8,
                             path = "vm") {
    set.seed(17)
    genObj <- make_wgCross_interval(n_lines, n_chr, n_mar)
    ids    <- paste0("L", seq_len(n_lines))

    SE_vals <- runif(n_lines, 0.1, 0.3)
    Vg_mock <- 0.6
    gebv_df <- data.frame(
        id       = ids,
        GEBV     = rnorm(n_lines, 0, 1),
        SE       = SE_vals,
        Accuracy = sqrt(pmax(0, 1 - SE_vals^2 / Vg_mock)),
        gen.H2   = 0.75,
        stringsAsFactors = FALSE
    )

    n_markers <- n_chr * n_mar
    raw_mat   <- matrix(rnorm(n_lines * n_markers), nrow = n_lines)
    rel_mat   <- tcrossprod(raw_mat) / n_markers
    rownames(rel_mat) <- colnames(rel_mat) <- ids

    GP <- list(
        gebv         = gebv_df,
        marker.effects = setNames(rnorm(n_markers), paste0("m", seq_len(n_markers))),
        gen.type     = "interval",
        path         = path,
        var.genetic  = 0.6,
        var.resid    = 0.4,
        heritability = 0.6,
        gen.H2       = 0.75,
        n.markers    = n_markers,
        rel.scale    = n_markers,
        rel.matrix   = rel_mat,
        genetic.term = "id"
    )

    vpar     <- c(0.6, 1.0); names(vpar) <- c("id", "R!variance")
    vpar.con <- c(0L, 0L); names(vpar.con) <- names(vpar)

    obj <- list(
        converge        = TRUE,
        loglik          = -80,
        sigma2          = 1.0,
        vparameters     = vpar,
        vparameters.con = vpar.con,
        coefficients    = list(fixed  = matrix(5, 1, 1,
                                               dimnames = list("(Intercept)", "effect")),
                               random = matrix(rnorm(n_lines), n_lines, 1,
                                               dimnames = list(paste0("id_", ids), "effect"))),
        formulae        = list(fixed = yld ~ 1, random = ~ id),
        mf              = data.frame(id = factor(ids), yld = rnorm(n_lines)),
        call            = call("asreml"),
        GP              = GP
    )
    class(obj) <- c("gpAim", "asreml")
    attr(obj, "genObj") <- genObj
    obj
}

# ---- 6. Convenience: build a fineMap object ----------------------------------

#' Build a minimal fineMap object without calling fineMap().
make_mock_fineMap <- function(n_qtl = 2, n_pos = 15) {
    set.seed(3)
    qtl_keys <- paste0("Chr.C", seq_len(n_qtl), ".2")
    result <- lapply(seq_len(n_qtl), function(k) {
        dists <- seq(0, 50, length.out = n_pos)
        lod   <- dnorm(dists, mean = 25, sd = 8) * 20 + rnorm(n_pos, 0, 0.2)
        df <- data.frame(mark   = paste0("m", seq_len(n_pos)),
                         dist   = dists,
                         pvalue = pchisq(lod * 4.6, 1, lower.tail = FALSE),
                         LOD    = lod,
                         stringsAsFactors = FALSE)
        attr(df, "qtl_pos") <- 25
        attr(df, "clipped") <- FALSE
        df
    })
    names(result) <- qtl_keys
    attr(result, "qtl.key") <- qtl_keys
    attr(result, "type")    <- "interval"
    class(result) <- "fineMap"
    result
}

# ---- 7. Shared minimal post-update ASReml mock (used by engine-asreml and
#         fineMap-mv tests) ----------------------------------------------------

#' Build a minimal post-update model that satisfies all engine reads.
#' n_int: total number of intervals/markers across all chromosomes.
make_updated_model <- function(merge.by = "id", n_int = 15, loglik = -50,
                                lrt_pass = TRUE, sigma2 = 1.0) {
    mbf_name <- paste0("mbf(ints)_", merge.by, "!mbf(ints)!var")
    vpar     <- c(0.3, 1.0); names(vpar) <- c(mbf_name, "R!variance")
    vpar.con <- c(0L, 0L);   names(vpar.con) <- names(vpar)
    ll       <- if (lrt_pass) loglik + 2 else loglik - 2
    int_nams <- paste0("mbf(ints)_", merge.by, "_", seq_len(n_int))
    rand_mat <- matrix(rnorm(n_int, 0, 0.1), n_int, 1,
                       dimnames = list(int_nams, "effect"))
    fix_mat  <- matrix(5.0, 1, 1, dimnames = list("(Intercept)", "effect"))
    list(
        converge        = TRUE,
        loglik          = ll,
        sigma2          = sigma2,
        vparameters     = vpar,
        vparameters.con = vpar.con,
        vcoeff          = list(
            random = setNames(rep(0.05, n_int), int_nams),
            fixed  = setNames(0.1, "(Intercept)")
        ),
        coefficients    = list(fixed = fix_mat, random = rand_mat),
        call            = call("asreml"),
        G.param         = list(),
        formulae        = list(fixed = yld ~ 1, random = ~ id),
        mf              = data.frame(id = factor(paste0("L", 1:5)), yld = 1:5)
    )
}

# ---- 8. MV mock qtlAim object (Trait path) -----------------------------------

#' Build a mock multivariate qtlAim object for fineMap MV tests.
#'
#' @param n_qtl       Number of QTL to simulate.
#' @param trait_name  Name of the Trait factor column (default "Trial").
#' @param n_trials    Number of trait/trial levels.
#' @param interact    Logical vector (length n_qtl): TRUE = QTL retained as
#'                    Trait:X.chr.idx interaction; FALSE = main effect only.
#'                    Recycled if shorter than n_qtl.
#' @param n_lines     Number of lines.
#' @param n_chr       Number of chromosomes.
#' @param n_mar       Markers per chromosome.
make_mock_mv_qtlAim <- function(n_qtl     = 2,
                                 trait_name = "Trial",
                                 n_trials   = 3,
                                 interact   = c(TRUE, FALSE),
                                 n_lines    = 40,
                                 n_chr      = 3,
                                 n_mar      = 5) {
    set.seed(77)
    interact <- rep_len(interact, n_qtl)
    genObj   <- make_wgCross_interval(n_lines, n_chr, n_mar)
    chrs     <- names(genObj$geno)
    trials   <- paste0("T", seq_len(n_trials))

    chr_idx  <- chrs[seq_len(n_qtl) %% n_chr + 1]
    int_idx  <- seq_len(n_qtl) + 1L
    qtl_keys  <- paste("Chr", chr_idx, int_idx, sep = ".")
    eff_names <- paste("X",   chr_idx, int_idx, sep = ".")

    # Fixed formula terms: interaction QTL -> "Trial:X.chr.idx"; main -> "X.chr.idx"
    final_terms <- ifelse(interact,
                          paste0(trait_name, ":", eff_names),
                          eff_names)

    # Fixed coefficient matrix: intercepts-by-trial + one row per final term
    trial_int_nams <- paste0(trait_name, "_", trials)
    fix_nams <- c("(Intercept)", trial_int_nams, final_terms)
    fix_mat  <- matrix(rnorm(length(fix_nams), 0.3, 0.1), length(fix_nams), 1,
                       dimnames = list(fix_nams, "effect"))

    n_ints_total <- n_chr * n_mar
    all_keys     <- paste("Chr", rep(chrs, each = n_mar),
                          rep(seq_len(n_mar), times = n_chr), sep = ".")
    state        <- setNames(rep(1L, n_ints_total), all_keys)

    effects  <- setNames(rnorm(n_qtl, 0.4, 0.1), eff_names)
    veffects <- setNames(runif(n_qtl, 0.02, 0.05), eff_names)

    coef.list  <- lapply(seq_len(n_qtl), function(k) effects[seq_len(k)])
    vcoef.list <- lapply(seq_len(n_qtl), function(k) veffects[seq_len(k)])

    lik_list <- lapply(seq_len(n_qtl), function(k)
        list(baseLogL = -50 + k, stat = 4 + k, pvalue = 0.02, pass = TRUE))
    lik.mat  <- matrix(
        c(sapply(lik_list, function(l) c(l$baseLogL, l$baseLogL + l$stat/2, l$stat, l$pvalue))),
        ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue"))
    )

    # Wald test result placeholder
    wald_df <- data.frame(
        "Wald Statistic" = round(runif(n_qtl, 5, 15), 4),
        "P-Value"        = round(runif(n_qtl, 0, 0.03), 4),
        check.names  = FALSE,
        row.names    = final_terms
    )

    QTL <- list(
        qtl           = qtl_keys,
        effects       = effects,
        veffects      = veffects,
        method        = "fixed",
        type          = "interval",
        selection     = "interval",
        TypeI         = 0.05,
        iterations    = n_qtl + 1L,
        breakout      = FALSE,
        Trait         = trait_name,
        trait.levels  = trials,
        wald.test     = wald_df,
        final.terms   = final_terms,
        is.interaction = interact,
        diag          = list(
            oint         = lapply(seq_len(n_qtl), function(k) {
                v <- runif(n_ints_total, 0, 2); names(v) <- all_keys; v }),
            blups        = lapply(seq_len(n_qtl), function(k) {
                m <- matrix(rnorm(n_ints_total * n_trials), n_ints_total, n_trials,
                            dimnames = list(all_keys, trials)); m }),
            lik          = lik_list,
            lik.mat      = lik.mat,
            state        = state,
            genetic.term = "id",
            rel.scale    = 1,
            coef.list    = coef.list,
            vcoef.list   = vcoef.list
        )
    )

    ids       <- paste0("L", seq_len(n_lines))
    phenoData <- data.frame(
        id    = factor(rep(ids, each = n_trials)),
        Trial = factor(rep(trials, times = n_lines)),
        yld   = rnorm(n_lines * n_trials, 5, 1),
        stringsAsFactors = FALSE
    )
    for (k in seq_len(n_qtl)) {
        xcol <- eff_names[k]
        phenoData[[xcol]] <- genObj$geno[[chr_idx[k]]]$interval.data[
            as.character(phenoData$id), int_idx[k]
        ]
    }

    vpar     <- setNames(c(0.3, 1.0), c("id", "R!variance"))
    vpar.con <- setNames(c(0L,  0L),  names(vpar))

    # Fixed formula contains the pruned terms
    fix_formula <- as.formula(
        paste("yld ~", trait_name, "+",
              paste(final_terms, collapse = " + "))
    )

    obj <- list(
        converge        = TRUE,
        loglik          = -45,
        sigma2          = 1.0,
        vparameters     = vpar,
        vparameters.con = vpar.con,
        coefficients    = list(fixed  = fix_mat,
                               random = matrix(rnorm(n_lines), n_lines, 1,
                                               dimnames = list(paste0("id_", ids),
                                                               "effect"))),
        formulae        = list(fixed = fix_formula, random = ~ id),
        mf              = phenoData,
        call            = call("asreml",
                               fixed  = fix_formula,
                               random = quote(~ id),
                               data   = quote(phenoData)),
        QTL             = QTL
    )
    environment(obj$call$fixed)    <- globalenv()
    environment(obj$formulae$fixed) <- globalenv()
    class(obj) <- c("qtlAim", "asreml")
    attr(obj, "genObj") <- genObj
    obj
}

# ---- 9. Raw panel data for checkPanel / filterPanel tests -------------------

#' Build a synthetic raw genotype matrix + map with controllable issues.
#'
#' @param n_lines      Number of lines.
#' @param n_chr        Number of chromosomes.
#' @param n_mar_chr    Markers per chromosome.
#' @param miss_rate    Baseline random missing rate.
#' @param n_dup_lines  Number of extra duplicate lines to inject.
#' @param n_dup_marks  Number of extra duplicate markers to inject.
#' @param n_mono       Number of monomorphic markers to inject.
#' @param n_highmiss_marks  Number of markers forced to high missingness.
#' @param n_highmiss_lines  Number of lines forced to high missingness.
#' @param n_highhet    Number of lines forced to high heterozygosity.
#' @param n_notinmap   Number of geno markers to exclude from the map.
#' @param encoding     `"012"` (default) or `"pm1"`.
make_raw_panel <- function(n_lines          = 60,
                            n_chr            = 2,
                            n_mar_chr        = 20,
                            miss_rate        = 0.03,
                            n_dup_lines      = 2,
                            n_dup_marks      = 2,
                            n_mono           = 3,
                            n_highmiss_marks = 2,
                            n_highmiss_lines = 2,
                            n_highhet        = 2,
                            n_notinmap       = 1,
                            encoding         = "012") {
    set.seed(101)
    n_markers <- n_chr * n_mar_chr
    chrs      <- paste0("Chr", seq_len(n_chr))
    line_ids  <- paste0("L", sprintf("%02d", seq_len(n_lines)))
    mar_ids   <- paste0(rep(chrs, each = n_mar_chr), "_M",
                        sprintf("%02d", rep(seq_len(n_mar_chr), times = n_chr)))

    # ---- Index regions so injected issues never overlap ----
    # Marker layout (all indices into 1..n_markers):
    #   [1..n_mono]                         : monomorphic (always 0)
    #   [n_mono+1 .. n_mono+n_highmiss_marks]: high-missingness markers
    #   [n_mono+n_highmiss_marks+1 .. -dup-notinmap]: polymorphic base markers
    #   [n_markers-n_notinmap-n_dup_marks+1 ..
    #    n_markers-n_notinmap]              : duplicate markers (in map)
    #   [n_markers-n_notinmap+1..n_markers] : not-in-map markers
    #
    # Line layout (1..n_lines):
    #   [1..n_highmiss_lines]               : high-missingness lines
    #   [n_highmiss_lines+1..+n_highhet]    : high-heterozygosity lines
    #   [n_lines-n_dup_lines+1..n_lines]    : duplicate lines
    #   source line for duplicates          : n_highmiss_lines+n_highhet+1

    src_line <- n_highmiss_lines + n_highhet + 1L   # the line we copy from
    src_mark <- n_mono + n_highmiss_marks + 1L      # the marker we copy from

    # Allele frequencies: polymorphic markers get realistic frequencies
    p_freq <- runif(n_markers, 0.2, 0.8)
    p_freq[seq_len(n_mono)] <- 0   # force first n_mono to monomorphic

    # Generate 0/1/2 base genotypes
    geno <- matrix(NA_real_, nrow = n_lines, ncol = n_markers,
                   dimnames = list(line_ids, mar_ids))
    for (j in seq_len(n_markers)) {
        pj <- p_freq[j]
        if (pj == 0) {
            geno[, j] <- 0
        } else {
            pr <- c((1 - pj)^2, 2 * pj * (1 - pj), pj^2)
            geno[, j] <- sample(0:2, n_lines, replace = TRUE, prob = pr)
        }
    }

    # Random background missing (only on polymorphic non-highhet lines)
    safe_lines <- seq(n_highmiss_lines + 1L,
                      n_lines - n_dup_lines)
    safe_marks <- seq(n_mono + 1L, n_markers - n_notinmap - n_dup_marks)
    if (miss_rate > 0 && length(safe_lines) && length(safe_marks)) {
        sub  <- geno[safe_lines, safe_marks, drop = FALSE]
        idx  <- sample(length(sub), round(miss_rate * length(sub)))
        sub[idx] <- NA
        geno[safe_lines, safe_marks] <- sub
    }

    # High-missingness markers (avoid mono markers)
    if (n_highmiss_marks > 0) {
        hi_m <- seq(n_mono + 1L, n_mono + n_highmiss_marks)
        for (j in hi_m)
            geno[sample(n_lines, round(0.35 * n_lines)), j] <- NA
    }

    # High-missingness lines (only on polymorphic markers to preserve mono=0)
    poly_marks <- seq(n_mono + 1L, n_markers)
    if (n_highmiss_lines > 0) {
        for (i in seq_len(n_highmiss_lines))
            geno[i, sample(poly_marks, round(0.35 * length(poly_marks)))] <- NA
    }

    # High-heterozygosity lines (only on polymorphic markers)
    if (n_highhet > 0) {
        het_idx <- seq(n_highmiss_lines + 1L, n_highmiss_lines + n_highhet)
        for (i in het_idx)
            geno[i, sample(poly_marks, round(0.60 * length(poly_marks)))] <- 1L
    }

    # Duplicate lines: copy src_line to last n_dup_lines positions
    if (n_dup_lines > 0 && src_line <= n_lines - n_dup_lines) {
        for (i in seq_len(n_dup_lines))
            geno[n_lines - n_dup_lines + i, ] <- geno[src_line, ]
    }

    # Duplicate markers: copy src_mark to positions just before the
    # not-in-map zone so they ARE in the map (and thus in common with geno)
    dup_mark_positions <- seq(n_markers - n_notinmap - n_dup_marks + 1L,
                               n_markers - n_notinmap)
    if (n_dup_marks > 0 && length(dup_mark_positions) > 0L) {
        for (j in dup_mark_positions)
            geno[, j] <- geno[, src_mark]
    }

    # Convert to pm1 if requested (after all injections in 012 space)
    if (encoding == "pm1")
        geno <- geno - 1

    # Build map: exclude the last n_notinmap markers (never duplicates)
    map_df <- data.frame(
        marker = mar_ids,
        chr    = rep(chrs, each = n_mar_chr),
        pos    = rep(seq(0, 100, length.out = n_mar_chr), times = n_chr),
        stringsAsFactors = FALSE
    )
    if (n_notinmap > 0)
        map_df <- map_df[seq_len(nrow(map_df) - n_notinmap), , drop = FALSE]

    list(geno = geno, map = map_df, encoding = encoding,
         n_lines = n_lines, n_markers = n_markers, n_chr = n_chr,
         n_mar_chr = n_mar_chr, n_mono = n_mono,
         n_dup_lines = n_dup_lines, n_dup_marks = n_dup_marks,
         n_highmiss_marks = n_highmiss_marks,
         n_highmiss_lines = n_highmiss_lines,
         n_highhet = n_highhet, n_notinmap = n_notinmap)
}
