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

    gebv_df <- data.frame(
        id   = ids,
        GEBV = rnorm(n_lines, 0, 1),
        SE   = runif(n_lines, 0.1, 0.3),
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
