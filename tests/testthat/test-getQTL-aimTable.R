# =============================================================================
# test-getQTL-aimTable.R
# Tests for getQTL(), aimTable(), and print.aimTable().
#
# Fixtures from helper-fixtures.R (loaded automatically by testthat):
#   make_wgCross_interval(), make_wgCross_marker(), make_wgPanel()
#
# Internal helper builders in this file produce self-consistent mock
# qtlAim / gwasAim objects whose naming conventions match what the
# real wgAim engine produces, so that getQTL() / summary.*() work without
# running ASReml-R.
#
# Key naming convention (engine_genodata.R):
#   Column keys  : paste("Chr", chr, int_index, sep = ".")  -> "Chr.C1.1"
#   Effect names : sub("^Chr", "X", key)                    -> "X.C1.1"
#   getQTL parses: substring(name, 3) -> "C1.1"; split "." -> chr="C1", i=1
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: build a summary-compatible mock qtlAim object.
#
# Keys follow the engine convention "Chr.<chrname>.<index>".
# Effect names follow "X.<chrname>.<index>".
# ---------------------------------------------------------------------------
build_summable_qtlAim <- function(n_qtl   = 2,
                                   n_lines = 40,
                                   n_chr   = 3,
                                   n_mar   = 5) {
  set.seed(77)
  genObj <- make_wgCross_interval(n_lines, n_chr, n_mar)
  chrs   <- names(genObj$geno)
  ids    <- paste0("L", seq_len(n_lines))

  # Build engine-style keys: "Chr.C1.1", "Chr.C1.2", ...
  gdat     <- lapply(genObj$geno, function(el) el$interval.data)
  nint_seq <- lapply(gdat, function(el) seq_len(ncol(el)))
  lint     <- unlist(lapply(nint_seq, length))
  all_keys <- paste("Chr",
                    rep(names(genObj$geno), times = lint),
                    unlist(nint_seq),
                    sep = ".")
  n_total  <- length(all_keys)

  # Mark the first n_qtl keys as detected QTL
  state    <- rep(0L, n_total); names(state) <- all_keys
  qtl_keys <- all_keys[seq_len(n_qtl)]
  for (k in seq_len(n_qtl)) state[qtl_keys[k]] <- 1L

  # Effect names: "X.C1.1", "X.C1.2", ...
  eff_names <- sub("^Chr", "X", qtl_keys)
  effects   <- setNames(rnorm(n_qtl, 0.4, 0.05), eff_names)
  veffects  <- setNames(runif(n_qtl, 0.02, 0.04), eff_names)

  oint_list  <- lapply(seq_len(n_qtl), function(k) {
    v <- runif(n_total, 0, 2); names(v) <- all_keys; v })
  blup_list  <- lapply(seq_len(n_qtl), function(k) {
    v <- rnorm(n_total); names(v) <- all_keys; v })
  lik_list   <- lapply(seq_len(n_qtl), function(k)
    list(baseLogL = -50 + k, stat = 4 + k, pvalue = 0.02, pass = TRUE))
  lik.mat    <- matrix(
    c(sapply(lik_list, function(l)
        c(l$baseLogL, l$baseLogL + l$stat / 2, l$stat, l$pvalue))),
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue"))
  )
  coef.list  <- lapply(seq_len(n_qtl), function(k) effects[seq_len(k)])
  vcoef.list <- lapply(seq_len(n_qtl), function(k) veffects[seq_len(k)])

  QTL <- list(
    qtl        = qtl_keys,
    effects    = effects,
    veffects   = veffects,
    method     = "fixed",
    type       = "interval",
    selection  = "interval",
    TypeI      = 0.05,
    iterations = n_qtl + 1L,
    breakout   = FALSE,
    diag       = list(
      oint         = oint_list,
      blups        = blup_list,
      lik          = lik_list,
      ochr         = NULL,
      lik.mat      = lik.mat,
      state        = state,
      genetic.term = "id",
      rel.scale    = 1L,
      coef.list    = coef.list,
      vcoef.list   = vcoef.list
    )
  )

  mf <- data.frame(id = factor(ids), yld = rnorm(n_lines, 5, 1),
                   stringsAsFactors = FALSE)
  # Attach covariate columns (summary.qtlAim builds genoSub from these levels)
  for (k in seq_len(n_qtl)) {
    parts  <- strsplit(qtl_keys[k], "\\.")[[1]]
    chr_k  <- parts[2]; int_k <- as.integer(parts[3])
    xcol   <- eff_names[k]
    mf[[xcol]] <- genObj$geno[[chr_k]]$interval.data[ids, int_k]
  }

  vpar     <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
  vpar.con <- c(0L, 0L);  names(vpar.con) <- names(vpar)

  obj <- list(
    converge        = TRUE,
    loglik          = -45,
    sigma2          = 1.0,
    vparameters     = vpar,
    vparameters.con = vpar.con,
    coefficients    = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf       = mf,
    call     = call("asreml", fixed = quote(yld ~ 1),
                    random = quote(~ id), data = quote(mf)),
    QTL      = QTL
  )
  class(obj) <- c("qtlAim", "asreml")
  list(obj = obj, genObj = genObj)
}

# ---------------------------------------------------------------------------
# Helper: build a summary-compatible mock gwasAim object (marker mode).
# Uses wgPanel from make_wgPanel(); key naming follows the engine convention.
# ---------------------------------------------------------------------------
build_summable_gwasAim <- function(n_qtl   = 2,
                                    n_lines = 50,
                                    n_chr   = 2,
                                    n_mar   = 10) {
  set.seed(88)
  genObj <- make_wgPanel(n_lines, n_chr, n_mar)
  chrs   <- names(genObj$geno)
  ids    <- paste0("S", seq_len(n_lines))

  # Engine naming for marker mode uses imputed.data columns
  gdat     <- lapply(genObj$geno, function(el) el$imputed.data)
  nint_seq <- lapply(gdat, function(el) seq_len(ncol(el)))
  lint     <- unlist(lapply(nint_seq, length))
  all_keys <- paste("Chr",
                    rep(names(genObj$geno), times = lint),
                    unlist(nint_seq),
                    sep = ".")
  n_total  <- length(all_keys)

  state    <- rep(0L, n_total); names(state) <- all_keys
  qtl_keys <- all_keys[seq_len(n_qtl)]
  for (k in seq_len(n_qtl)) state[qtl_keys[k]] <- 1L

  eff_names <- sub("^Chr", "X", qtl_keys)
  effects   <- setNames(rnorm(n_qtl, 0.3, 0.05), eff_names)
  veffects  <- setNames(runif(n_qtl, 0.01, 0.03), eff_names)

  oint_list <- lapply(seq_len(n_qtl), function(k) {
    v <- runif(n_total, 0, 2); names(v) <- all_keys; v })
  blup_list <- lapply(seq_len(n_qtl), function(k) {
    v <- rnorm(n_total); names(v) <- all_keys; v })
  lik_list  <- lapply(seq_len(n_qtl), function(k)
    list(baseLogL = -60 + k, stat = 5 + k, pvalue = 0.01, pass = TRUE))
  lik.mat   <- matrix(
    c(sapply(lik_list, function(l)
        c(l$baseLogL, l$baseLogL + l$stat / 2, l$stat, l$pvalue))),
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, c("L0", "L1", "Statistic", "Pvalue"))
  )
  coef.list  <- lapply(seq_len(n_qtl), function(k) effects[seq_len(k)])
  vcoef.list <- lapply(seq_len(n_qtl), function(k) veffects[seq_len(k)])

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
      oint         = oint_list,
      blups        = blup_list,
      lik          = lik_list,
      ochr         = NULL,
      lik.mat      = lik.mat,
      state        = state,
      genetic.term = "id",
      rel.scale    = 1L,
      coef.list    = coef.list,
      vcoef.list   = vcoef.list
    )
  )

  mf <- data.frame(id = factor(ids), yld = rnorm(n_lines, 10, 2),
                   stringsAsFactors = FALSE)
  for (k in seq_len(n_qtl)) {
    parts <- strsplit(qtl_keys[k], "\\.")[[1]]
    chr_k <- parts[2]; int_k <- as.integer(parts[3])
    xcol  <- eff_names[k]
    mf[[xcol]] <- genObj$geno[[chr_k]]$imputed.data[ids, int_k]
  }

  vpar     <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
  vpar.con <- c(0L, 0L);  names(vpar.con) <- names(vpar)

  obj <- list(
    converge        = TRUE,
    loglik          = -55,
    sigma2          = 1.0,
    vparameters     = vpar,
    vparameters.con = vpar.con,
    coefficients    = list(
      fixed  = matrix(10, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf       = mf,
    call     = call("asreml", fixed = quote(yld ~ 1),
                    random = quote(~ id), data = quote(mf)),
    QTL      = QTL
  )
  class(obj) <- c("gwasAim", "asreml")
  list(obj = obj, genObj = genObj)
}

# =============================================================================
# getQTL() — interval type
# =============================================================================

test_that("getQTL on interval-type qtlAim returns an 8-column character matrix", {
  fix  <- build_summable_qtlAim(n_qtl = 2)
  qtlm <- getQTL(fix$obj, fix$genObj)

  expect_true(is.matrix(qtlm))
  expect_equal(ncol(qtlm), 8L)
  expect_equal(nrow(qtlm), 2L)
  expect_true(is.character(qtlm))
})

test_that("getQTL interval-type: column 1 is valid chromosome, cols 3/5/7 are non-empty strings", {
  fix  <- build_summable_qtlAim(n_qtl = 2)
  qtlm <- getQTL(fix$obj, fix$genObj)

  chrs <- names(fix$genObj$geno)
  # Column 1: chromosome names should all be in the genObj chromosome list
  expect_true(all(qtlm[, 1] %in% chrs))
  # Columns 3 (left marker), 5 (inferred marker), 7 (right marker): non-empty strings
  for (col in c(3L, 5L, 7L))
    expect_true(all(nchar(qtlm[, col]) > 0),
                info = paste("Column", col, "should contain non-empty marker names"))
})

test_that("getQTL interval-type: columns 4, 6, 8 are numeric-like cM strings", {
  fix  <- build_summable_qtlAim(n_qtl = 2)
  qtlm <- getQTL(fix$obj, fix$genObj)

  for (col in c(4L, 6L, 8L)) {
    vals <- suppressWarnings(as.numeric(qtlm[, col]))
    expect_false(anyNA(vals),
                 info = paste("Column", col, "should be numeric-like cM values"))
  }
})

# =============================================================================
# getQTL() — marker type
# =============================================================================

test_that("getQTL on marker-type qtlAim returns a 4-column character matrix", {
  set.seed(55)
  genObj <- make_wgCross_marker(n_lines = 30, n_chr = 2, n_mar = 5)
  chrs   <- names(genObj$geno)
  ids    <- paste0("L", seq_len(30))

  gdat     <- lapply(genObj$geno, function(el) el$imputed.data)
  nint_seq <- lapply(gdat, function(el) seq_len(ncol(el)))
  lint     <- unlist(lapply(nint_seq, length))
  all_keys <- paste("Chr", rep(chrs, times = lint), unlist(nint_seq), sep = ".")
  n_total  <- length(all_keys)
  state    <- rep(0L, n_total); names(state) <- all_keys
  n_qtl    <- 2L
  qtl_keys <- all_keys[seq_len(n_qtl)]
  for (k in seq_len(n_qtl)) state[qtl_keys[k]] <- 1L

  eff_names <- sub("^Chr", "X", qtl_keys)
  effects  <- setNames(rnorm(n_qtl, 0.3, 0.05), eff_names)
  veffects <- setNames(runif(n_qtl, 0.02, 0.04), eff_names)

  QTL <- list(
    qtl = qtl_keys, effects = effects, veffects = veffects,
    method = "fixed", type = "marker", selection = "interval",
    TypeI = 0.05, iterations = n_qtl + 1L, breakout = FALSE,
    diag = list(state = state, genetic.term = "id", rel.scale = 1L,
                oint = list(), blups = list(), lik = list(),
                ochr = NULL, lik.mat = NULL,
                coef.list = list(), vcoef.list = list())
  )
  vpar <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
  vpar.con <- c(0L, 0L); names(vpar.con) <- names(vpar)
  mf <- data.frame(id = factor(ids), yld = rnorm(30), stringsAsFactors = FALSE)
  obj <- list(
    converge = TRUE, loglik = -45, sigma2 = 1.0,
    vparameters = vpar, vparameters.con = vpar.con,
    coefficients = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(30), 30, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id), mf = mf,
    call = call("asreml"), QTL = QTL
  )
  class(obj) <- c("qtlAim", "asreml")

  qtlm <- getQTL(obj, genObj)
  expect_equal(ncol(qtlm), 4L)
  expect_equal(nrow(qtlm), n_qtl)
  expect_true(is.character(qtlm))
})

test_that("getQTL marker-type: column 3 is a marker name, column 4 is cM-like", {
  set.seed(56)
  genObj <- make_wgCross_marker(n_lines = 20, n_chr = 2, n_mar = 4)
  chrs   <- names(genObj$geno)
  ids    <- paste0("L", seq_len(20))

  gdat     <- lapply(genObj$geno, function(el) el$imputed.data)
  nint_seq <- lapply(gdat, function(el) seq_len(ncol(el)))
  lint     <- unlist(lapply(nint_seq, length))
  all_keys <- paste("Chr", rep(chrs, times = lint), unlist(nint_seq), sep = ".")
  n_total  <- length(all_keys)
  state    <- rep(0L, n_total); names(state) <- all_keys
  qtl_keys <- all_keys[1L]
  state[qtl_keys] <- 1L
  eff_names <- sub("^Chr", "X", qtl_keys)
  effects  <- setNames(0.35, eff_names)
  veffects <- setNames(0.03, eff_names)

  QTL <- list(
    qtl = qtl_keys, effects = effects, veffects = veffects,
    method = "fixed", type = "marker", selection = "interval",
    TypeI = 0.05, iterations = 2L, breakout = FALSE,
    diag = list(state = state, genetic.term = "id", rel.scale = 1L,
                oint = list(), blups = list(), lik = list(),
                ochr = NULL, lik.mat = NULL,
                coef.list = list(), vcoef.list = list())
  )
  vpar <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
  vpar.con <- c(0L, 0L); names(vpar.con) <- names(vpar)
  mf <- data.frame(id = factor(ids), yld = rnorm(20), stringsAsFactors = FALSE)
  obj <- list(
    converge = TRUE, loglik = -40, sigma2 = 1.0,
    vparameters = vpar, vparameters.con = vpar.con,
    coefficients = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(20), 20, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf = mf, call = call("asreml"), QTL = QTL
  )
  class(obj) <- c("qtlAim", "asreml")

  qtlm <- getQTL(obj, genObj)
  # Column 3: non-empty string (marker name)
  expect_true(is.character(qtlm[1, 3]) && nchar(qtlm[1, 3]) > 0)
  # Column 4: numeric-like cM value
  expect_false(is.na(suppressWarnings(as.numeric(qtlm[1, 4]))))
})

# =============================================================================
# getQTL() — single-marker chromosome: lhmark == rhmark == imark
# =============================================================================

test_that("getQTL single-marker chr: left, inferred, and right marker names/positions are equal", {
  set.seed(57)
  n_lines <- 20
  ids     <- paste0("L", seq_len(n_lines))

  # Build an interval wgCross where C1 has exactly 1 marker
  genObj <- make_wgCross_interval(n_lines = n_lines, n_chr = 2, n_mar = 4)

  # Shrink C1 to a single marker
  single_map  <- genObj$geno[["C1"]]$map[1L]
  single_name <- names(single_map)
  genObj$geno[["C1"]]$map          <- single_map
  genObj$geno[["C1"]]$imputed.data <- genObj$geno[["C1"]]$imputed.data[, 1L, drop = FALSE]
  genObj$geno[["C1"]]$interval.data <- genObj$geno[["C1"]]$interval.data[, 1L, drop = FALSE]
  colnames(genObj$geno[["C1"]]$interval.data) <- single_name
  genObj$geno[["C1"]]$inferred.map <- setNames(single_map[1L], single_name)

  # Build an engine-style key for the single marker on C1 (index 1)
  qtl_key  <- "Chr.C1.1"
  eff_name <- "X.C1.1"

  # Also include C2 keys in the state vector for completeness
  c2_gdat <- lapply(list(genObj$geno[["C2"]]), function(el) el$interval.data)
  c2_nint <- seq_len(ncol(c2_gdat[[1]]))
  c2_keys <- paste("Chr", "C2", c2_nint, sep = ".")
  all_keys <- c(qtl_key, c2_keys)
  state    <- rep(0L, length(all_keys)); names(state) <- all_keys
  state[qtl_key] <- 1L

  effects  <- setNames(0.4, eff_name)
  veffects <- setNames(0.03, eff_name)

  QTL <- list(
    qtl = qtl_key, effects = effects, veffects = veffects,
    method = "fixed", type = "interval", selection = "interval",
    TypeI = 0.05, iterations = 2L, breakout = FALSE,
    diag = list(state = state, genetic.term = "id", rel.scale = 1L,
                oint = list(), blups = list(), lik = list(),
                ochr = NULL, lik.mat = NULL,
                coef.list = list(), vcoef.list = list())
  )
  vpar <- c(0.3, 1.0); names(vpar) <- c("id", "R!variance")
  vpar.con <- c(0L, 0L); names(vpar.con) <- names(vpar)
  mf  <- data.frame(id = factor(ids), yld = rnorm(n_lines), stringsAsFactors = FALSE)
  obj <- list(
    converge = TRUE, loglik = -42, sigma2 = 1.0,
    vparameters = vpar, vparameters.con = vpar.con,
    coefficients = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf = mf, call = call("asreml"), QTL = QTL
  )
  class(obj) <- c("qtlAim", "asreml")

  qtlm <- getQTL(obj, genObj)
  # For a single-marker chromosome: lhmark == rhmark == imark (col 3 = col 5 = col 7)
  expect_equal(qtlm[1, 3], qtlm[1, 5])
  expect_equal(qtlm[1, 3], qtlm[1, 7])
  # Positions (cols 4, 6, 8) should also be equal
  expect_equal(qtlm[1, 4], qtlm[1, 6])
  expect_equal(qtlm[1, 4], qtlm[1, 8])
})

# =============================================================================
# aimTable() — single qtlAim object
# =============================================================================

test_that("aimTable with one qtlAim returns aimTable data.frame with Trait column", {
  fix <- build_summable_qtlAim(n_qtl = 2)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "Grain weight", LOD = TRUE)

  expect_s3_class(at, "aimTable")
  expect_s3_class(at, "data.frame")
  expect_true("Trait" %in% names(at))
  expect_gte(nrow(at), 1L)
})

test_that("aimTable single qtlAim row count equals n_qtl", {
  fix <- build_summable_qtlAim(n_qtl = 2)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "Y")

  expect_equal(nrow(at), 2L)
})

# =============================================================================
# aimTable() — two qtlAim objects stacked
# =============================================================================

test_that("aimTable with two qtlAim objects stacks rows with two distinct Trait values", {
  fix1 <- build_summable_qtlAim(n_qtl = 1)
  fix2 <- build_summable_qtlAim(n_qtl = 2)

  at <- aimTable(fix1$obj, fix2$obj,
                 genObj = list(fix1$genObj, fix2$genObj),
                 labels = c("Trait1", "Trait2"))

  expect_s3_class(at, "aimTable")
  expect_equal(length(unique(at$Trait)), 2L)
  expect_equal(nrow(at), 1L + 2L)
})

# =============================================================================
# aimTable() — labels= overrides default trait labels
# =============================================================================

test_that("aimTable labels= argument overrides automatic trait labels", {
  fix      <- build_summable_qtlAim(n_qtl = 2)
  my_label <- "My Custom Trait"
  at       <- aimTable(fix$obj, genObj = fix$genObj, labels = my_label)

  expect_true(all(at$Trait == my_label))
})

# =============================================================================
# aimTable() — columns= subsets columns
# =============================================================================

test_that("aimTable columns= numeric vector results in Trait + requested N cols", {
  fix <- build_summable_qtlAim(n_qtl = 2)

  # columns = 1:2 means keep the first 2 columns of the summary output
  at_sub <- aimTable(fix$obj, genObj = fix$genObj,
                     labels = "T", LOD = FALSE, columns = 1:2)

  # Result: Trait column + 2 summary columns = 3 total
  expect_equal(ncol(at_sub), 3L)
})

# =============================================================================
# aimTable() — LOD= flag
# =============================================================================

test_that("aimTable LOD=TRUE appends a LOD column", {
  fix <- build_summable_qtlAim(n_qtl = 2)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "T", LOD = TRUE)

  expect_true("LOD" %in% names(at))
})

test_that("aimTable LOD=FALSE omits the LOD column", {
  fix <- build_summable_qtlAim(n_qtl = 2)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "T", LOD = FALSE)

  expect_false("LOD" %in% names(at))
})

# =============================================================================
# aimTable() — mixing qtlAim and gwasAim stops with error
# =============================================================================

test_that("aimTable stops when qtlAim and gwasAim objects are mixed", {
  fix_qtl  <- build_summable_qtlAim(n_qtl = 1)
  fix_gwas <- build_summable_gwasAim(n_qtl = 1)

  expect_error(
    aimTable(fix_qtl$obj, fix_gwas$obj,
             genObj = list(fix_qtl$genObj, fix_gwas$genObj)),
    regexp = "same class"
  )
})

# =============================================================================
# aimTable() — single gwasAim object
# =============================================================================

test_that("aimTable with gwasAim returns aimTable; obj.class attr is 'gwasAim'", {
  fix <- build_summable_gwasAim(n_qtl = 2)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "GWAS trait")

  expect_s3_class(at, "aimTable")
  expect_equal(attr(at, "obj.class"), "gwasAim")
  expect_gte(nrow(at), 1L)
  expect_true("Trait" %in% names(at))
})

# =============================================================================
# aimTable() — no objects supplied stops
# =============================================================================

test_that("aimTable stops when called with no model objects", {
  expect_error(
    aimTable(genObj = make_wgCross_interval()),
    regexp = "at least one"
  )
})

# =============================================================================
# aimTable() — missing genObj stops
# =============================================================================

test_that("aimTable stops when genObj is missing", {
  fix <- build_summable_qtlAim(n_qtl = 1)

  expect_error(
    aimTable(fix$obj, labels = "T"),
    regexp = "genObj.*required|required.*genObj"
  )
})

# =============================================================================
# aimTable() — object with no QTL: message emitted; that trait dropped
# =============================================================================

test_that("aimTable emits message and drops the trait when QTL$effects is NULL", {
  fix    <- build_summable_qtlAim(n_qtl = 2)
  # Clone and remove QTL effects to simulate no detection
  no_qtl <- fix$obj
  no_qtl$QTL$effects <- NULL

  expect_message(
    at <- aimTable(fix$obj, no_qtl,
                   genObj = list(fix$genObj, fix$genObj),
                   labels = c("WithQTL", "NoQTL")),
    regexp = "excluded|no detected QTL"
  )
  expect_true( "WithQTL" %in% at$Trait)
  expect_false("NoQTL"   %in% at$Trait)
})

# =============================================================================
# aimTable() — all objects have no QTL: empty aimTable returned with message
# =============================================================================

test_that("aimTable returns 0-row aimTable with message when all models have no QTL", {
  fix    <- build_summable_qtlAim(n_qtl = 1)
  empty1 <- fix$obj; empty1$QTL$effects <- NULL
  empty2 <- fix$obj; empty2$QTL$effects <- NULL

  expect_message(
    at <- aimTable(empty1, empty2,
                   genObj = fix$genObj,
                   labels = c("T1", "T2")),
    regexp = "[Nn]o model"
  )
  expect_s3_class(at, "aimTable")
  expect_equal(nrow(at), 0L)
})

# =============================================================================
# aimTable() — attr obj.class is 'qtlAim' for qtlAim inputs
# =============================================================================

test_that("aimTable attr obj.class is 'qtlAim' for qtlAim objects", {
  fix <- build_summable_qtlAim(n_qtl = 1)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "T")

  expect_equal(attr(at, "obj.class"), "qtlAim")
})

# =============================================================================
# print.aimTable() — populated table: returns invisible(x)
# =============================================================================

test_that("print.aimTable returns x invisibly", {
  fix <- build_summable_qtlAim(n_qtl = 2)
  at  <- aimTable(fix$obj, genObj = fix$genObj, labels = "T")

  result <- withVisible(print(at))
  expect_false(result$visible)
  expect_identical(result$value, at)
})

# =============================================================================
# print.aimTable() — 0-row table prints "(no QTL detected)" message
# =============================================================================

test_that("print.aimTable on 0-row table prints no-QTL notice", {
  empty_at           <- data.frame(Trait = character(0L))
  class(empty_at)    <- c("aimTable", "data.frame")
  attr(empty_at, "obj.class") <- "qtlAim"

  output <- capture.output(print(empty_at))
  expect_true(any(grepl("no QTL", output, ignore.case = TRUE)))
})

# =============================================================================
# print.aimTable() — class is preserved after printing
# =============================================================================

test_that("print.aimTable preserves the aimTable class of the returned value", {
  fix      <- build_summable_qtlAim(n_qtl = 2)
  at       <- aimTable(fix$obj, genObj = fix$genObj, labels = "T")
  returned <- print(at)

  expect_s3_class(returned, "aimTable")
})

# =============================================================================
# print.aimTable() — header contains analysis type and trait count
# =============================================================================

test_that("print.aimTable header contains obj.class and correct trait count", {
  fix1 <- build_summable_qtlAim(n_qtl = 1)
  fix2 <- build_summable_qtlAim(n_qtl = 1)
  at   <- aimTable(fix1$obj, fix2$obj,
                   genObj = list(fix1$genObj, fix2$genObj),
                   labels = c("TraitA", "TraitB"))

  output   <- capture.output(print(at))
  combined <- paste(output, collapse = " ")
  expect_true(grepl("qtlAim", combined))
  expect_true(grepl("2", combined))     # two distinct traits
})
