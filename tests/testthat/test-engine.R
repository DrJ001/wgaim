# =============================================================================
# test-engine.R
# Tests for internal engine functions:
#   .validateModel  (engine_validate.R)
#   .buildGenoData  (engine_genodata.R)
#   .lrtTest        (engine_select.R)
#   .packResults    (engine_results.R)
#   .mergeEffect    (engine_effects.R)
# =============================================================================

# helper-fixtures.R is loaded automatically by testthat

# Aliases for internal functions so tests are readable
.validateModel <- wgAim:::.validateModel
.buildGenoData <- wgAim:::.buildGenoData
.lrtTest       <- wgAim:::.lrtTest
.packResults   <- wgAim:::.packResults
.mergeEffect   <- wgAim:::.mergeEffect

# =============================================================================
# .validateModel
# =============================================================================

# Helper: build a mock model whose call$data evaluates to `data_val` directly.
# .validateModel calls eval(baseModel$call$data) inside the package, so we
# must embed the data frame as a literal object in the call rather than as a
# symbol name.  We do this by storing a self-evaluating quoted form via
# as.call(list(`function`() NULL)), or more simply by replacing call$data with
# the data frame wrapped so eval() returns it immediately.
#
# The cleanest approach: replace call$data with a call to `local(data_val)`
# where data_val is spliced in via bquote — this evaluates to the data frame
# regardless of the calling frame.
.make_validate_mock <- function(merge.by = "id") {
  df <- data.frame(
    id  = factor(paste0("L", 1:10)),
    yld = rnorm(10),
    stringsAsFactors = FALSE
  )
  m <- make_mock_asreml(data_name = "phenoData", data_val = df,
                        merge.by = merge.by)
  # Replace call$data with `local(.(df))` so eval() returns df unconditionally
  m$call$data <- bquote(local(.(df)))
  m
}

test_that(".validateModel: converged model passes and returns correct slots", {
  m <- .make_validate_mock()

  result <- .validateModel(m, merge.by = "id", method = "fixed",
                           selection = "interval", breakout = -1)

  expect_type(result, "list")
  expect_named(result, c("baseModel", "asremlEnv", "phenoData"))
  expect_true(is.list(result$asremlEnv))
  expect_true(is.data.frame(result$phenoData))
})

test_that(".validateModel: merge.by = NULL stops with informative message", {
  m <- .make_validate_mock()

  expect_error(
    .validateModel(m, merge.by = NULL, method = "fixed",
                   selection = "interval", breakout = -1),
    regexp = "merge\\.by"
  )
})

test_that(".validateModel: bad method stops with informative message", {
  m <- .make_validate_mock()

  expect_error(
    .validateModel(m, merge.by = "id", method = "bad",
                   selection = "interval", breakout = -1),
    regexp = "method must be"
  )
})

test_that(".validateModel: bad selection stops with informative message", {
  m <- .make_validate_mock()

  expect_error(
    .validateModel(m, merge.by = "id", method = "fixed",
                   selection = "bad", breakout = -1),
    regexp = "selection must be"
  )
})

test_that(".validateModel: breakout = 0 stops with informative message", {
  m <- .make_validate_mock()

  expect_error(
    .validateModel(m, merge.by = "id", method = "fixed",
                   selection = "interval", breakout = 0),
    regexp = "breakout must be"
  )
})

test_that(".validateModel: breakout = -2 stops with informative message", {
  m <- .make_validate_mock()

  expect_error(
    .validateModel(m, merge.by = "id", method = "fixed",
                   selection = "interval", breakout = -2),
    regexp = "breakout must be"
  )
})

test_that(".validateModel: fractional breakout (0.5) passes (not rejected by guard)", {
  # Re-check the actual condition:
  #   !is.numeric(breakout) | breakout < -1 | breakout == 0
  # 0.5 is numeric, not < -1, not == 0 → does NOT error.
  # The guard only rejects 0 and values strictly less than -1.
  m <- .make_validate_mock()

  expect_no_error(
    .validateModel(m, merge.by = "id", method = "fixed",
                   selection = "interval", breakout = 0.5)
  )
})

test_that(".validateModel: breakout = 2 (positive integer) passes", {
  m <- .make_validate_mock()

  expect_no_error(
    .validateModel(m, merge.by = "id", method = "fixed",
                   selection = "interval", breakout = 2)
  )
})

test_that(".validateModel: breakout = -1 (default sentinel) passes", {
  m <- .make_validate_mock()

  result <- .validateModel(m, merge.by = "id", method = "random",
                           selection = "chromosome", breakout = -1)
  expect_type(result, "list")
})

# =============================================================================
# .buildGenoData
# =============================================================================

test_that(".buildGenoData (interval): returns list with correct slots and matrix", {
  genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 4)
  glines <- plines <- paste0("L", 1:20)

  result <- .buildGenoData(genObj, gen.type = "interval",
                           glines = glines, plines = plines)

  expect_type(result, "list")
  expect_named(result, c("genoData", "mnams", "state"))
  expect_true(is.matrix(result$genoData))
})

test_that(".buildGenoData (marker): uses imputed.data; column names start with 'Chr.'", {
  genObj <- make_wgCross_marker(n_lines = 20, n_chr = 2, n_mar = 4)
  glines <- plines <- paste0("L", 1:20)

  result <- .buildGenoData(genObj, gen.type = "marker",
                           glines = glines, plines = plines)

  expect_true(is.matrix(result$genoData))
  # Column names are built as "Chr.<chrname>.<index>"
  expect_true(all(startsWith(colnames(result$genoData), "Chr.")))
})

test_that(".buildGenoData: state vector length equals ncol(genoData) and all values 1", {
  genObj <- make_wgCross_interval(n_lines = 30, n_chr = 3, n_mar = 5)
  glines <- plines <- paste0("L", 1:30)

  result <- .buildGenoData(genObj, gen.type = "interval",
                           glines = glines, plines = plines)

  expect_equal(length(result$state), ncol(result$genoData))
  expect_true(all(result$state == 1L))
})

test_that(".buildGenoData: lines not in plines are excluded from genoData rows", {
  genObj  <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 4)
  glines  <- paste0("L", 1:20)
  # Only keep first 10 lines in phenotypic data
  plines  <- paste0("L", 1:10)

  result <- .buildGenoData(genObj, gen.type = "interval",
                           glines = glines, plines = plines)

  expect_equal(nrow(result$genoData), 10L)
  expect_true(all(rownames(result$genoData) %in% plines))
})

test_that(".buildGenoData (interval): mnams length equals ncol(genoData)", {
  genObj <- make_wgCross_interval(n_lines = 20, n_chr = 2, n_mar = 4)
  glines <- plines <- paste0("L", 1:20)

  result <- .buildGenoData(genObj, gen.type = "interval",
                           glines = glines, plines = plines)

  expect_equal(length(result$mnams), ncol(result$genoData))
})

# =============================================================================
# .lrtTest
# =============================================================================

# Helper: build a minimal model-like list with only $loglik
.make_ll <- function(loglik) list(loglik = loglik)

test_that(".lrtTest: stat > threshold → pass = TRUE", {
  TypeI    <- 0.05
  thresh   <- qchisq(1 - 2 * TypeI, 1)
  base_ll  <- -50
  # stat = 2*(qtl_ll - base_ll); choose stat > thresh
  qtl_ll   <- base_ll + (thresh + 5) / 2

  res <- .lrtTest(.make_ll(qtl_ll), .make_ll(base_ll), TypeI)

  expect_true(res$pass)
  expect_gt(res$stat, thresh)
})

test_that(".lrtTest: stat < threshold → pass = FALSE", {
  TypeI   <- 0.05
  thresh  <- qchisq(1 - 2 * TypeI, 1)
  base_ll <- -50
  # stat = 2*(qtl_ll - base_ll); choose stat < thresh
  qtl_ll  <- base_ll + (thresh - 1) / 2

  res <- .lrtTest(.make_ll(qtl_ll), .make_ll(base_ll), TypeI)

  expect_false(res$pass)
})

test_that(".lrtTest: pvalue equals (1 - pchisq(stat, 1)) / 2", {
  TypeI   <- 0.05
  base_ll <- -60
  qtl_ll  <- -55   # stat = 10
  expected_stat <- 2 * (qtl_ll - base_ll)

  res <- .lrtTest(.make_ll(qtl_ll), .make_ll(base_ll), TypeI)

  expect_equal(res$stat, expected_stat)
  expect_equal(res$pvalue, (1 - pchisq(expected_stat, 1)) / 2)
})

test_that(".lrtTest: stat = 0 gives pvalue = 0.5 and pass = FALSE", {
  TypeI   <- 0.05
  base_ll <- -50
  qtl_ll  <- -50   # stat = 0

  res <- .lrtTest(.make_ll(qtl_ll), .make_ll(base_ll), TypeI)

  expect_equal(res$stat,   0)
  expect_equal(res$pvalue, 0.5)
  expect_false(res$pass)
})

test_that(".lrtTest: baseLogL slot equals base model loglik", {
  base_ll <- -42.7
  qtl_ll  <- -40.0

  res <- .lrtTest(.make_ll(qtl_ll), .make_ll(base_ll), TypeI = 0.05)

  expect_equal(res$baseLogL, base_ll)
})

# =============================================================================
# .packResults
# =============================================================================

# Build minimal inputs for .packResults
.make_pack_inputs <- function(n_qtl = 2, n_ints = 6) {
  set.seed(99)
  qtl_keys  <- paste("Chr", paste0("C", seq_len(n_qtl)),
                     seq_len(n_qtl) + 1L, sep = ".")
  eff_names <- paste("X",  paste0("C", seq_len(n_qtl)),
                     seq_len(n_qtl) + 1L, sep = ".")
  effects   <- setNames(rnorm(n_qtl, 0.5, 0.1), eff_names)
  veffects  <- setNames(runif(n_qtl, 0.02, 0.05), eff_names)

  coef.list  <- lapply(seq_len(n_qtl), function(k) effects[seq_len(k)])
  vcoef.list <- lapply(seq_len(n_qtl), function(k) veffects[seq_len(k)])

  all_keys <- paste("Chr", rep(paste0("C", seq_len(3)), each = 2),
                    rep(1:2, times = 3), sep = ".")
  state <- rep(1L, n_ints); names(state) <- all_keys

  lik_entries <- lapply(seq_len(n_qtl), function(k)
    list(baseLogL = -50 + k, stat = 4 + k, pvalue = 0.02, pass = TRUE))
  # Non-significant final iteration
  lik_entries[[n_qtl + 1]] <- list(baseLogL = -48, stat = 0.5, pvalue = 0.3,
                                   pass = FALSE)

  ldiag <- list(lik = lik_entries)

  list(
    qtl        = qtl_keys,
    coef.list  = coef.list,
    vcoef.list = vcoef.list,
    ldiag      = ldiag,
    state      = state,
    iter       = n_qtl,
    breakout   = -1L,
    cov.env    = NULL,
    genetic.term = "id",
    method     = "fixed",
    gen.type   = "interval",
    selection  = "interval",
    TypeI      = 0.05
  )
}

test_that(".packResults: non-empty qtl returns list with required slots", {
  inp <- .make_pack_inputs(n_qtl = 2, n_ints = 6)

  out <- do.call(.packResults, inp)
  res <- out$qtl.list          # unwrap: .packResults returns list(qtl.list, qtlModel.pruned)

  expect_type(out, "list")
  expect_true(!is.null(res$qtl))
  expect_true(!is.null(res$effects))
  expect_true(!is.null(res$veffects))
  expect_true(!is.null(res$diag$lik.mat))
  expect_true(!is.null(res$diag$state))
})

test_that(".packResults: lik.mat has four columns with correct names", {
  inp <- .make_pack_inputs(n_qtl = 2, n_ints = 6)
  out <- do.call(.packResults, inp)
  res <- out$qtl.list

  expect_equal(ncol(res$diag$lik.mat), 4L)
  expect_equal(colnames(res$diag$lik.mat),
               c("L0", "L1", "Statistic", "Pvalue"))
})

test_that(".packResults: empty qtl returns list with settings only, no effects", {
  inp       <- .make_pack_inputs(n_qtl = 2, n_ints = 6)
  inp$qtl   <- character(0)            # empty -> no QTL detected
  # iter must still be 1 for a degenerate run
  inp$iter  <- 1L
  # ldiag needs at least one entry
  inp$ldiag <- list(lik = list(list(baseLogL = -50, stat = 0.5,
                                    pvalue = 0.3, pass = FALSE)))

  out <- do.call(.packResults, inp)
  res <- out$qtl.list

  expect_named(res, c("selection", "method", "type", "TypeI", "diag",
                       "iterations"),
               ignore.order = TRUE)
  expect_null(res$effects)
  expect_null(res$qtl)
})

test_that(".packResults: breakout slot is logical (FALSE when breakout == -1)", {
  inp <- .make_pack_inputs(n_qtl = 2, n_ints = 6)
  out <- do.call(.packResults, inp)
  res <- out$qtl.list

  expect_type(res$breakout, "logical")
  expect_false(res$breakout)   # breakout = -1 -> breakout != -1 is FALSE
})

test_that(".packResults: breakout slot is TRUE when breakout > 0", {
  inp          <- .make_pack_inputs(n_qtl = 2, n_ints = 6)
  inp$breakout <- 3L
  out          <- do.call(.packResults, inp)
  res          <- out$qtl.list

  expect_true(res$breakout)
})

test_that(".packResults: settings slots are correctly passed through", {
  inp <- .make_pack_inputs(n_qtl = 2, n_ints = 6)
  out <- do.call(.packResults, inp)
  res <- out$qtl.list

  expect_equal(res$method,    "fixed")
  expect_equal(res$type,      "interval")
  expect_equal(res$selection, "interval")
  expect_equal(res$TypeI,     0.05)
})

# =============================================================================
# .mergeEffect
# =============================================================================

test_that(".mergeEffect: returns list with $phenoData and $qtl.x", {
  set.seed(42)
  ids      <- paste0("L", 1:20)
  phenoData <- data.frame(
    id  = factor(ids),
    yld = rnorm(20),
    stringsAsFactors = FALSE
  )
  genoData        <- matrix(sample(c(-1, 1), 20 * 5, replace = TRUE),
                            nrow = 20, dimnames = list(ids, paste0("Chr.C1.", 1:5)))
  qtl_name        <- "Chr.C1.3"

  res <- .mergeEffect(phenoData, genoData, qtl = qtl_name, merge.by = "id")

  expect_type(res, "list")
  expect_named(res, c("phenoData", "qtl.x"))
})

test_that(".mergeEffect: qtl.x has 'X.' prefix replacing 'Chr.'", {
  set.seed(42)
  ids      <- paste0("L", 1:20)
  phenoData <- data.frame(
    id  = factor(ids),
    yld = rnorm(20),
    stringsAsFactors = FALSE
  )
  genoData   <- matrix(sample(c(-1, 1), 20 * 5, replace = TRUE),
                       nrow = 20, dimnames = list(ids, paste0("Chr.C1.", 1:5)))
  qtl_name   <- "Chr.C1.3"

  res <- .mergeEffect(phenoData, genoData, qtl = qtl_name, merge.by = "id")

  expect_equal(res$qtl.x, "X.C1.3")
})

test_that(".mergeEffect: new column added to phenoData matches genoData values", {
  set.seed(42)
  ids       <- paste0("L", 1:20)
  phenoData <- data.frame(
    id  = factor(ids),
    yld = rnorm(20),
    stringsAsFactors = FALSE
  )
  geno_col  <- sample(c(-1, 1), 20, replace = TRUE)
  genoData  <- matrix(geno_col, nrow = 20,
                      dimnames = list(ids, "Chr.C1.1"))
  qtl_name  <- "Chr.C1.1"

  res <- .mergeEffect(phenoData, genoData, qtl = qtl_name, merge.by = "id")

  new_col <- "X.C1.1"
  expect_true(new_col %in% names(res$phenoData))

  # Values should match for lines present in both
  merged <- res$phenoData
  matched_vals <- merged[[new_col]][match(ids, as.character(merged$id))]
  expect_equal(matched_vals, geno_col, ignore_attr = TRUE)
})

test_that(".mergeEffect: phenoData row order is preserved after merge", {
  set.seed(7)
  ids       <- paste0("L", 1:15)
  phenoData <- data.frame(
    id  = factor(ids),
    yld = rnorm(15),
    stringsAsFactors = FALSE
  )
  genoData  <- matrix(sample(c(-1, 1), 15, replace = TRUE),
                      nrow = 15, dimnames = list(ids, "Chr.C2.1"))

  res <- .mergeEffect(phenoData, genoData, qtl = "Chr.C2.1", merge.by = "id")

  # Original row order must be retained (ord column is dropped but used for sort)
  expect_equal(as.character(res$phenoData$id), ids)
})
