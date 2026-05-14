# =============================================================================
# test-linkMap.R
# Comprehensive tests for linkMap.cross(), linkMap.qtlAim(),
# linkMap.gwasAim(), linkMap.default(), theme_scatter(), and
# fineMap() guard clauses (no ASReml required).
#
# Fixtures supplied by helper-fixtures.R:
#   make_wgCross_interval(), make_wgCross_marker(), make_wgPanel(),
#   make_mock_qtlAim(), make_mock_gwasAim()
# =============================================================================

testthat::local_edition(3)

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixture instances (built once per test file load)
# ─────────────────────────────────────────────────────────────────────────────

wgcross_int    <- make_wgCross_interval()   # class: wgCross / bc / cross
wgcross_marker <- make_wgCross_marker()     # class: wgCross / bc / cross
panel          <- make_wgPanel()            # class: wgPanel / wgCross / cross

chr_names_int    <- names(wgcross_int$geno)    # "C1" "C2" "C3"
chr_names_panel  <- names(panel$geno)          # "P1" "P2"

# ─────────────────────────────────────────────────────────────────────────────
# 1. linkMap.cross — basic call (all chromosomes)
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: basic call with wgCross returns ggplot", {
  result <- linkMap(wgcross_int)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 2. linkMap.cross — chr= filter to one chromosome
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: chr= filter to one chromosome returns ggplot", {
  result <- linkMap(wgcross_int, chr = chr_names_int[1])
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 3. linkMap.cross — marker.names = "none"
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: marker.names='none' returns ggplot without labels", {
  result <- linkMap(wgcross_int, marker.names = "none")
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 4. linkMap.cross — marker.names = "markers"
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: marker.names='markers' returns ggplot", {
  result <- linkMap(wgcross_int, marker.names = "markers")
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 5. linkMap.cross — nrow = 2 (multiple chromosomes)
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: nrow=2 returns ggplot", {
  # wgcross_int has 3 chromosomes, so nrow=2 gives a 2-row layout
  result <- linkMap(wgcross_int, nrow = 2)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 6. linkMap.cross — nrow = "auto"
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: nrow='auto' returns ggplot", {
  result <- linkMap(wgcross_int, nrow = "auto")
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 7. linkMap.cross — row.chr list layout
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: row.chr list returns ggplot", {
  chrs <- chr_names_int            # "C1" "C2" "C3"
  rc   <- list(row1 = chrs[1:2], row2 = chrs[3])
  result <- linkMap(wgcross_int, row.chr = rc)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 8. linkMap.cross — chr.dist with valid $start / $end
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: chr.dist clips to cM range and returns ggplot", {
  # Map spans 0–80 cM; restrict to 10–60 cM
  result <- linkMap(wgcross_int, chr.dist = list(start = 10, end = 60))
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 9. linkMap.cross — bad chr name stops with error
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.cross: non-existent chr name stops with error", {
  expect_error(linkMap(wgcross_int, chr = "NOTACHR"),
               regexp = "do not match chromosome names")
})

# ─────────────────────────────────────────────────────────────────────────────
# 10. linkMap.qtlAim — single model, interval type → ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.qtlAim: single interval-type model returns ggplot", {
  qtl   <- make_mock_qtlAim()
  genObj <- attr(qtl, "genObj")
  result <- linkMap(qtl, genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 11. linkMap.qtlAim — no QTL effects → warning + ggplot (map only)
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.qtlAim: no QTL effects emits warning and returns ggplot", {
  qtl_noqtl         <- make_mock_qtlAim()
  genObj            <- attr(qtl_noqtl, "genObj")
  qtl_noqtl$QTL$effects <- NULL
  qtl_noqtl$QTL$qtl    <- NULL

  expect_warning(
    {
      result <- linkMap(qtl_noqtl, genObj = genObj)
    },
    regexp = "No significant QTL"
  )
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 12. linkMap.qtlAim — marker-type wgCross → ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.qtlAim: marker-type model returns ggplot", {
  qtl    <- make_mock_qtlAim(gen.type = "marker")
  # Override QTL type to "marker" and rebuild QTL keys to be consistent
  qtl$QTL$type <- "marker"
  genObj <- make_wgCross_marker()

  # Re-derive effects so they reference marker columns that exist in genObj
  chrs      <- names(genObj$geno)
  chr_idx   <- chrs[(seq_len(qtl$QTL$iterations - 1L) - 1L) %% length(chrs) + 1L]
  int_idx   <- seq_len(qtl$QTL$iterations - 1L) + 1L
  eff_names <- paste("X", chr_idx, int_idx, sep = ".")
  qtl$QTL$effects  <- stats::setNames(rnorm(length(eff_names), 0.4, 0.1), eff_names)
  qtl$QTL$veffects <- stats::setNames(runif(length(eff_names), 0.02, 0.05), eff_names)
  qtl$QTL$qtl      <- paste("Chr", chr_idx, int_idx, sep = ".")

  result <- linkMap(qtl, genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 13. linkMap.qtlAim — multi-model (two qtlAim objects) → ggplot with legend
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.qtlAim: multi-model returns ggplot", {
  qtl1   <- make_mock_qtlAim()
  qtl2   <- make_mock_qtlAim()
  # Give qtl2 a different trait name so labels differ
  qtl2$call$fixed <- quote(tgw ~ 1)
  genObj <- attr(qtl1, "genObj")

  result <- linkMap(qtl1, qtl2, genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 14. linkMap.qtlAim — flanking = FALSE → ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.qtlAim: flanking=FALSE returns ggplot", {
  qtl    <- make_mock_qtlAim()
  genObj <- attr(qtl, "genObj")
  result <- linkMap(qtl, genObj = genObj, flanking = FALSE)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 15. linkMap.gwasAim — single model → ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.gwasAim: single model returns ggplot", {
  gwas   <- make_mock_gwasAim()
  genObj <- attr(gwas, "genObj")
  result <- linkMap(gwas, genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 16. linkMap.gwasAim — multi-model (two gwasAim objects) → ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.gwasAim: multi-model returns ggplot", {
  gwas1  <- make_mock_gwasAim()
  gwas2  <- make_mock_gwasAim()
  gwas2$call$fixed <- quote(protein ~ 1)
  genObj <- attr(gwas1, "genObj")

  result <- linkMap(gwas1, gwas2, genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 17. linkMap.gwasAim — no significant markers → warning + ggplot (map only)
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.gwasAim: no significant markers emits warning and returns ggplot", {
  gwas_noqtl              <- make_mock_gwasAim()
  genObj                  <- attr(gwas_noqtl, "genObj")
  gwas_noqtl$QTL$effects  <- NULL
  gwas_noqtl$QTL$qtl      <- NULL

  expect_warning(
    {
      result <- linkMap(gwas_noqtl, genObj = genObj)
    },
    regexp = "No significant GWAS markers"
  )
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 18. linkMap.default — list of two qtlAim objects → dispatches and returns ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.default: list of two qtlAim objects returns ggplot", {
  qtl1   <- make_mock_qtlAim()
  qtl2   <- make_mock_qtlAim()
  qtl2$call$fixed <- quote(protein ~ 1)
  genObj <- attr(qtl1, "genObj")

  result <- linkMap(list(qtl1, qtl2), genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 19. linkMap.default — list of two gwasAim objects → dispatches and returns ggplot
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.default: list of two gwasAim objects returns ggplot", {
  gwas1  <- make_mock_gwasAim()
  gwas2  <- make_mock_gwasAim()
  gwas2$call$fixed <- quote(protein ~ 1)
  genObj <- attr(gwas1, "genObj")

  result <- linkMap(list(gwas1, gwas2), genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# ─────────────────────────────────────────────────────────────────────────────
# 20. linkMap.default — empty list → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.default: empty list stops with error", {
  expect_error(linkMap(list()), regexp = "empty list")
})

# ─────────────────────────────────────────────────────────────────────────────
# 21. linkMap.default — non-list → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.default: non-list (numeric) stops with error", {
  expect_error(linkMap(42), regexp = "No linkMap method")
})

# ─────────────────────────────────────────────────────────────────────────────
# 22. linkMap.default — mixed classes (one qtlAim + one gwasAim) → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.default: mixed qtlAim and gwasAim in list stops with error", {
  qtl  <- make_mock_qtlAim()
  gwas <- make_mock_gwasAim()
  expect_error(linkMap(list(qtl, gwas)), regexp = "same class")
})

# ─────────────────────────────────────────────────────────────────────────────
# 23. linkMap.default — wrong class in list (not qtlAim/gwasAim) → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("linkMap.default: unrecognised class in list stops with error", {
  bad_obj         <- list(QTL = list())
  class(bad_obj)  <- "unknownClass"
  # Two identical bad objects so the single-class check passes, then dispatch fails
  expect_error(linkMap(list(bad_obj, bad_obj)), regexp = "No multi-model linkMap method")
})

# ─────────────────────────────────────────────────────────────────────────────
# 24. theme_scatter — returns a ggplot2 theme object
# ─────────────────────────────────────────────────────────────────────────────

test_that("theme_scatter: returns an object of class 'theme'", {
  th <- wgAim:::theme_scatter()
  expect_s3_class(th, "theme")
})

# ─────────────────────────────────────────────────────────────────────────────
# 25. theme_scatter — base_size argument is accepted without error
# ─────────────────────────────────────────────────────────────────────────────

test_that("theme_scatter: base_size argument is accepted", {
  th <- wgAim:::theme_scatter(base_size = 14)
  expect_s3_class(th, "theme")
})

# ─────────────────────────────────────────────────────────────────────────────
# 26. fineMap guard — non-qtlAim/gwasAim object → stops with informative error
# ─────────────────────────────────────────────────────────────────────────────

test_that("fineMap: non-qtlAim/gwasAim object stops with error", {
  bad <- list(QTL = list())
  class(bad) <- "notAnAim"
  expect_error(fineMap(bad, genObj = wgcross_int),
               regexp = "qtlAim.*gwasAim|qtlAim|gwasAim")
})

# ─────────────────────────────────────────────────────────────────────────────
# 27. fineMap guard — object with no QTL ($QTL$qtl = NULL) → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("fineMap: object with no QTL stops with error", {
  qtl_noqtl      <- make_mock_qtlAim()
  qtl_noqtl$QTL$qtl <- NULL
  expect_error(fineMap(qtl_noqtl, genObj = attr(qtl_noqtl, "genObj")),
               regexp = "No significant QTL")
})

# ─────────────────────────────────────────────────────────────────────────────
# 28. fineMap guard — wrong genObj class for qtlAim (non-wgCross passed) → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("fineMap: non-wgCross genObj passed to qtlAim object stops with error", {
  qtl      <- make_mock_qtlAim()
  # wgPanel inherits wgCross, so fabricate a genObj that truly is NOT wgCross
  bad_genObj       <- list(geno = list(), pheno = data.frame())
  class(bad_genObj) <- "notWgCross"
  expect_error(fineMap(qtl, genObj = bad_genObj),
               regexp = "wgCross")
})

# ─────────────────────────────────────────────────────────────────────────────
# 29. fineMap guard — wrong genObj class for gwasAim (wgCross passed) → stops
# ─────────────────────────────────────────────────────────────────────────────

test_that("fineMap: wgCross genObj passed to gwasAim object stops with error", {
  gwas <- make_mock_gwasAim()
  expect_error(fineMap(gwas, genObj = wgcross_int),
               regexp = "wgPanel")
})

# ─────────────────────────────────────────────────────────────────────────────
# Additional edge-case tests
# ─────────────────────────────────────────────────────────────────────────────

# Test that linkMap.cross handles a wgPanel (which inherits "cross") correctly
test_that("linkMap.cross: wgPanel (inherits cross) returns ggplot", {
  result <- linkMap(panel)
  expect_s3_class(result, "ggplot")
})

# Test that chr.dist with only $start is accepted
test_that("linkMap.cross: chr.dist with only $start clips map and returns ggplot", {
  result <- linkMap(wgcross_int, chr.dist = list(start = 20))
  expect_s3_class(result, "ggplot")
})

# Test that chr.dist with only $end is accepted
test_that("linkMap.cross: chr.dist with only $end clips map and returns ggplot", {
  result <- linkMap(wgcross_int, chr.dist = list(end = 50))
  expect_s3_class(result, "ggplot")
})

# Test that chr.dist with bad element names errors
test_that("linkMap.cross: chr.dist with bad element name stops with error", {
  expect_error(
    linkMap(wgcross_int, chr.dist = list(begin = 0, finish = 80)),
    regexp = "chr.dist.*start.*end|start.*end"
  )
})

# Test that row.chr with no valid chromosomes errors
test_that("linkMap.cross: row.chr with all invalid chromosomes stops with error", {
  expect_error(
    linkMap(wgcross_int, row.chr = list(row1 = c("BADCHR1", "BADCHR2"))),
    regexp = "None of the chromosomes"
  )
})

# Test that row.chr with a non-character-vector element errors
test_that("linkMap.cross: row.chr with non-character element stops with error", {
  expect_error(
    linkMap(wgcross_int, row.chr = list(row1 = 1:3)),
    regexp = "list of character vectors"
  )
})

# linkMap.qtlAim — missing genObj → stops
test_that("linkMap.qtlAim: missing genObj stops with error", {
  qtl <- make_mock_qtlAim()
  expect_error(linkMap(qtl), regexp = "genObj.*required|required.*genObj")
})

# linkMap.qtlAim — wrong genObj class (not wgCross) → stops
test_that("linkMap.qtlAim: non-wgCross genObj stops with error", {
  qtl           <- make_mock_qtlAim()
  # wgPanel inherits wgCross, so use an object with no wgCross ancestry
  bad_genObj       <- list(geno = list(), pheno = data.frame())
  class(bad_genObj) <- "notWgCross"
  expect_error(linkMap(qtl, genObj = bad_genObj), regexp = "wgCross")
})

# linkMap.gwasAim — missing genObj → stops
test_that("linkMap.gwasAim: missing genObj stops with error", {
  gwas <- make_mock_gwasAim()
  expect_error(linkMap(gwas), regexp = "genObj.*required|required.*genObj")
})

# linkMap.gwasAim — wrong genObj class (wgCross) → stops
test_that("linkMap.gwasAim: wgCross passed as genObj stops with error", {
  gwas <- make_mock_gwasAim()
  expect_error(linkMap(gwas, genObj = wgcross_int), regexp = "wgPanel")
})

# linkMap.qtlAim — chr= filter reduces displayed chromosomes, still ggplot
test_that("linkMap.qtlAim: chr= argument returns ggplot for a subset of chromosomes", {
  qtl    <- make_mock_qtlAim()
  genObj <- attr(qtl, "genObj")
  # Show only the first chromosome from the QTL-bearing set
  qtl_chrs <- names(genObj$geno)[1]
  result <- linkMap(qtl, genObj = genObj, chr = qtl_chrs)
  expect_s3_class(result, "ggplot")
})

# linkMap.gwasAim — marker.names = "none" → ggplot
test_that("linkMap.gwasAim: marker.names='none' returns ggplot", {
  gwas   <- make_mock_gwasAim()
  genObj <- attr(gwas, "genObj")
  result <- linkMap(gwas, genObj = genObj, marker.names = "none")
  expect_s3_class(result, "ggplot")
})

# linkMap.gwasAim — marker.names = "dist" → ggplot
test_that("linkMap.gwasAim: marker.names='dist' returns ggplot", {
  gwas   <- make_mock_gwasAim()
  genObj <- attr(gwas, "genObj")
  result <- linkMap(gwas, genObj = genObj, marker.names = "dist")
  expect_s3_class(result, "ggplot")
})

# linkMap.qtlAim single-model: custom trait label → ggplot
test_that("linkMap.qtlAim: custom trait.labels returned in ggplot", {
  qtl    <- make_mock_qtlAim()
  genObj <- attr(qtl, "genObj")
  result <- linkMap(qtl, genObj = genObj, trait.labels = "MyTrait")
  expect_s3_class(result, "ggplot")
})

# linkMap.gwasAim single-model: custom trait label → ggplot
test_that("linkMap.gwasAim: custom trait.labels returned in ggplot", {
  gwas   <- make_mock_gwasAim()
  genObj <- attr(gwas, "genObj")
  result <- linkMap(gwas, genObj = genObj, trait.labels = "MyTrait")
  expect_s3_class(result, "ggplot")
})

# linkMap.default — single-element qtlAim list (edge case: no extras)
test_that("linkMap.default: single-element qtlAim list returns ggplot", {
  qtl    <- make_mock_qtlAim()
  genObj <- attr(qtl, "genObj")
  result <- linkMap(list(qtl), genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# linkMap.default — single-element gwasAim list (edge case: no extras)
test_that("linkMap.default: single-element gwasAim list returns ggplot", {
  gwas   <- make_mock_gwasAim()
  genObj <- attr(gwas, "genObj")
  result <- linkMap(list(gwas), genObj = genObj)
  expect_s3_class(result, "ggplot")
})

# fineMap — empty QTL vector (length 0) also stops
test_that("fineMap: zero-length QTL$qtl stops with error", {
  qtl_empty       <- make_mock_qtlAim()
  qtl_empty$QTL$qtl <- character(0L)
  expect_error(fineMap(qtl_empty, genObj = attr(qtl_empty, "genObj")),
               regexp = "No significant QTL")
})
