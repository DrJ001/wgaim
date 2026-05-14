# =============================================================================
# test-print-summary.R
# Tests for S3 print/summary/aimTrace methods:
#   print.qtlAim      (print.qtlAim.R)
#   summary.qtlAim    (summary.qtlAim.R)
#   print.gwasAim     (print.gwasAim.R)
#   summary.gwasAim   (summary.gwasAim.R)
#   print.gpAim       (print.gpAim.R)
#   summary.gpAim     (summary.gpAim.R)
#   aimTrace.qtlAim   (aimTrace.qtlAim.R)
#   aimTrace.gwasAim  (aimTrace.gwasAim.R)
# =============================================================================

# helper-fixtures.R is loaded automatically by testthat

# =============================================================================
# print.qtlAim
# =============================================================================

test_that("print.qtlAim: missing genObj stops with informative message", {
  obj <- make_mock_qtlAim(n_qtl = 2)
  expect_error(print(obj), regexp = "genObj")
})

test_that("print.qtlAim: wrong class for genObj stops", {
  obj    <- make_mock_qtlAim(n_qtl = 2)
  bad_gn <- list(geno = list(), pheno = data.frame())   # not wgCross
  expect_error(print(obj, genObj = bad_gn), regexp = "wgCross")
})

test_that("print.qtlAim: no QTL (effects = NULL) prints 'no significant' message", {
  obj             <- make_mock_qtlAim(n_qtl = 2)
  obj$QTL$effects <- NULL
  genObj          <- attr(obj, "genObj")

  out <- capture.output(print(obj, genObj = genObj))
  expect_true(any(grepl("no significant", tolower(out))))
})

test_that("print.qtlAim (interval type): output contains chromosome name", {
  obj    <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
  genObj <- attr(obj, "genObj")

  out <- capture.output(print(obj, genObj = genObj))
  # Should mention at least one chromosome
  chr_names <- names(genObj$geno)
  expect_true(any(sapply(chr_names, function(cn) any(grepl(cn, out)))))
})

test_that("print.qtlAim (marker type): output contains chromosome name", {
  obj_m  <- make_mock_qtlAim(n_qtl = 1, gen.type = "marker")
  # For marker type we need a marker-typed genObj; rebuild QTL$type only.
  obj_m$QTL$type <- "marker"
  genObj <- attr(obj_m, "genObj")

  out <- capture.output(print(obj_m, genObj = genObj))
  chr_names <- names(genObj$geno)
  expect_true(any(sapply(chr_names, function(cn) any(grepl(cn, out)))))
})

# =============================================================================
# summary.qtlAim
# =============================================================================

test_that("summary.qtlAim: missing genObj stops", {
  obj <- make_mock_qtlAim(n_qtl = 2)
  expect_error(summary(obj), regexp = "genObj")
})

test_that("summary.qtlAim: no QTL returns invisible NULL and prints message", {
  obj             <- make_mock_qtlAim(n_qtl = 2)
  obj$QTL$effects <- NULL
  genObj          <- attr(obj, "genObj")

  out <- capture.output(result <- summary(obj, genObj = genObj))
  expect_null(result)
  expect_true(length(out) > 0)
})

test_that("summary.qtlAim: LOD = TRUE adds a LOD column", {
  obj    <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = TRUE)
  expect_true("LOD" %in% names(tab))
})

test_that("summary.qtlAim: LOD = FALSE omits the LOD column", {
  obj    <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  expect_false("LOD" %in% names(tab))
})

test_that("summary.qtlAim (method = 'fixed'): pvalue column is character, values in [0,1]", {
  obj         <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval", method = "fixed")
  obj$QTL$method <- "fixed"
  genObj      <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  # pvalue may be stored as character ("<0.0001") or numeric
  pv <- tab$Prob
  numeric_pv <- suppressWarnings(as.numeric(pv))
  non_special <- !is.na(numeric_pv)
  if (any(non_special))
    expect_true(all(numeric_pv[non_special] >= 0 & numeric_pv[non_special] <= 1))
})

test_that("summary.qtlAim (method = 'random'): pvalue column exists and is character/numeric", {
  obj            <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval", method = "random")
  obj$QTL$method <- "random"
  genObj         <- attr(obj, "genObj")

  # For method = "random", summary.qtlAim greps vparameters for "X\\." to find
  # QTL variance estimates.  We need to add those entries to the mock.
  qtl_eff_names <- names(obj$QTL$effects)             # e.g. "X.C2.2" "X.C3.3"
  qtl_vpar      <- setNames(rep(0.1, length(qtl_eff_names)), qtl_eff_names)
  obj$vparameters     <- c(obj$vparameters, qtl_vpar)
  obj$vparameters.con <- c(obj$vparameters.con,
                           setNames(rep(0L, length(qtl_eff_names)), qtl_eff_names))

  tab <- suppressWarnings(summary(obj, genObj = genObj, LOD = FALSE))
  expect_true("Prob" %in% names(tab))
})

test_that("summary.qtlAim (interval type): data.frame has at least 10 columns (without LOD)", {
  obj    <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  # Interval type: Chromosome, Left Marker, dist(cM), Infer. Marker, dist(cM),
  #                Right Marker, dist(cM), Size, Prob, Perc.Var  = 10 cols
  expect_gte(ncol(tab), 10L)
})

test_that("summary.qtlAim (marker type): data.frame has 6 columns without LOD", {
  obj         <- make_mock_qtlAim(n_qtl = 2, gen.type = "interval")
  obj$QTL$type <- "marker"
  genObj      <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  # Marker type: Chromosome, Marker, dist(cM), Size, Prob, Perc.Var = 6 cols
  expect_equal(ncol(tab), 6L)
})

test_that("summary.qtlAim: table is sorted by chromosome then cM position", {
  # Build object with 2+ QTL spanning different chromosomes
  obj    <- make_mock_qtlAim(n_qtl = 3, gen.type = "interval")
  genObj <- attr(obj, "genObj")

  # suppressWarnings: mock vparameters may have mismatched lengths in var calc
  tab <- suppressWarnings(summary(obj, genObj = genObj, LOD = FALSE))

  if (nrow(tab) >= 2) {
    chr_lead <- as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1",
                               as.character(tab$Chromosome)))
    expect_equal(chr_lead, sort(chr_lead))
  }
})

test_that("summary.qtlAim: returned object is a data.frame", {
  obj    <- make_mock_qtlAim(n_qtl = 2)
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj)
  expect_s3_class(tab, "data.frame")
})

# =============================================================================
# print.gwasAim
# =============================================================================

test_that("print.gwasAim: missing genObj stops", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  expect_error(print(obj), regexp = "genObj")
})

test_that("print.gwasAim: wrong class for genObj stops", {
  obj   <- make_mock_gwasAim(n_qtl = 2)
  bad_g <- list(geno = list(), pheno = data.frame())  # not wgPanel
  expect_error(print(obj, genObj = bad_g), regexp = "wgPanel")
})

test_that("print.gwasAim: output contains 'TypeI' and 'markers'", {
  obj    <- make_mock_gwasAim(n_qtl = 2)
  genObj <- attr(obj, "genObj")

  out <- capture.output(print(obj, genObj = genObj))
  expect_true(any(grepl("TypeI",   out, ignore.case = TRUE)))
  expect_true(any(grepl("marker",  out, ignore.case = TRUE)))
})

test_that("print.gwasAim: no QTL output contains 'No significant'", {
  obj             <- make_mock_gwasAim(n_qtl = 2)
  obj$QTL$effects <- NULL
  genObj          <- attr(obj, "genObj")

  out <- capture.output(print(obj, genObj = genObj))
  expect_true(any(grepl("No significant", out, ignore.case = TRUE)))
})

test_that("print.gwasAim: with QTL output contains chromosome name", {
  obj    <- make_mock_gwasAim(n_qtl = 2)
  genObj <- attr(obj, "genObj")

  out      <- capture.output(print(obj, genObj = genObj))
  chr_names <- names(genObj$geno)
  expect_true(any(sapply(chr_names, function(cn) any(grepl(cn, out)))))
})

# =============================================================================
# summary.gwasAim
# =============================================================================

test_that("summary.gwasAim: missing genObj stops", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  expect_error(summary(obj), regexp = "genObj")
})

test_that("summary.gwasAim: LOD = TRUE returns data.frame with LOD column", {
  obj    <- make_mock_gwasAim(n_qtl = 2)
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = TRUE)
  expect_s3_class(tab, "data.frame")
  expect_true("LOD" %in% names(tab))
})

test_that("summary.gwasAim: LOD = FALSE omits LOD column", {
  obj    <- make_mock_gwasAim(n_qtl = 2)
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  expect_false("LOD" %in% names(tab))
})

test_that("summary.gwasAim: correct column names for marker type", {
  obj    <- make_mock_gwasAim(n_qtl = 2)
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  expected_base <- c("Chromosome", "Marker", "dist(cM)", "Size", "Pvalue", "Perc.Var")
  expect_equal(names(tab), expected_base)
})

test_that("summary.gwasAim: no QTL returns invisible NULL and prints message", {
  obj             <- make_mock_gwasAim(n_qtl = 2)
  obj$QTL$effects <- NULL
  genObj          <- attr(obj, "genObj")

  out <- capture.output(result <- summary(obj, genObj = genObj))
  expect_null(result)
  expect_true(length(out) > 0)
})

test_that("summary.gwasAim: result has correct number of rows (one per QTL)", {
  n_qtl  <- 2
  obj    <- make_mock_gwasAim(n_qtl = n_qtl)
  genObj <- attr(obj, "genObj")

  tab <- summary(obj, genObj = genObj, LOD = FALSE)
  expect_equal(nrow(tab), n_qtl)
})

# =============================================================================
# print.gpAim
# =============================================================================

test_that("print.gpAim: output contains 'Genomic' and 'GEBV'", {
  obj <- make_mock_gpAim()
  out <- capture.output(print(obj))
  expect_true(any(grepl("Genomic",    out, ignore.case = TRUE)))
  expect_true(any(grepl("GEBV|gebv",  out, ignore.case = TRUE)))
})

test_that("print.gpAim: returns the object invisibly", {
  obj    <- make_mock_gpAim()
  result <- withVisible(print(obj))
  # invisible() means visible = FALSE
  expect_false(result$visible)
  expect_identical(result$value, obj)
})

test_that("print.gpAim: output contains heritability figure", {
  obj <- make_mock_gpAim()
  out <- capture.output(print(obj))
  expect_true(any(grepl("erit|h2|h\\^2", out, ignore.case = TRUE)))
})

test_that("print.gpAim: output contains variance component labels", {
  obj <- make_mock_gpAim()
  out <- capture.output(print(obj))
  expect_true(any(grepl("Var|variance", out, ignore.case = TRUE)))
})

# =============================================================================
# summary.gpAim
# =============================================================================

test_that("summary.gpAim: returns a data.frame sorted descending by GEBV", {
  obj <- make_mock_gpAim(n_lines = 30)
  tab <- summary(obj)

  expect_s3_class(tab, "data.frame")
  # Check descending order
  expect_equal(tab$GEBV, sort(tab$GEBV, decreasing = TRUE))
})

test_that("summary.gpAim: first row has the highest GEBV", {
  obj <- make_mock_gpAim(n_lines = 20)
  tab <- summary(obj)

  expect_equal(tab$GEBV[1], max(tab$GEBV))
})

test_that("summary.gpAim: returned data.frame contains GEBV and SE columns", {
  obj <- make_mock_gpAim()
  tab <- summary(obj)

  expect_true("GEBV" %in% names(tab))
  expect_true("SE"   %in% names(tab))
})

test_that("summary.gpAim: number of rows equals number of genotyped lines", {
  n_lines <- 25
  obj     <- make_mock_gpAim(n_lines = n_lines)
  tab     <- summary(obj)

  expect_equal(nrow(tab), n_lines)
})

test_that("summary.gpAim: printed output contains heritability and counts", {
  obj <- make_mock_gpAim(n_lines = 20, n_mar = 4)
  out <- capture.output(summary(obj))

  expect_true(any(grepl("h2|herit", out, ignore.case = TRUE)))
  expect_true(any(grepl("lines|markers", out, ignore.case = TRUE)))
})

# =============================================================================
# aimTrace.qtlAim
# =============================================================================

test_that("aimTrace.qtlAim: plot = FALSE returns invisible NULL", {
  obj    <- make_mock_qtlAim(n_qtl = 2)
  result <- aimTrace(obj, plot = FALSE, lik.out = FALSE)
  expect_null(result)
})

test_that("aimTrace.qtlAim: lik.out = TRUE prints the LRT table", {
  obj <- make_mock_qtlAim(n_qtl = 2)
  expect_output(
    aimTrace(obj, plot = FALSE, lik.out = TRUE),
    regexp = "Likelihood Ratio"
  )
})

test_that("aimTrace.qtlAim: always prints incremental p-value matrix", {
  obj <- make_mock_qtlAim(n_qtl = 2)
  expect_output(
    aimTrace(obj, plot = FALSE, lik.out = FALSE),
    regexp = "P-value Matrix|P-Value Matrix|p-value"
  )
})

test_that("aimTrace.qtlAim: plot = 'lrt' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  obj    <- make_mock_qtlAim(n_qtl = 2)
  result <- suppressWarnings(aimTrace(obj, plot = "lrt", lik.out = FALSE))
  expect_s3_class(result, "ggplot")
})

test_that("aimTrace.qtlAim: plot = 'stability' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  obj    <- make_mock_qtlAim(n_qtl = 2)
  result <- suppressWarnings(aimTrace(obj, plot = "stability", lik.out = FALSE))
  expect_s3_class(result, "ggplot")
})

test_that("aimTrace.qtlAim: plot = 'both' returns a named list with $lrt and $stability", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  obj    <- make_mock_qtlAim(n_qtl = 2)
  result <- suppressWarnings(aimTrace(obj, plot = "both", lik.out = FALSE))

  expect_type(result, "list")
  expect_named(result, c("lrt", "stability"))
  expect_s3_class(result$lrt,       "ggplot")
  expect_s3_class(result$stability, "ggplot")
})

test_that("aimTrace.qtlAim: iter beyond range produces warning and drops invalid iter", {
  obj    <- make_mock_qtlAim(n_qtl = 2)
  # valid iterations are 1:2; pass 99 as well
  expect_warning(
    aimTrace(obj, iter = c(1, 99), plot = FALSE, lik.out = FALSE),
    regexp = "iter"
  )
})

test_that("aimTrace.qtlAim: iter = all iterations (default) produces output without error", {
  obj <- make_mock_qtlAim(n_qtl = 3)
  expect_no_error(
    capture.output(aimTrace(obj, plot = FALSE, lik.out = FALSE))
  )
})

# =============================================================================
# aimTrace.gwasAim
# =============================================================================

test_that("aimTrace.gwasAim: plot = FALSE returns invisible NULL", {
  obj    <- make_mock_gwasAim(n_qtl = 2)
  result <- aimTrace(obj, plot = FALSE, lik.out = FALSE)
  expect_null(result)
})

test_that("aimTrace.gwasAim: output header contains 'GWAS'", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  expect_output(
    aimTrace(obj, plot = FALSE, lik.out = FALSE),
    regexp = "GWAS"
  )
})

test_that("aimTrace.gwasAim: lik.out = TRUE prints LRT table", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  expect_output(
    aimTrace(obj, plot = FALSE, lik.out = TRUE),
    regexp = "Likelihood Ratio"
  )
})

test_that("aimTrace.gwasAim: plot = 'lrt' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  obj    <- make_mock_gwasAim(n_qtl = 2)
  result <- suppressWarnings(aimTrace(obj, plot = "lrt", lik.out = FALSE))
  expect_s3_class(result, "ggplot")
})

test_that("aimTrace.gwasAim: plot = 'stability' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  obj    <- make_mock_gwasAim(n_qtl = 2)
  result <- suppressWarnings(aimTrace(obj, plot = "stability", lik.out = FALSE))
  expect_s3_class(result, "ggplot")
})

test_that("aimTrace.gwasAim: plot = 'both' returns list with $lrt and $stability", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  obj    <- make_mock_gwasAim(n_qtl = 2)
  result <- suppressWarnings(aimTrace(obj, plot = "both", lik.out = FALSE))

  expect_type(result, "list")
  expect_named(result, c("lrt", "stability"))
  expect_s3_class(result$lrt,       "ggplot")
  expect_s3_class(result$stability, "ggplot")
})

test_that("aimTrace.gwasAim: bad iter produces a warning", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  expect_warning(
    aimTrace(obj, iter = c(1, 999), plot = FALSE, lik.out = FALSE),
    regexp = "iter"
  )
})

test_that("aimTrace.gwasAim: incremental marker p-value matrix is always printed", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  expect_output(
    aimTrace(obj, plot = FALSE, lik.out = FALSE),
    regexp = "Marker P-value Matrix|P-value Matrix"
  )
})

test_that("aimTrace.gwasAim: TypeI and n.markers appear in console output", {
  obj <- make_mock_gwasAim(n_qtl = 2)
  out <- capture.output(aimTrace(obj, plot = FALSE, lik.out = FALSE))

  expect_true(any(grepl("TypeI",   out, ignore.case = TRUE)))
  expect_true(any(grepl("marker",  out, ignore.case = TRUE)))
})
