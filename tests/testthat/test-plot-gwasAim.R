# =============================================================================
# test-plot-gwasAim.R
# Tests for plot.gwasAim().
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_mock_gwasAim()       — fully populated gwasAim mock object
#   make_wgPanel()            — correct genObj class for gwasAim
#   make_wgCross_interval()   — wrong genObj class (triggers error guard)
#
# Strategy
# --------
# Type group A — "manhattan", "outlier", "blups", "heatmap"
#   Pure diagnostic types reading only $QTL$diag.  No ASReml calls.
#   Tested directly.
#
# Type group B — "effects", "contrast"
#   "effects":  calls summary.gwasAim()          (in wgAim namespace → mock)
#   "contrast": calls summary.gwasAim() AND predict(object, ...)
#               predict.asreml is NOT in wgAim namespace, so we mock the
#               entire .build_contrast_df() helper instead.
#
# Error-guard tests are pure input-validation.
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: build a structurally correct summary.gwasAim result.
#         Marker-type columns: Chromosome | Marker | dist(cM) | Size |
#         Pvalue | Perc.Var [| LOD]
# ---------------------------------------------------------------------------
.make_marker_summ_gwas <- function(gwas_fit, genObj) {
  eff_names <- names(gwas_fit$QTL$effects)
  parts     <- strsplit(sub("^X\\.", "", eff_names), "\\.")
  eff_chrs  <- sapply(parts, "[", 1)
  eff_idxs  <- as.integer(sapply(parts, "[", 2))
  n_qtl     <- length(eff_names)

  mk_pos  <- mapply(function(ch, idx) as.numeric(genObj$geno[[ch]]$map[idx]),
                    eff_chrs, eff_idxs)
  mk_nms  <- mapply(function(ch, idx) names(genObj$geno[[ch]]$map)[idx],
                    eff_chrs, eff_idxs)

  data.frame(
    Chromosome  = eff_chrs,
    Marker      = mk_nms,
    "dist(cM)"  = round(mk_pos, 2),
    Size        = round(gwas_fit$QTL$effects, 4),
    Pvalue      = 0.001,
    "Perc.Var"  = c(6.1, 3.8, 2.5)[seq_len(n_qtl)],
    LOD         = 3.0,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Helper: minimal contrast data-frame list for .plot_contrast() — GWAS path.
#         allele_class uses discrete "0", "1", "2" levels.
# ---------------------------------------------------------------------------
.make_contrast_list_gwas <- function(gwas_fit) {
  n_qtl   <- length(gwas_fit$QTL$effects)
  n_lines <- nrow(gwas_fit$mf)
  ids     <- as.character(gwas_fit$mf$id)

  allele_cols <- c("0" = "steelblue", "1" = "grey60", "2" = "firebrick")

  lapply(seq_len(n_qtl), function(k) {
    scores <- sample(c(-1L, 0L, 1L), n_lines, replace = TRUE)
    cls    <- as.character(scores + 1L)   # -1 -> "0", 0 -> "1", 1 -> "2"
    df <- data.frame(
      line          = ids,
      total_genetic = rnorm(n_lines),
      score         = as.numeric(scores),
      allele_class  = factor(cls, levels = names(allele_cols)),
      facet_label   = factor(
        paste0("Chr P", k, "  \u00b7  ", k * 5, " cM"),
        levels = paste0("Chr P", k, "  \u00b7  ", k * 5, " cM")),
      effect_txt    = "Effect: +0.3 \u00b1 0.1  (6%%  var)   p = 0.001",
      stringsAsFactors = FALSE
    )
    df$allele_cols_map <- list(allele_cols)
    df
  })
}


# =============================================================================
# Group A: Diagnostic plot types
# =============================================================================

test_that("plot.gwasAim type='manhattan' returns a ggplot", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "manhattan")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim type='outlier' returns a ggplot", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "outlier")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim type='blups' returns a ggplot", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "blups")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim type='heatmap' returns a ggplot", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "heatmap")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim type='heatmap' with cap= returns a ggplot", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "heatmap", cap = 4.0)

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# Group B: Results plot types — effects and contrast
# =============================================================================

test_that("plot.gwasAim type='effects' returns a ggplot (summary.gwasAim mocked)", {
  gwas_fit  <- make_mock_gwasAim()
  genObj    <- attr(gwas_fit, "genObj")

  mock_summ <- .make_marker_summ_gwas(gwas_fit, genObj)

  with_mocked_bindings(
    `summary.gwasAim` = function(object, genObj, ...) mock_summ,
    .package = "wgAim",
    {
      gp <- plot.gwasAim(gwas_fit, genObj, type = "effects")
      expect_s3_class(gp, "ggplot")
    }
  )
})

test_that("plot.gwasAim type='contrast' with data=NULL errors", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  expect_error(
    plot.gwasAim(gwas_fit, genObj, type = "contrast", data = NULL),
    "data is required"
  )
})

test_that("plot.gwasAim type='contrast' with phenoData returns ggplot (.build_contrast_df mocked)", {
  gwas_fit  <- make_mock_gwasAim()
  genObj    <- attr(gwas_fit, "genObj")
  phenoData <- gwas_fit$mf

  cdf_mock  <- .make_contrast_list_gwas(gwas_fit)
  attr(cdf_mock, "is_gwas") <- TRUE

  with_mocked_bindings(
    `.build_contrast_df` = function(object, genObj, data) cdf_mock,
    .package = "wgAim",
    {
      gp <- plot.gwasAim(gwas_fit, genObj, type = "contrast", data = phenoData)
      expect_s3_class(gp, "ggplot")
    }
  )
})

# =============================================================================
# Iteration and chromosome filters
# =============================================================================

test_that("plot.gwasAim chr= filter: single chromosome returns a ggplot (outlier)", {
  gwas_fit  <- make_mock_gwasAim()
  genObj    <- attr(gwas_fit, "genObj")
  first_chr <- names(genObj$geno)[1]

  gp <- plot.gwasAim(gwas_fit, genObj, type = "outlier", chr = first_chr)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim chr= filter: single chromosome returns a ggplot (manhattan)", {
  gwas_fit  <- make_mock_gwasAim()
  genObj    <- attr(gwas_fit, "genObj")
  first_chr <- names(genObj$geno)[1]

  gp <- plot.gwasAim(gwas_fit, genObj, type = "manhattan", chr = first_chr)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim iter= filter: iter=1 returns a ggplot", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "outlier", iter = 1L)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim iter= beyond range stops with error", {
  gwas_fit  <- make_mock_gwasAim()
  genObj    <- attr(gwas_fit, "genObj")
  n_stored  <- length(gwas_fit$QTL$diag$oint)

  expect_error(
    plot.gwasAim(gwas_fit, genObj, type = "outlier", iter = n_stored + 10L),
    "iter contains values greater than"
  )
})

# =============================================================================
# Input validation / error guards
# =============================================================================

test_that("plot.gwasAim errors when genObj is missing", {
  gwas_fit <- make_mock_gwasAim()

  expect_error(
    plot.gwasAim(gwas_fit),
    "genObj is a required argument"
  )
})

test_that("plot.gwasAim errors when genObj is wrong class (wgCross not wgPanel)", {
  gwas_fit <- make_mock_gwasAim()
  wrong    <- make_wgCross_interval()

  expect_error(
    plot.gwasAim(gwas_fit, wrong),
    "\"wgPanel\""
  )
})

test_that("plot.gwasAim errors when no significant markers (effects type)", {
  gwas_fit             <- make_mock_gwasAim()
  gwas_fit$QTL$effects <- NULL
  genObj               <- attr(gwas_fit, "genObj")

  expect_error(
    plot.gwasAim(gwas_fit, genObj, type = "effects"),
    "[Nn]o significant"
  )
})

test_that("plot.gwasAim errors when no significant markers (contrast type)", {
  gwas_fit             <- make_mock_gwasAim()
  gwas_fit$QTL$effects <- NULL
  genObj               <- attr(gwas_fit, "genObj")

  expect_error(
    plot.gwasAim(gwas_fit, genObj, type = "contrast"),
    "[Nn]o significant"
  )
})

test_that("plot.gwasAim errors for unknown chromosome name", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  expect_error(
    plot.gwasAim(gwas_fit, genObj, chr = "GHOST_CHR_999"),
    "[Nn]ot found in genObj"
  )
})

# =============================================================================
# Aesthetic arguments
# =============================================================================

test_that("plot.gwasAim chr.lines=TRUE returns a ggplot (outlier type)", {
  gwas_fit <- make_mock_gwasAim(n_chr = 2)
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(gwas_fit, genObj, type = "outlier", chr.lines = TRUE)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim custom sig.col and pt.col accepted (manhattan)", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(
    gwas_fit, genObj,
    type    = "manhattan",
    sig.col = "purple",
    pt.col  = c("navy", "goldenrod")
  )

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gwasAim band.col argument accepted (manhattan)", {
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")

  gp <- plot.gwasAim(
    gwas_fit, genObj,
    type     = "manhattan",
    band.col = c("lightyellow", "white")
  )

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# Internal: .plot_gwas_manhattan directly
# =============================================================================

test_that(".plot_gwas_manhattan produces a ggplot with a facet specification", {
  build_cumpos         <- getFromNamespace(".build_cumpos",         "wgAim")
  plot_gwas_manhattan  <- getFromNamespace(".plot_gwas_manhattan",  "wgAim")
  gwas_fit <- make_mock_gwasAim()
  genObj   <- attr(gwas_fit, "genObj")
  chrs     <- names(genObj$geno)
  iter     <- seq_len(length(gwas_fit$QTL$diag$oint))
  cp       <- build_cumpos(genObj, "marker", chrs)

  gp <- plot_gwas_manhattan(
    gwas_fit, genObj, iter, chrs, cp,
    "red", c("steelblue", "grey60"), c("grey92", "white"), 0.6
  )

  expect_s3_class(gp, "ggplot")
  # Must use facet_wrap (not FacetNull)
  expect_false(inherits(gp$facet, "FacetNull"))
})
