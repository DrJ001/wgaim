# =============================================================================
# test-plot-qtlAim.R
# Tests for plot.qtlAim() and shared internal engine helpers.
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_mock_qtlAim()         — fully populated qtlAim mock object
#   make_wgCross_interval()    — minimal wgCross (interval type)
#   make_wgCross_marker()      — minimal wgCross (marker type)
#   make_wgPanel()             — wrong class for qtlAim (triggers error)
#
# Strategy
# --------
# Type group A — pure diagnostic types: "outlier", "blups", "chr", "heatmap"
#   Only read $QTL$diag slot; no ASReml calls needed.
#   Tested directly with the mock object.
#
# Type group B — results types: "effects", "contrast"
#   .build_effects_df()  calls summary(object, genObj)     [summary.qtlAim]
#   .build_contrast_df() calls summary(object, genObj) AND predict(object, ...)
#
#   For "effects": mock summary.qtlAim (in wgAim namespace) with
#   with_mocked_bindings().
#
#   For "contrast": additionally need predict.asreml — not in wgAim namespace,
#   so we mock the entire .build_contrast_df() function directly, returning a
#   minimal per-QTL data-frame list that .plot_contrast() can consume.
#
# Error-guard tests are pure input-validation — no full render required.
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: build a structurally correct summary.qtlAim result for n_qtl QTLs
#         drawn from a given qtlAim mock and its genObj.
#         (interval type columns: Chr | LMark | dist | IMark | dist | RMark |
#          dist | Size | Prob | Perc.Var [ | LOD ])
# ---------------------------------------------------------------------------
.make_interval_summ <- function(qtl_fit, genObj) {
  eff_names <- names(qtl_fit$QTL$effects)
  parts     <- strsplit(sub("^X\\.", "", eff_names), "\\.")
  eff_chrs  <- sapply(parts, "[", 1)
  eff_idxs  <- as.integer(sapply(parts, "[", 2))
  n_qtl     <- length(eff_names)

  imarks  <- mapply(function(ch, idx) genObj$geno[[ch]]$inferred.map[idx],
                    eff_chrs, eff_idxs, SIMPLIFY = TRUE)
  lh_pos  <- mapply(function(ch, idx) {
    pos <- genObj$geno[[ch]]$map
    pos[findInterval(imarks[match(paste(ch, idx, sep="."),
                                 paste(eff_chrs, eff_idxs, sep="."))],
                     pos)]
  }, eff_chrs, eff_idxs, SIMPLIFY = TRUE)

  # Simpler: just use map[idx] as the left flank position
  lh_marks <- mapply(function(ch, idx) as.numeric(genObj$geno[[ch]]$map[idx]),
                     eff_chrs, eff_idxs)
  rh_marks <- lh_marks + diff(range(genObj$geno[[eff_chrs[1]]]$map)) / 5

  data.frame(
    Chromosome       = eff_chrs,
    "Left Marker"    = paste0("mk_L", seq_len(n_qtl)),
    "dist(cM).lh"    = round(lh_marks, 2),
    "Infer. Marker"  = paste0("int", seq_len(n_qtl)),
    "dist(cM).im"    = round(as.numeric(imarks), 2),
    "Right Marker"   = paste0("mk_R", seq_len(n_qtl)),
    "dist(cM).rh"    = round(rh_marks, 2),
    Size             = round(qtl_fit$QTL$effects, 4),
    Prob             = 0.01,
    "Perc.Var"       = c(5.2, 4.1, 3.0)[seq_len(n_qtl)],
    LOD              = 2.5,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Helper: build a minimal contrast data-frame list for .plot_contrast().
#         is_gwas = FALSE (biparental path: allele_class "B (-1)" / "A (+1)").
# ---------------------------------------------------------------------------
.make_contrast_list <- function(qtl_fit, genObj) {
  n_qtl   <- length(qtl_fit$QTL$effects)
  n_lines <- nrow(qtl_fit$mf)
  ids     <- as.character(qtl_fit$mf$id)

  allele_cols <- c("B (\u22121)" = "steelblue", "A (+1)" = "firebrick")

  lapply(seq_len(n_qtl), function(k) {
    scores <- sample(c(-1, 1), n_lines, replace = TRUE)
    df <- data.frame(
      line          = ids,
      total_genetic = rnorm(n_lines),
      score         = scores,
      allele_class  = factor(
        ifelse(scores < 0, "B (\u22121)", "A (+1)"),
        levels = names(allele_cols)),
      facet_label   = factor(
        paste0("Chr C", k, "  \u00b7  ", k * 10, " cM"),
        levels = paste0("Chr C", k, "  \u00b7  ", k * 10, " cM")),
      effect_txt    = sprintf("Effect: +0.4 \u00b1 0.1  (5%%  var)   p = 0.01"),
      stringsAsFactors = FALSE
    )
    df$allele_cols_map <- list(allele_cols)
    df
  })
}


# =============================================================================
# Group A: Diagnostic plot types
# =============================================================================

test_that("plot.qtlAim type='outlier' returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "outlier")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim type='blups' returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "blups")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim type='chr' returns a ggplot (ochr populated)", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  # Fixture always builds ochr
  expect_false(is.null(qtl_fit$QTL$diag$ochr[[1]]))

  gp <- plot.qtlAim(qtl_fit, genObj, type = "chr")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim type='chr' errors when ochr is NULL/empty", {
  qtl_fit <- make_mock_qtlAim()
  # Blank out all ochr entries to NULL
  qtl_fit$QTL$diag$ochr <- vector("list", length(qtl_fit$QTL$diag$ochr))
  genObj  <- attr(qtl_fit, "genObj")

  expect_error(
    plot.qtlAim(qtl_fit, genObj, type = "chr"),
    "[Cc]hromosome outlier statistics"
  )
})

test_that("plot.qtlAim type='heatmap' returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "heatmap")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim type='heatmap' with cap= argument returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "heatmap", cap = 3.0)

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# Group B: Results plot types — effects and contrast
# =============================================================================

test_that("plot.qtlAim type='effects' returns a ggplot (summary.qtlAim mocked)", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  mock_summ <- .make_interval_summ(qtl_fit, genObj)

  with_mocked_bindings(
    `summary.qtlAim` = function(object, genObj, ...) mock_summ,
    .package = "wgAim",
    {
      gp <- plot.qtlAim(qtl_fit, genObj, type = "effects")
      expect_s3_class(gp, "ggplot")
    }
  )
})

test_that("plot.qtlAim type='contrast' with data=NULL stops with error", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  expect_error(
    plot.qtlAim(qtl_fit, genObj, type = "contrast", data = NULL),
    "data is required"
  )
})

test_that("plot.qtlAim type='contrast' with phenoData returns a ggplot (.build_contrast_df mocked)", {
  qtl_fit   <- make_mock_qtlAim()
  genObj    <- attr(qtl_fit, "genObj")
  phenoData <- qtl_fit$mf

  cdf_mock <- .make_contrast_list(qtl_fit, genObj)
  attr(cdf_mock, "is_gwas") <- FALSE

  with_mocked_bindings(
    `.build_contrast_df` = function(object, genObj, data) cdf_mock,
    .package = "wgAim",
    {
      gp <- plot.qtlAim(qtl_fit, genObj, type = "contrast", data = phenoData)
      expect_s3_class(gp, "ggplot")
    }
  )
})

# =============================================================================
# Iteration and chromosome filters
# =============================================================================

test_that("plot.qtlAim chr= filter: single chromosome returns a ggplot", {
  qtl_fit   <- make_mock_qtlAim()
  genObj    <- attr(qtl_fit, "genObj")
  first_chr <- names(genObj$geno)[1]

  gp <- plot.qtlAim(qtl_fit, genObj, type = "outlier", chr = first_chr)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim iter= filter: iter=1 returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "outlier", iter = 1L)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim iter= filter: multiple iterations returns a ggplot", {
  qtl_fit <- make_mock_qtlAim(n_qtl = 3)
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "blups", iter = c(1L, 2L))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim iter= beyond range stops with informative error", {
  qtl_fit  <- make_mock_qtlAim()
  genObj   <- attr(qtl_fit, "genObj")
  n_stored <- length(qtl_fit$QTL$diag$oint)

  expect_error(
    plot.qtlAim(qtl_fit, genObj, type = "outlier", iter = n_stored + 5L),
    "iter contains values greater than"
  )
})

# =============================================================================
# Input validation / error guards
# =============================================================================

test_that("plot.qtlAim errors when genObj is missing", {
  qtl_fit <- make_mock_qtlAim()

  expect_error(
    plot.qtlAim(qtl_fit),
    "genObj is a required argument"
  )
})

test_that("plot.qtlAim errors when genObj is wrong class (plain list, not wgCross)", {
  qtl_fit <- make_mock_qtlAim()
  # wgPanel inherits wgCross, so we need something that does NOT inherit wgCross
  wrong   <- list(geno = list(), pheno = data.frame())

  expect_error(
    plot.qtlAim(qtl_fit, wrong),
    "\"wgCross\""
  )
})

test_that("plot.qtlAim errors when QTL$effects is NULL", {
  qtl_fit             <- make_mock_qtlAim()
  qtl_fit$QTL$effects <- NULL
  genObj              <- attr(qtl_fit, "genObj")

  expect_error(
    plot.qtlAim(qtl_fit, genObj),
    "[Nn]o significant"
  )
})

test_that("plot.qtlAim errors for unknown chromosome name", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  expect_error(
    plot.qtlAim(qtl_fit, genObj, chr = "GHOST_CHR_XYZ"),
    "[Nn]ot found in genObj"
  )
})

# =============================================================================
# Aesthetic arguments
# =============================================================================

test_that("plot.qtlAim chr.lines=TRUE returns a ggplot", {
  qtl_fit <- make_mock_qtlAim(n_chr = 3)
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "outlier", chr.lines = TRUE)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim sig.col='blue' is accepted and returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "outlier", sig.col = "blue")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.qtlAim custom pt.col vector returns a ggplot", {
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")

  gp <- plot.qtlAim(qtl_fit, genObj, type = "blups",
                    pt.col = c("navy", "goldenrod", "forestgreen"))

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# Marker-type wgCross (non-interval)
# =============================================================================

test_that("plot.qtlAim works with marker-type qtlAim (no inferred.map)", {
  set.seed(99)
  genObj    <- make_wgCross_marker(n_chr = 2, n_mar = 6)
  chrs      <- names(genObj$geno)
  n_lines   <- nrow(genObj$pheno)
  ids       <- paste0("L", seq_len(n_lines))
  n_qtl     <- 2L

  chr_idx   <- chrs[c(1, 2)]
  int_idx   <- c(2L, 3L)
  qtl_keys  <- paste("Chr", chr_idx, int_idx, sep = ".")
  eff_names <- paste("X",   chr_idx, int_idx, sep = ".")
  effects   <- setNames(rnorm(n_qtl, 0.4, 0.1), eff_names)
  veffects  <- setNames(runif(n_qtl, 0.02, 0.05), eff_names)

  n_total  <- 2 * 6
  all_keys <- paste("Chr", rep(chrs, each = 6),
                    rep(seq_len(6), times = 2), sep = ".")
  oint_list <- lapply(seq_len(n_qtl), function(k) {
    v <- runif(n_total, 0, 2); names(v) <- all_keys; v
  })
  blup_list <- lapply(seq_len(n_qtl), function(k) {
    v <- rnorm(n_total); names(v) <- all_keys; v
  })

  QTL <- list(
    qtl        = qtl_keys,
    effects    = effects,
    veffects   = veffects,
    method     = "fixed",
    type       = "marker",
    selection  = "interval",
    TypeI      = 0.05,
    iterations = n_qtl + 1L,
    breakout   = FALSE,
    diag       = list(
      oint         = oint_list,
      blups        = blup_list,
      lik          = lapply(seq_len(n_qtl), function(k)
        list(baseLogL = -50, stat = 5, pvalue = 0.02, pass = TRUE)),
      ochr         = NULL,
      state        = setNames(rep(1, n_total), all_keys),
      genetic.term = "id",
      rel.scale    = 1
    )
  )

  obj <- list(
    converge        = TRUE,
    loglik          = -45,
    sigma2          = 1.0,
    vparameters     = c(id = 0.3, "R!variance" = 1.0),
    vparameters.con = c(id = 0L, "R!variance" = 0L),
    coefficients    = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
    call     = call("asreml"),
    QTL      = QTL
  )
  class(obj) <- c("qtlAim", "asreml")

  gp <- plot.qtlAim(obj, genObj, type = "outlier")
  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# Internal helpers accessed via :::
# =============================================================================

test_that(".build_cumpos (interval) returns correct named slots", {
  build_cumpos <- getFromNamespace(".build_cumpos", "wgAim")
  genObj <- make_wgCross_interval(n_chr = 2, n_mar = 4)
  chrs   <- names(genObj$geno)

  cp <- build_cumpos(genObj, "interval", chrs)

  expect_named(cp, c("pos_lookup", "chr.mid", "chr.end"))
  expect_type(cp$pos_lookup, "double")
  expect_length(cp$chr.mid, length(chrs))
  expect_length(cp$chr.end, length(chrs))
  # All keys start with "Chr."
  expect_true(all(grepl("^Chr\\.", names(cp$pos_lookup))))
})

test_that(".build_cumpos (marker) uses map (not inferred.map)", {
  build_cumpos <- getFromNamespace(".build_cumpos", "wgAim")
  genObj <- make_wgCross_marker(n_chr = 3, n_mar = 5)
  chrs   <- names(genObj$geno)

  cp <- build_cumpos(genObj, "marker", chrs)

  expect_type(cp$pos_lookup, "double")
  expect_length(cp$chr.mid, 3)
  # Cumulative positions are non-decreasing
  expect_true(all(diff(sort(cp$pos_lookup)) >= 0))
})

test_that(".build_stat_df returns long-format data frame with required columns", {
  build_cumpos  <- getFromNamespace(".build_cumpos",  "wgAim")
  build_stat_df <- getFromNamespace(".build_stat_df", "wgAim")
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")
  chrs    <- names(genObj$geno)
  cp      <- build_cumpos(genObj, "interval", chrs)

  df <- build_stat_df(qtl_fit, "oint", 1L, chrs, cp)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("values", "chr", "dist", "iteration") %in% names(df)))
  expect_s3_class(df$iteration, "factor")
})

test_that(".build_stat_df errors when no matching values found", {
  build_cumpos  <- getFromNamespace(".build_cumpos",  "wgAim")
  build_stat_df <- getFromNamespace(".build_stat_df", "wgAim")
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")
  chrs    <- names(genObj$geno)
  cp      <- build_cumpos(genObj, "interval", chrs)

  expect_error(
    build_stat_df(qtl_fit, "oint", 1L, "GHOST_CHR", cp),
    "[Nn]o statistic values"
  )
})

test_that(".build_heatmap_df returns heat_df and sig_df", {
  build_cumpos     <- getFromNamespace(".build_cumpos",     "wgAim")
  build_heatmap_df <- getFromNamespace(".build_heatmap_df", "wgAim")
  qtl_fit <- make_mock_qtlAim()
  genObj  <- attr(qtl_fit, "genObj")
  chrs    <- names(genObj$geno)
  cp      <- build_cumpos(genObj, "interval", chrs)
  iters   <- seq_len(length(qtl_fit$QTL$diag$oint))

  hd <- build_heatmap_df(qtl_fit, iters, chrs, cp)

  expect_named(hd, c("heat_df", "sig_df"))
  expect_s3_class(hd$heat_df, "data.frame")
  expect_true(all(c("dist", "iter_num", "fill_val", "tile_w") %in% names(hd$heat_df)))
})
