# =============================================================================
# test-plot-gpAim.R
# Tests for plot.gpAim().
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_mock_gpAim()          — fully populated gpAim mock object
#   make_wgCross_interval()    — wgCross (correct class for blups type)
#   make_wgPanel()             — wgPanel (also acceptable for blups type)
#
# Strategy
# --------
# All four plot types work purely from $GP data — no ASReml calls.
# Tests exercise the three main branches:
#
#   caterpillar  — .plot_gp_caterpillar()
#   heatmap      — .plot_gp_heatmap()  (warns for n > 500)
#   density      — .plot_gp_density()
#   blups        — .plot_gp_blups()   (requires genObj; reads marker.effects)
#
# The fixture's $GP$marker.effects is a named numeric vector, which does NOT
# match the real gpAim output (a data.frame with columns "marker" and
# "effect").  Tests for type = "blups" therefore build a local gpAim mock with
# the correct data.frame structure.
#
# For the large-n heatmap warning test we build a gpAim mock with n_lines=600.
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: build a gpAim mock whose GP$marker.effects is a properly structured
# data.frame with columns "marker" and "effect".  Used for blups tests.
# ---------------------------------------------------------------------------
.make_gpAim_for_blups <- function(n_lines = 40, n_chr = 2, n_mar = 8) {
  set.seed(21)
  genObj <- make_wgCross_interval(n_lines, n_chr, n_mar)
  chrs   <- names(genObj$geno)
  ids    <- paste0("L", seq_len(n_lines))

  n_markers <- n_chr * n_mar

  # All marker keys in "Chr.CHRNAME.IDX" format — same format .build_cumpos uses
  all_keys <- paste("Chr", rep(chrs, each = n_mar),
                    rep(seq_len(n_mar), times = n_chr), sep = ".")

  marker_effects_df <- data.frame(
    marker = all_keys,
    effect = rnorm(n_markers),
    stringsAsFactors = FALSE
  )

    SE_vals <- runif(n_lines, 0.1, 0.3)
    gebv_df <- data.frame(
      id       = ids,
      GEBV     = rnorm(n_lines, 0, 1),
      SE       = SE_vals,
      Accuracy = sqrt(pmax(0, 1 - SE_vals^2 / 0.6)),
      gen.H2   = 0.75,
      stringsAsFactors = FALSE
    )
    raw_mat <- matrix(rnorm(n_lines * n_markers), nrow = n_lines)
    rel_mat <- tcrossprod(raw_mat) / n_markers
    rownames(rel_mat) <- colnames(rel_mat) <- ids

    GP <- list(
      gebv           = gebv_df,
      marker.effects = marker_effects_df,
      gen.type       = "interval",
      path           = "vm",
      var.genetic    = 0.6,
      var.resid      = 0.4,
      heritability   = 0.6,
      gen.H2         = 0.75,
      n.markers      = n_markers,
      rel.scale      = n_markers,
      rel.matrix     = rel_mat,
      genetic.term   = "id"
    )

  obj <- list(
    converge        = TRUE,
    loglik          = -80,
    sigma2          = 1.0,
    vparameters     = c(id = 0.6, "R!variance" = 1.0),
    vparameters.con = c(id = 0L,  "R!variance" = 0L),
    coefficients    = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
    call     = call("asreml"),
    GP       = GP
  )
  class(obj) <- c("gpAim", "asreml")
  attr(obj, "genObj") <- genObj
  obj
}

# ---------------------------------------------------------------------------
# Helper: build a large gpAim mock (n_lines > 500) for the heatmap warning.
# ---------------------------------------------------------------------------
.make_gpAim_large <- function(n_lines = 600) {
  set.seed(77)
  ids     <- paste0("L", seq_len(n_lines))
  n_marks <- 10L
  rel_mat <- diag(n_lines)  # identity G for speed
  rownames(rel_mat) <- colnames(rel_mat) <- ids

  se_large <- runif(n_lines, 0.1, 0.3)
  GP <- list(
    gebv = data.frame(
      id       = ids,
      GEBV     = rnorm(n_lines),
      SE       = se_large,
      Accuracy = sqrt(pmax(0, 1 - se_large^2 / 0.5)),
      gen.H2   = 0.70,
      stringsAsFactors = FALSE
    ),
    marker.effects = data.frame(
      marker = paste0("m", seq_len(n_marks)),
      effect = rnorm(n_marks),
      stringsAsFactors = FALSE
    ),
    gen.type     = "marker",
    path         = "vm",
    var.genetic  = 0.5,
    var.resid    = 0.5,
    heritability = 0.5,
    gen.H2       = 0.70,
    n.markers    = n_marks,
    rel.scale    = n_marks,
    rel.matrix   = rel_mat,
    genetic.term = "id"
  )

  obj <- list(
    converge        = TRUE,
    loglik          = -200,
    sigma2          = 1.0,
    vparameters     = c(id = 0.5, "R!variance" = 1.0),
    vparameters.con = c(id = 0L,  "R!variance" = 0L),
    coefficients    = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
    call     = call("asreml"),
    GP       = GP
  )
  class(obj) <- c("gpAim", "asreml")
  obj
}

# =============================================================================
# type = "caterpillar"
# =============================================================================

test_that("plot.gpAim type='caterpillar' returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  # suppressWarnings: aes_string() deprecation in package source (not tests)
  gp <- suppressWarnings(plot.gpAim(gp_fit, type = "caterpillar"))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim caterpillar: prop.select=0.1 returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  gp <- suppressWarnings(plot.gpAim(gp_fit, type = "caterpillar", prop.select = 0.1))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim caterpillar: prop.select=0.3 returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  gp <- suppressWarnings(plot.gpAim(gp_fit, type = "caterpillar", prop.select = 0.3))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim caterpillar: threshold as raw GEBV value returns a ggplot", {
  gp_fit <- make_mock_gpAim()
  # Use median GEBV as a raw threshold (outside 0-1 quantile interpretation
  # only when the value itself is outside (0,1); ensure we use a value > 1)
  raw_thr <- max(gp_fit$GP$gebv$GEBV) - 0.5

  gp <- suppressWarnings(plot.gpAim(gp_fit, type = "caterpillar", threshold = raw_thr))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim caterpillar: threshold=0.9 treated as 90th-percentile quantile", {
  gp_fit <- make_mock_gpAim()

  # threshold in (0,1) is interpreted as a quantile — should succeed
  gp <- suppressWarnings(plot.gpAim(gp_fit, type = "caterpillar", threshold = 0.9))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim caterpillar: top.n capped to floor(n/2) when top.n > n/2", {
  gp_fit <- make_mock_gpAim(n_lines = 10)
  n      <- nrow(gp_fit$GP$gebv)

  # top.n = n (> n/2 = 5) — should be silently capped to 5, not error.
  # suppressWarnings: aes_string() deprecation in package source (not tests)
  expect_no_error(
    suppressWarnings(plot.gpAim(gp_fit, type = "caterpillar", top.n = n))
  )
})

# =============================================================================
# type = "heatmap"
# =============================================================================

test_that("plot.gpAim type='heatmap' returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  # suppressWarnings: aes_string() deprecation in package source (not tests)
  gp <- suppressWarnings(plot.gpAim(gp_fit, type = "heatmap"))

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim type='heatmap' with n > 500 lines issues a warning about n", {
  large_fit <- .make_gpAim_large(n_lines = 600)

  # Capture all warnings — do NOT wrap in suppressWarnings() or the handler
  # will not see them.  At least one must mention the line count (600).
  # The aes_string() deprecation warning may also be raised; we only care
  # that the large-n warning is present.
  warns <- character(0)
  withCallingHandlers(
    plot.gpAim(large_fit, type = "heatmap"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("600", warns)))
})

test_that("plot.gpAim heatmap: custom pt.col returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  gp <- suppressWarnings(
    plot.gpAim(gp_fit, type = "heatmap", pt.col = c("darkred", "lightblue"))
  )

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# type = "density"
# =============================================================================

test_that("plot.gpAim type='density' returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  gp <- plot.gpAim(gp_fit, type = "density")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim density: prop.select=0.1 returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  gp <- plot.gpAim(gp_fit, type = "density", prop.select = 0.1)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim density: threshold=0.5 (quantile) returns a ggplot", {
  gp_fit <- make_mock_gpAim()

  gp <- plot.gpAim(gp_fit, type = "density", threshold = 0.5)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim density: threshold as raw value outside (0,1) returns a ggplot", {
  gp_fit  <- make_mock_gpAim()
  raw_thr <- max(gp_fit$GP$gebv$GEBV) - 0.5  # guaranteed > 0, could be > 1

  gp <- plot.gpAim(gp_fit, type = "density", threshold = raw_thr)

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# type = "blups"
# =============================================================================

test_that("plot.gpAim type='blups' with genObj=NULL stops with error", {
  gp_fit <- make_mock_gpAim()

  expect_error(
    plot.gpAim(gp_fit, genObj = NULL, type = "blups"),
    "genObj is required"
  )
})

test_that("plot.gpAim type='blups' with wrong genObj class stops with error", {
  gp_fit <- make_mock_gpAim()
  # Pass a plain list — not wgCross or wgPanel
  bad_gen <- list(geno = list())

  expect_error(
    plot.gpAim(gp_fit, genObj = bad_gen, type = "blups"),
    "\"wgCross\" or \"wgPanel\""
  )
})

test_that("plot.gpAim type='blups' with valid wgCross genObj returns a ggplot", {
  gp_fit <- .make_gpAim_for_blups(n_lines = 40, n_chr = 2, n_mar = 8)
  genObj <- attr(gp_fit, "genObj")

  gp <- plot.gpAim(gp_fit, genObj = genObj, type = "blups")

  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim type='blups' with wgPanel genObj returns a ggplot", {
  # Build a gpAim mock that uses wgPanel (marker gen.type)
  set.seed(31)
  n_lines  <- 30L
  n_chr    <- 2L
  n_mar    <- 6L
  genObj   <- make_wgPanel(n_lines, n_chr, n_mar)
  chrs     <- names(genObj$geno)
  ids      <- paste0("S", seq_len(n_lines))
  n_markers <- n_chr * n_mar

  all_keys <- paste("Chr", rep(chrs, each = n_mar),
                    rep(seq_len(n_mar), times = n_chr), sep = ".")

  marker_effects_df <- data.frame(
    marker = all_keys,
    effect = rnorm(n_markers),
    stringsAsFactors = FALSE
  )
  rel_mat <- diag(n_lines)
  rownames(rel_mat) <- colnames(rel_mat) <- ids

  se_wgp <- runif(n_lines, 0.1, 0.3)
  GP <- list(
    gebv           = data.frame(id = ids, GEBV = rnorm(n_lines),
                                SE  = se_wgp,
                                Accuracy = sqrt(pmax(0, 1 - se_wgp^2 / 0.5)),
                                gen.H2   = 0.68),
    marker.effects = marker_effects_df,
    gen.type       = "marker",
    path           = "vm",
    var.genetic    = 0.5,
    var.resid      = 0.5,
    heritability   = 0.5,
    gen.H2         = 0.68,
    n.markers      = n_markers,
    rel.scale      = n_markers,
    rel.matrix     = rel_mat,
    genetic.term   = "id"
  )

  obj <- list(
    converge        = TRUE,
    loglik          = -100,
    sigma2          = 1.0,
    vparameters     = c(id = 0.5, "R!variance" = 1.0),
    vparameters.con = c(id = 0L,  "R!variance" = 0L),
    coefficients    = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))
    ),
    formulae = list(fixed = yld ~ 1, random = ~ id),
    mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
    call     = call("asreml"),
    GP       = GP
  )
  class(obj) <- c("gpAim", "asreml")

  gp <- plot.gpAim(obj, genObj = genObj, type = "blups")
  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim type='blups' errors when marker.effects is NULL", {
  gp_fit <- .make_gpAim_for_blups()
  gp_fit$GP$marker.effects <- NULL
  genObj <- attr(gp_fit, "genObj")

  expect_error(
    plot.gpAim(gp_fit, genObj = genObj, type = "blups"),
    "[Nn]o marker effects"
  )
})

# =============================================================================
# .gp_threshold helper (internal)
# =============================================================================

test_that(".gp_threshold: value in (0,1) treated as quantile", {
  gp_threshold <- getFromNamespace(".gp_threshold", "wgAim")
  set.seed(1)
  vals <- rnorm(100)

  thr <- gp_threshold(vals, threshold = 0.9, prop.select = 0.1)

  expect_equal(thr, quantile(vals, 0.9, names = FALSE))
})

test_that(".gp_threshold: value outside (0,1) returned as-is", {
  gp_threshold <- getFromNamespace(".gp_threshold", "wgAim")
  vals <- rnorm(100)

  thr <- gp_threshold(vals, threshold = 2.5, prop.select = 0.1)

  expect_equal(thr, 2.5)
})

test_that(".gp_threshold: NULL threshold uses prop.select quantile", {
  gp_threshold <- getFromNamespace(".gp_threshold", "wgAim")
  set.seed(2)
  vals <- rnorm(100)

  thr <- gp_threshold(vals, threshold = NULL, prop.select = 0.2)

  expect_equal(thr, quantile(vals, 1 - 0.2, names = FALSE))
})

# =============================================================================
# MV plot types (caterpillar_mv and density_mv)
# Tests that the factor Trial column and gen.H2 per-trial work correctly.
# =============================================================================

.make_mock_gpAim_mv <- function(n_lines = 30, trials = c("T1", "T2")) {
  set.seed(55)
  ids    <- paste0("L", seq_len(n_lines))
  ntrait <- length(trials)
  SE_mv  <- runif(n_lines * ntrait, 0.05, 0.20)
  G_diag <- setNames(c(0.8, 0.5)[seq_len(ntrait)], trials)
  gen.H2_vals <- setNames(c(0.85, 0.78)[seq_len(ntrait)], trials)

  gebv_df <- data.frame(
    id       = factor(rep(ids, ntrait)),
    Trial    = factor(rep(trials, each = n_lines), levels = trials),
    GEBV     = rnorm(n_lines * ntrait),
    SE       = SE_mv,
    Accuracy = sqrt(pmax(0, 1 - SE_mv^2 /
                         rep(G_diag[trials], each = n_lines))),
    gen.H2   = rep(gen.H2_vals[trials], each = n_lines),
    stringsAsFactors = FALSE
  )
  names(gebv_df)[1] <- "id"

  n_markers <- 10L
  rel_mat   <- diag(n_lines)
  rownames(rel_mat) <- colnames(rel_mat) <- ids

  Ga   <- diag(G_diag)
  sds  <- sqrt(G_diag)
  Gcor <- Ga / outer(sds, sds); diag(Gcor) <- 1
  dimnames(Ga) <- dimnames(Gcor) <- list(trials, trials)

  GP <- list(
    gebv         = gebv_df,
    marker.effects = NULL,
    gen.type     = "marker",
    path         = "vm",
    var.genetic  = G_diag,
    var.resid    = setNames(rep(1, ntrait), trials),
    heritability = setNames(G_diag / (G_diag + 1), trials),
    gen.H2       = gen.H2_vals,
    n.markers    = n_markers,
    rel.scale    = n_markers,
    rel.matrix   = rel_mat,
    genetic.term = "id",
    Trait        = "Trial",
    trait.levels = trials,
    Ga           = Ga,
    Gcor         = Gcor
  )

  obj <- list(
    converge        = TRUE, loglik = -100, sigma2 = 1.0,
    vparameters     = c(id = 0.6, "R!variance" = 1.0),
    vparameters.con = c(id = 0L,  "R!variance" = 0L),
    coefficients    = list(
      fixed  = matrix(5, 1, 1, dimnames = list("(Intercept)", "effect")),
      random = matrix(rnorm(n_lines), n_lines, 1,
                      dimnames = list(paste0("id_", ids), "effect"))),
    formulae = list(fixed = yld ~ 1, random = ~ Trial:id),
    mf       = data.frame(id = factor(ids), yld = rnorm(n_lines)),
    call     = call("asreml"),
    GP       = GP
  )
  class(obj) <- c("gpAim", "asreml")
  obj
}

test_that("plot.gpAim MV caterpillar: Trial factor column produces faceted ggplot", {
  gp_fit <- .make_mock_gpAim_mv()
  gp <- plot.gpAim(gp_fit, type = "caterpillar")
  expect_s3_class(gp, "ggplot")
  # Should have a FacetWrap layer (one panel per trial)
  expect_true(inherits(gp$facet, "FacetWrap"))
})

test_that("plot.gpAim MV caterpillar: Trial levels preserved in factor column", {
  gp_fit <- .make_mock_gpAim_mv(trials = c("T1", "T2"))
  # Confirm the Trait column is a factor before plotting
  expect_s3_class(gp_fit$GP$gebv$Trial, "factor")
  expect_equal(levels(gp_fit$GP$gebv$Trial), c("T1", "T2"))
  gp <- plot.gpAim(gp_fit, type = "caterpillar")
  expect_s3_class(gp, "ggplot")
})

test_that("plot.gpAim MV density: returns faceted ggplot with gen.H2 annotation", {
  gp_fit <- .make_mock_gpAim_mv()
  gp <- plot.gpAim(gp_fit, type = "density", prop.select = 0.1)
  expect_s3_class(gp, "ggplot")
  expect_true(inherits(gp$facet, "FacetWrap"))
})
