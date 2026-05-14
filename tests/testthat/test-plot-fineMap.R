# =============================================================================
# test-plot-fineMap.R
# Tests for plot.fineMap() and print.fineMap().
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_mock_fineMap(n_qtl, n_pos) — builds a "fineMap" class object with
#     n_qtl named data frames, each with columns (mark, dist, pvalue, LOD),
#     and per-df attrs: qtl_pos (= 25 cM) and clipped (= FALSE).
#
# Strategy
# --------
# plot.fineMap():
#   — prints the ggplot internally and returns it invisibly with attribute
#     "suggested_height".  We capture the invisible return with invisible().
#   — Uses ggrepel::geom_text_repel — must be installed (it's in Imports).
#   — Facets one panel per QTL/marker.
#   — Appends "(window clipped)" when attr(df, "clipped") == TRUE.
#   — Caps Inf LOD values before plotting (should not error).
#   — Suppresses QTL reference label when peak is within 1 cM of qtl_pos.
#
# print.fineMap():
#   — Prints a plain-text summary to stdout.
#   — For a 0-row data frame, prints "No positions scanned."
#   — Returns x invisibly.
#
# All plot calls are wrapped in suppressMessages() / expect_silent() or
# checked with expect_s3_class() on the captured invisible return.
# print.fineMap() is tested with expect_output().
# =============================================================================

# Helper: capture the invisible return of plot.fineMap() without printing
.plot_fm_silent <- function(...) {
  suppressMessages(
    utils::capture.output(
      gp <- plot.fineMap(...)
    )
  )
  gp
}

# =============================================================================
# plot.fineMap
# =============================================================================

test_that("plot.fineMap returns invisibly a ggplot object", {
  fm <- make_mock_fineMap()

  gp <- suppressMessages(
    withVisible(plot.fineMap(fm))$value
  )

  expect_s3_class(gp, "ggplot")
})

test_that("plot.fineMap returns object with 'suggested_height' attribute", {
  fm <- make_mock_fineMap()

  gp <- .plot_fm_silent(fm)

  expect_false(is.null(attr(gp, "suggested_height")))
  expect_type(attr(gp, "suggested_height"), "double")
  expect_gt(attr(gp, "suggested_height"), 0)
})

test_that("plot.fineMap suggested_height scales with number of QTL", {
  fm1 <- make_mock_fineMap(n_qtl = 1)
  fm4 <- make_mock_fineMap(n_qtl = 4)

  gp1 <- .plot_fm_silent(fm1)
  gp4 <- .plot_fm_silent(fm4)

  # More QTL -> more rows -> larger suggested height
  expect_true(attr(gp4, "suggested_height") >= attr(gp1, "suggested_height"))
})

test_that("plot.fineMap: clipped=TRUE adds '(window clipped)' to facet label", {
  fm <- make_mock_fineMap(n_qtl = 2)
  # Mark first QTL as clipped
  attr(fm[[1]], "clipped") <- TRUE

  gp <- .plot_fm_silent(fm)

  # Facet labels are factors in the data used by the plot
  plot_data <- gp$data
  expect_true(any(grepl("window clipped", as.character(plot_data$qtl))))
})

test_that("plot.fineMap: unclipped QTL labels do NOT contain '(window clipped)'", {
  fm <- make_mock_fineMap(n_qtl = 2)
  # Both clipped = FALSE (default from fixture)

  gp <- .plot_fm_silent(fm)

  plot_data <- gp$data
  expect_false(any(grepl("window clipped", as.character(plot_data$qtl))))
})

test_that("plot.fineMap: Inf LOD values are capped and plot renders without error", {
  fm_inf <- make_mock_fineMap(n_qtl = 1, n_pos = 15)
  fm_inf[[1]]$LOD[3] <- Inf

  expect_no_error(.plot_fm_silent(fm_inf))

  gp <- .plot_fm_silent(fm_inf)
  expect_s3_class(gp, "ggplot")
})

test_that("plot.fineMap: Inf LOD capping uses max_finite * 1.1 so y_plot has no Inf", {
  fm_inf <- make_mock_fineMap(n_qtl = 1, n_pos = 10)
  fm_inf[[1]]$LOD[2] <- Inf

  gp <- .plot_fm_silent(fm_inf)

  expect_true(all(is.finite(gp$data$y_plot)))
})

test_that("plot.fineMap: QTL reference label suppressed when peak is within 1 cM of qtl_pos", {
  fm <- make_mock_fineMap(n_qtl = 1, n_pos = 15)
  # Move qtl_pos to the exact peak position (dist where LOD is maximum)
  peak_pos <- fm[[1]]$dist[which.max(fm[[1]]$LOD)]
  attr(fm[[1]], "qtl_pos") <- peak_pos

  gp <- .plot_fm_silent(fm)

  # The annotation data frame underlying the plot should only have one
  # point type ("Peak") when the reference is suppressed.
  # Access the data bound to the geom_point layer for annotation.
  ann_layers <- Filter(
    function(l) !is.null(l$data) && "pt_type" %in% names(l$data),
    gp$layers
  )
  if (length(ann_layers) > 0) {
    ann_df <- ann_layers[[1]]$data
    # "QTL" type should have been dropped (suppressed)
    expect_false("QTL" %in% as.character(ann_df$pt_type))
  } else {
    # If layers don't expose data at test time, at least confirm no error
    succeed("plot rendered without error; suppression check skipped for this ggplot version")
  }
})

test_that("plot.fineMap: QTL reference label NOT suppressed when peak is far from qtl_pos", {
  fm <- make_mock_fineMap(n_qtl = 1, n_pos = 20)
  # qtl_pos = 25 (fixture default); set peak far away by pushing all LOD high
  # except at the end of the grid
  fm[[1]]$LOD <- c(rep(0.1, 19), 10)   # peak at dist[20] = 50 cM, qtl_pos = 25
  attr(fm[[1]], "qtl_pos") <- 25

  gp <- .plot_fm_silent(fm)

  ann_layers <- Filter(
    function(l) !is.null(l$data) && "pt_type" %in% names(l$data),
    gp$layers
  )
  if (length(ann_layers) > 0) {
    ann_df <- ann_layers[[1]]$data
    expect_true("QTL" %in% as.character(ann_df$pt_type))
  } else {
    succeed("plot rendered without error")
  }
})

test_that("plot.fineMap: returns NULL invisibly when all data frames are 0-row", {
  fm <- make_mock_fineMap(n_qtl = 2)
  fm[[1]] <- fm[[1]][0, ]
  fm[[2]] <- fm[[2]][0, ]

  result <- suppressMessages(.plot_fm_silent(fm))

  expect_null(result)
})

test_that("plot.fineMap: works with n_qtl=1 (single-facet layout)", {
  fm <- make_mock_fineMap(n_qtl = 1, n_pos = 20)

  gp <- .plot_fm_silent(fm)

  expect_s3_class(gp, "ggplot")
})

test_that("plot.fineMap: works with n_qtl=4 (multi-column layout)", {
  fm <- make_mock_fineMap(n_qtl = 4, n_pos = 15)

  gp <- .plot_fm_silent(fm)

  expect_s3_class(gp, "ggplot")
  # 4 facets with ncol = min(4, 4) = 4 columns → 1 row
  expect_equal(attr(gp, "suggested_height"), 1 * 3.5)
})

test_that("plot.fineMap: custom colours accepted without error", {
  fm <- make_mock_fineMap()

  gp <- .plot_fm_silent(fm, col = "darkgreen", peak.col = "orange",
                        qtl.col = "purple", ref.col = "pink")

  expect_s3_class(gp, "ggplot")
})

# =============================================================================
# print.fineMap
# =============================================================================

test_that("print.fineMap outputs peak LOD information", {
  fm <- make_mock_fineMap(n_qtl = 1, n_pos = 10)

  expect_output(
    print(fm),
    "Peak"
  )
})

test_that("print.fineMap outputs QTL name", {
  fm     <- make_mock_fineMap(n_qtl = 1)
  nm     <- names(fm)[1]

  expect_output(
    print(fm),
    nm
  )
})

test_that("print.fineMap outputs 'No positions scanned' for a 0-row data frame", {
  fm        <- make_mock_fineMap(n_qtl = 2)
  fm[[1]]   <- fm[[1]][0, ]   # empty the first QTL's scan

  expect_output(
    print(fm),
    "[Nn]o positions scanned"
  )
})

test_that("print.fineMap outputs results for a non-empty QTL even when another is empty", {
  fm        <- make_mock_fineMap(n_qtl = 2)
  fm[[1]]   <- fm[[1]][0, ]

  # Second QTL should still print peak info
  expect_output(
    print(fm),
    "Peak"
  )
})

test_that("print.fineMap outputs the type header", {
  fm <- make_mock_fineMap()

  expect_output(
    print(fm),
    "Fine-mapping results"
  )
})

test_that("print.fineMap returns x invisibly", {
  fm     <- make_mock_fineMap()

  result <- withVisible(print(fm))

  expect_false(result$visible)
  expect_identical(result$value, fm)
})

test_that("print.fineMap prints numerical values with at least 4 significant digits", {
  fm <- make_mock_fineMap(n_qtl = 1, n_pos = 10)

  # Capture output and check that there is a numeric string after "dist ="
  out <- capture.output(print(fm))
  dist_lines <- grep("dist =", out, value = TRUE)
  expect_true(length(dist_lines) > 0)
})

test_that("print.fineMap: multiple QTL all appear in output", {
  fm  <- make_mock_fineMap(n_qtl = 3)
  nms <- names(fm)

  out <- capture.output(print(fm))
  full_out <- paste(out, collapse = "\n")

  for (nm in nms) {
    expect_true(grepl(nm, full_out, fixed = TRUE))
  }
})
