# =============================================================================
# test-filterPanel.R
# Tests for filterPanel() and the primePanel() KNN / filteredPanel paths.
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_raw_panel() — synthetic geno + map with known issues.
#
# Strategy
# --------
# filterPanel() tests verify:
#   - Correct class and slots returned ($geno, $map, $history, $n.original,
#     $n.final).
#   - Each filter step removes exactly the right items.
#   - Custom step list order is respected.
#   - filterPanel.checkPanel() S3 dispatch works.
#   - filterPanel.filteredPanel() S3 dispatch appends a new pass to $history.
#   - Steps with NULL / FALSE values are skipped.
#   - 0-row / 0-col warnings fire when over-filtering.
#   - print.filteredPanel() produces correct output for single and multi-pass
#     objects.
#
# primePanel() tests verify:
#   - impute = "knn" produces a complete (no-NA) wgPanel.
#   - knn.k is respected.
#   - primePanel.filteredPanel() S3 dispatch works end-to-end.
#   - maf argument no longer accepted.
# =============================================================================

testthat::local_edition(3)

# ---- shared fixtures --------------------------------------------------------
p     <- make_raw_panel()
p_min <- make_raw_panel(n_dup_lines = 0, n_dup_marks = 0, n_mono = 0,
                         n_highmiss_marks = 0, n_highmiss_lines = 0,
                         n_highhet = 0, n_notinmap = 0, miss_rate = 0)

# Convenience: default steps list (mirrors formals(filterPanel.default)$steps)
default_steps <- list(
    map         = TRUE,
    miss.marker = 0.20,
    miss.line   = 0.20,
    het.line    = NULL,
    het.marker  = NULL,
    dup.lines   = TRUE,
    dup.markers = FALSE,
    maf         = 0.05
)

# =============================================================================
# 1. Return structure
# =============================================================================

test_that("filterPanel returns object of class 'filteredPanel'", {
    fp <- filterPanel(p$geno, p$map)
    expect_s3_class(fp, "filteredPanel")
})

test_that("filteredPanel has all expected slots", {
    fp <- filterPanel(p$geno, p$map)
    expect_true(all(c("geno", "map", "encoding", "id", "map.id", "map.chr",
                      "map.pos", "history", "n.original",
                      "n.final") %in% names(fp)))
})

test_that("filteredPanel does NOT have old $removed or $thresholds slots", {
    fp <- filterPanel(p$geno, p$map)
    expect_false("removed"    %in% names(fp))
    expect_false("thresholds" %in% names(fp))
})

test_that("filteredPanel$geno is a matrix with row and column names", {
    fp <- filterPanel(p$geno, p$map)
    expect_true(is.matrix(fp$geno))
    expect_false(is.null(rownames(fp$geno)))
    expect_false(is.null(colnames(fp$geno)))
})

test_that("filteredPanel$history is a list with one element after single call", {
    fp <- filterPanel(p$geno, p$map)
    expect_type(fp$history, "list")
    expect_equal(length(fp$history), 1L)
})

test_that("filteredPanel history pass 1 has expected fields", {
    fp   <- filterPanel(p$geno, p$map)
    pass <- fp$history[[1L]]
    expect_true(all(c("pass", "steps", "removed", "n.before",
                      "n.after") %in% names(pass)))
    expect_equal(pass$pass, 1L)
})

test_that("filteredPanel n.original records pre-filter dimensions", {
    fp <- filterPanel(p$geno, p$map)
    expect_equal(unname(fp$n.original["lines"]),   nrow(p$geno))
    expect_equal(unname(fp$n.original["markers"]), ncol(p$geno))
})

test_that("filteredPanel n.before in pass 1 equals n.original", {
    fp <- filterPanel(p$geno, p$map)
    expect_equal(fp$history[[1L]]$n.before, fp$n.original)
})

# =============================================================================
# 2. Individual filter steps (via steps = list(...))
# =============================================================================

test_that("map consistency always runs first and removes markers not in map", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE))
    expect_equal(length(fp$history[[1L]]$removed$map_consistency),
                 p$n_notinmap)
    expect_false(any(fp$history[[1L]]$removed$map_consistency %in%
                         colnames(fp$geno)))
})

test_that("miss.marker step removes high-missing markers", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, miss.marker = 0.20))
    rem <- fp$history[[1L]]$removed$miss.marker
    expect_true(length(rem) >= p$n_highmiss_marks)
    mm <- colMeans(is.na(fp$geno))
    expect_true(all(mm <= 0.20))
})

test_that("miss.line step removes high-missing lines", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, miss.line = 0.20))
    rem <- fp$history[[1L]]$removed$miss.line
    expect_true(length(rem) >= p$n_highmiss_lines)
    ml <- rowMeans(is.na(fp$geno))
    expect_true(all(ml <= 0.20))
})

test_that("dup.lines = TRUE removes duplicate lines", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, dup.lines = TRUE))
    rem <- fp$history[[1L]]$removed$dup.lines
    expect_equal(length(rem), p$n_dup_lines)
    expect_false(any(rem %in% rownames(fp$geno)))
})

test_that("dup.lines = FALSE skips duplicate line removal", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, dup.lines = FALSE))
    expect_equal(length(fp$history[[1L]]$removed$dup.lines), 0L)
})

test_that("dup.markers = TRUE removes duplicate markers", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, dup.markers = TRUE))
    rem <- fp$history[[1L]]$removed$dup.markers
    expect_true(length(rem) >= p$n_dup_marks)
    expect_false(any(rem %in% colnames(fp$geno)))
})

test_that("maf step removes monomorphic and low-MAF markers", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, maf = 0.05))
    rem <- fp$history[[1L]]$removed$maf
    expect_true(length(rem) >= p$n_mono)
    col_means <- colMeans(fp$geno, na.rm = TRUE)
    pvals     <- col_means / 2
    mafs      <- pmin(pvals, 1 - pvals)
    expect_true(all(mafs >= 0.05 - 1e-9))
})

test_that("het.line step removes high-heterozygosity lines", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, het.line = 0.55))
    rem <- fp$history[[1L]]$removed$het.line
    expect_true(length(rem) >= p$n_highhet)
    het_vals <- rowMeans(fp$geno == 1L, na.rm = TRUE)
    expect_true(all(het_vals <= 0.55))
})

test_that("het.marker step removes high-heterozygosity markers", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, het.marker = 0.55))
    het_vals <- colMeans(fp$geno == 1L, na.rm = TRUE)
    expect_true(all(het_vals <= 0.55))
})

test_that("het.line = NULL and het.marker = NULL are skipped", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, het.line = NULL,
                                   het.marker = NULL))
    expect_equal(length(fp$history[[1L]]$removed$het.line),   0L)
    expect_equal(length(fp$history[[1L]]$removed$het.marker), 0L)
})

# =============================================================================
# 3. Marker missingness before line missingness (order matters)
# =============================================================================

test_that("marker missingness applied before line missingness by default", {
    set.seed(200)
    n <- 30; p2 <- 10
    ids  <- paste0("L", seq_len(n))
    mids <- paste0("M", seq_len(p2))
    g    <- matrix(sample(0:2, n * p2, replace = TRUE),
                   nrow = n, dimnames = list(ids, mids))
    g[, 1] <- NA   # 100% missing marker inflates all line miss rates
    mp <- data.frame(marker = mids, chr = "C1",
                     pos = seq_len(p2), stringsAsFactors = FALSE)

    # Default order: miss.marker before miss.line
    fp <- filterPanel(g, mp,
                      steps = list(map         = TRUE,
                                   miss.marker = 0.50,
                                   miss.line   = 0.20))
    expect_true("M1" %in% fp$history[[1L]]$removed$miss.marker)
    expect_equal(length(fp$history[[1L]]$removed$miss.line), 0L)
})

test_that("custom step order is respected: miss.line before miss.marker", {
    set.seed(201)
    n <- 30; p2 <- 10
    ids  <- paste0("L", seq_len(n))
    mids <- paste0("M", seq_len(p2))
    g    <- matrix(sample(0:2, n * p2, replace = TRUE),
                   nrow = n, dimnames = list(ids, mids))
    g[, 1] <- NA
    mp <- data.frame(marker = mids, chr = "C1",
                     pos = seq_len(p2), stringsAsFactors = FALSE)

    # Reversed order: miss.line runs first (while bad marker is still there)
    fp <- filterPanel(g, mp,
                      steps = list(map       = TRUE,
                                   miss.line   = 0.20,
                                   miss.marker = 0.50))
    step_names <- names(fp$history[[1L]]$steps)
    miss_line_pos   <- which(step_names == "miss.line")
    miss_marker_pos <- which(step_names == "miss.marker")
    expect_true(miss_line_pos < miss_marker_pos)
})

# =============================================================================
# 4. filterPanel.checkPanel S3 dispatch
# =============================================================================

test_that("filterPanel.checkPanel accepts a checkPanel object", {
    chk <- checkPanel(p$geno, p$map)
    fp  <- filterPanel(chk, steps = list(map = TRUE, maf = 0.05))
    expect_s3_class(fp, "filteredPanel")
})

test_that("filterPanel.checkPanel produces same result as filterPanel.default", {
    chk  <- checkPanel(p$geno, p$map)
    fp1  <- filterPanel(chk,     steps = default_steps)
    fp2  <- filterPanel(p$geno, p$map, steps = default_steps)
    expect_equal(dim(fp1$geno), dim(fp2$geno))
    expect_equal(fp1$n.final,   fp2$n.final)
})

# =============================================================================
# 5. filterPanel.filteredPanel S3 dispatch — multi-pass / history
# =============================================================================

test_that("filterPanel on a filteredPanel returns a filteredPanel", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_s3_class(fp2, "filteredPanel")
})

test_that("second pass appends to $history (length 2)", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_equal(length(fp2$history), 2L)
    expect_equal(fp2$history[[2L]]$pass, 2L)
})

test_that("n.original is preserved across multiple passes", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_equal(fp2$n.original,
                 c(lines = nrow(p$geno), markers = ncol(p$geno)))
})

test_that("n.before in pass 2 equals n.after from pass 1", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_equal(fp2$history[[2L]]$n.before, fp1$n.final)
})

test_that("n.final reflects state after both passes", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_equal(fp2$n.final,
                 c(lines   = nrow(fp2$geno),
                   markers = ncol(fp2$geno)))
})

test_that("three-pass filtering accumulates three history entries", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1,   steps = list(miss.line = 0.20))
    fp3 <- filterPanel(fp2,   steps = list(maf = 0.05))
    expect_equal(length(fp3$history), 3L)
    expect_equal(fp3$history[[3L]]$pass, 3L)
})

# =============================================================================
# 6. Unknown step name causes an error
# =============================================================================

test_that("unknown step name in steps list causes an error", {
    expect_error(
        filterPanel(p$geno, p$map, steps = list(map = TRUE, foo = 0.1)),
        regexp = "Unknown step name"
    )
})

# =============================================================================
# 7. Map consistency always runs first even if omitted from custom steps list
# =============================================================================

test_that("map consistency is prepended even when not first in steps list", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(miss.marker = 0.20))
    step_names <- names(fp$history[[1L]]$steps)
    expect_equal(step_names[1L], "map")
})

# =============================================================================
# 8. NULL / FALSE thresholds skip steps
# =============================================================================

test_that("miss.marker = NULL skips that step", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, miss.marker = NULL))
    expect_equal(length(fp$history[[1L]]$removed$miss.marker), 0L)
})

test_that("maf = NULL skips MAF step", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, maf = NULL))
    expect_equal(length(fp$history[[1L]]$removed$maf), 0L)
})

test_that("dup.markers = FALSE skips duplicate marker removal", {
    fp <- filterPanel(p$geno, p$map,
                      steps = list(map = TRUE, dup.markers = FALSE))
    expect_equal(length(fp$history[[1L]]$removed$dup.markers), 0L)
})

# =============================================================================
# 9. Over-filtering warnings
# =============================================================================

test_that("filterPanel warns when all lines are removed", {
    g_all_miss        <- p$geno
    g_all_miss[, ]    <- NA
    expect_warning(
        filterPanel(g_all_miss, p$map,
                    steps = list(map = TRUE, miss.line = 0.50)),
        "All lines"
    )
})

test_that("filterPanel warns when all markers are removed", {
    expect_warning(
        filterPanel(p$geno, p$map,
                    steps = list(map = TRUE, maf = 0.99)),
        "All markers"
    )
})

# =============================================================================
# 10. print.filteredPanel
# =============================================================================

test_that("print.filteredPanel produces output without error", {
    fp <- filterPanel(p$geno, p$map)
    expect_output(print(fp), "Filter History")
})

test_that("print.filteredPanel shows After line for single pass", {
    fp <- filterPanel(p$geno, p$map)
    expect_output(print(fp), "After")
})

test_that("print.filteredPanel shows Remaining line for single pass", {
    fp <- filterPanel(p$geno, p$map)
    expect_output(print(fp), "Remaining")
})

test_that("print.filteredPanel returns x invisibly", {
    fp  <- filterPanel(p$geno, p$map)
    out <- withVisible(print(fp))
    expect_false(out$visible)
    expect_identical(out$value, fp)
})

test_that("print.filteredPanel shows pass headers for multi-pass object", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_output(print(fp2), "Pass 1")
    expect_output(print(fp2), "Pass 2")
})

test_that("print.filteredPanel shows Final line for multi-pass object", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_output(print(fp2), "Final")
})

test_that("print.filteredPanel shows [2 passes] in header for multi-pass", {
    fp1 <- filterPanel(p$geno, p$map,
                       steps = list(map = TRUE, miss.marker = 0.20))
    fp2 <- filterPanel(fp1, steps = list(maf = 0.05))
    expect_output(print(fp2), "2 passes")
})

# =============================================================================
# 11. primePanel() — KNN imputation path
# =============================================================================

test_that("primePanel impute='knn' returns wgPanel with no NAs", {
    set.seed(55)
    n <- 20; m <- 8
    ids  <- paste0("L", seq_len(n))
    mids <- paste0("M", seq_len(m))
    g    <- matrix(sample(0:2, n * m, replace = TRUE),
                   nrow = n, dimnames = list(ids, mids))
    g[sample(length(g), round(0.05 * length(g)))] <- NA
    mp  <- data.frame(marker = mids, chr = "C1",
                      pos = seq_len(m), stringsAsFactors = FALSE)
    panel <- suppressMessages(
        suppressWarnings(primePanel(g, mp, encoding = "012", impute = "knn"))
    )
    imp <- panel$geno[["C1"]]$imputed.data
    expect_false(anyNA(imp))
})

test_that("primePanel knn.k argument is respected (no error for k=1 or k=10)", {
    set.seed(56)
    n <- 15; m <- 6
    g <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                dimnames = list(paste0("L", seq_len(n)),
                                paste0("M", seq_len(m))))
    g[1, 1] <- NA
    mp <- data.frame(marker = paste0("M", seq_len(m)), chr = "C1",
                     pos = seq_len(m), stringsAsFactors = FALSE)
    expect_no_error(suppressMessages(
        primePanel(g, mp, encoding = "012", impute = "knn", knn.k = 1L)))
    expect_no_error(suppressMessages(
        primePanel(g, mp, encoding = "012", impute = "knn", knn.k = 10L)))
})

test_that("primePanel impute='knn' imputed values are in [-1, 1]", {
    set.seed(57)
    n <- 20; m <- 10
    g <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                dimnames = list(paste0("L", seq_len(n)),
                                paste0("M", seq_len(m))))
    g[sample(length(g), 5)] <- NA
    mp <- data.frame(marker = paste0("M", seq_len(m)), chr = "C1",
                     pos = seq_len(m), stringsAsFactors = FALSE)
    panel <- suppressMessages(
        primePanel(g, mp, encoding = "012", impute = "knn"))
    imp <- panel$geno[["C1"]]$imputed.data
    expect_true(all(imp >= -1 & imp <= 1, na.rm = TRUE))
})

test_that("primePanel impute='knn': data slot retains original NAs", {
    set.seed(58)
    n <- 15; m <- 6
    g <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n,
                dimnames = list(paste0("L", seq_len(n)),
                                paste0("M", seq_len(m))))
    g[1, 1] <- NA
    mp <- data.frame(marker = paste0("M", seq_len(m)), chr = "C1",
                     pos = seq_len(m), stringsAsFactors = FALSE)
    panel <- suppressMessages(
        primePanel(g, mp, encoding = "012", impute = "knn"))
    expect_true( anyNA(panel$geno[["C1"]]$data))
    expect_false(anyNA(panel$geno[["C1"]]$imputed.data))
})

# =============================================================================
# 12. primePanel() — filteredPanel dispatch
# =============================================================================

test_that("primePanel.filteredPanel returns wgPanel from a filteredPanel", {
    fp    <- filterPanel(p_min$geno, p_min$map)
    panel <- suppressMessages(primePanel(fp, impute = "none"))
    expect_s3_class(panel, "wgPanel")
})

test_that("primePanel.filteredPanel preserves line count", {
    fp    <- filterPanel(p_min$geno, p_min$map)
    panel <- suppressMessages(primePanel(fp, impute = "none"))
    expect_equal(nrow(panel$pheno), unname(fp$n.final["lines"]))
})

test_that("primePanel.filteredPanel preserves marker count", {
    fp    <- filterPanel(p_min$geno, p_min$map)
    panel <- suppressMessages(primePanel(fp, impute = "none"))
    n_markers <- sum(sapply(panel$geno, function(el) ncol(el$imputed.data)))
    expect_equal(n_markers, unname(fp$n.final["markers"]))
})

test_that("primePanel.filteredPanel works with impute='knn'", {
    g_na       <- p_min$geno
    g_na[1, 1] <- NA
    fp    <- filterPanel(g_na, p_min$map,
                         steps = list(map = TRUE))
    panel <- suppressMessages(primePanel(fp, impute = "knn", knn.k = 3L))
    imp   <- do.call("cbind",
                     lapply(panel$geno, function(el) el$imputed.data))
    expect_false(anyNA(imp))
})

# multi-pass filteredPanel flows through to primePanel
test_that("primePanel accepts a multi-pass filteredPanel", {
    fp1   <- filterPanel(p_min$geno, p_min$map,
                         steps = list(map = TRUE, miss.marker = 0.50))
    fp2   <- filterPanel(fp1, steps = list(maf = 0.01))
    panel <- suppressMessages(primePanel(fp2, impute = "none"))
    expect_s3_class(panel, "wgPanel")
})

# =============================================================================
# 13. primePanel() — maf argument removed
# =============================================================================

test_that("primePanel no longer accepts maf argument", {
    expect_error(
        primePanel(p_min$geno, p_min$map, maf = 0.05),
        "unused argument"
    )
})

test_that("primePanel does NOT filter by MAF (monomorphic markers kept)", {
    g       <- p_min$geno
    g[, 1]  <- 0L
    panel   <- suppressMessages(
        suppressWarnings(primePanel(g, p_min$map, impute = "mean")))
    all_marks <- unlist(lapply(panel$geno,
                               function(el) colnames(el$imputed.data)))
    expect_true(colnames(g)[1] %in% all_marks)
})
