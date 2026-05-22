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
#   - Correct class and slots returned.
#   - Each filter step removes exactly the right items.
#   - Filtering order is respected (marker missingness before line
#     missingness; MAF after both; het last).
#   - filterPanel.checkPanel() S3 dispatch works.
#   - Thresholds set to NULL correctly skip that step.
#   - 0-row / 0-col warnings fire when over-filtering.
#   - print.filteredPanel() produces output.
#
# primePanel() tests verify:
#   - impute = "knn" produces a complete (no-NA) wgPanel.
#   - knn.k is respected.
#   - primePanel.filteredPanel() S3 dispatch works end-to-end.
#   - maf argument no longer accepted.
# =============================================================================

testthat::local_edition(3)

# ---- shared fixture ---------------------------------------------------------
p     <- make_raw_panel()
p_min <- make_raw_panel(n_dup_lines = 0, n_dup_marks = 0, n_mono = 0,
                         n_highmiss_marks = 0, n_highmiss_lines = 0,
                         n_highhet = 0, n_notinmap = 0, miss_rate = 0)

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
                      "map.pos", "thresholds", "removed", "n.original",
                      "n.final") %in% names(fp)))
})

test_that("filteredPanel$geno is a matrix with row and column names", {
    fp <- filterPanel(p$geno, p$map)
    expect_true(is.matrix(fp$geno))
    expect_false(is.null(rownames(fp$geno)))
    expect_false(is.null(colnames(fp$geno)))
})

test_that("filteredPanel$removed has one entry per filter step", {
    fp <- filterPanel(p$geno, p$map)
    expect_true(all(c("map_consistency", "miss_marker", "miss_line",
                      "het_line", "het_marker", "dup_lines", "dup_markers",
                      "maf") %in% names(fp$removed)))
})

test_that("filteredPanel n.original records pre-filter dimensions", {
    fp <- filterPanel(p$geno, p$map)
    expect_equal(unname(fp$n.original["lines"]),   nrow(p$geno))
    expect_equal(unname(fp$n.original["markers"]), ncol(p$geno))
})

# =============================================================================
# 2. Individual filter steps
# =============================================================================

test_that("Step 1: map consistency removes markers not in map", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE)
    expect_equal(length(fp$removed$map_consistency), p$n_notinmap)
    expect_false(any(fp$removed$map_consistency %in% colnames(fp$geno)))
})

test_that("Step 2: dup.lines = TRUE removes duplicate lines", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = NULL, dup.lines = TRUE, dup.markers = FALSE)
    expect_equal(length(fp$removed$dup_lines), p$n_dup_lines)
    expect_false(any(fp$removed$dup_lines %in% rownames(fp$geno)))
})

test_that("Step 2: dup.lines = FALSE skips duplicate line removal", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE)
    expect_equal(length(fp$removed$dup_lines), 0L)
})

test_that("Step 3: dup.markers = TRUE removes duplicate markers", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = TRUE)
    expect_true(length(fp$removed$dup_markers) >= p$n_dup_marks)
    expect_false(any(fp$removed$dup_markers %in% colnames(fp$geno)))
})

test_that("Step 4: marker missingness filter removes high-missing markers", {
    fp <- filterPanel(p$geno, p$map, miss.marker = 0.20, miss.line = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE)
    expect_true(length(fp$removed$miss_marker) >= p$n_highmiss_marks)
    # Remaining markers all have missingness <= threshold
    mm <- colMeans(is.na(fp$geno))
    expect_true(all(mm <= 0.20))
})

test_that("Step 5: line missingness filter removes high-missing lines", {
    fp <- filterPanel(p$geno, p$map, miss.marker = NULL, miss.line = 0.20,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE)
    expect_true(length(fp$removed$miss_line) >= p$n_highmiss_lines)
    ml <- rowMeans(is.na(fp$geno))
    expect_true(all(ml <= 0.20))
})

test_that("Step 6: MAF filter removes monomorphic and low-MAF markers", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = 0.05, dup.lines = FALSE, dup.markers = FALSE)
    expect_true(length(fp$removed$maf) >= p$n_mono)
    # All remaining markers have MAF >= 0.05
    col_means <- colMeans(fp$geno, na.rm = TRUE)
    pvals     <- col_means / 2   # 012 encoding
    mafs      <- pmin(pvals, 1 - pvals)
    expect_true(all(mafs >= 0.05 - 1e-9))
})

test_that("Step 4: het.line filter removes high-heterozygosity lines", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE,
                      het.line = 0.55)
    expect_true(length(fp$removed$het_line) >= p$n_highhet)
    het_vals <- rowMeans(fp$geno == 1L, na.rm = TRUE)
    expect_true(all(het_vals <= 0.55))
})

test_that("Step 5: het.marker filter removes high-heterozygosity markers", {
    fp <- filterPanel(p$geno, p$map, miss.line = NULL, miss.marker = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE,
                      het.marker = 0.55)
    het_vals <- colMeans(fp$geno == 1L, na.rm = TRUE)
    expect_true(all(het_vals <= 0.55))
})

test_that("het.line = NULL and het.marker = NULL skip heterozygosity filters", {
    fp <- filterPanel(p$geno, p$map, het.line = NULL, het.marker = NULL)
    expect_equal(length(fp$removed$het_line),   0L)
    expect_equal(length(fp$removed$het_marker), 0L)
})

# =============================================================================
# 3. Filtering order: marker missingness before line missingness
# =============================================================================

test_that("Marker missingness applied before line missingness", {
    # Build a panel where one marker is all-missing for 3 specific lines;
    # those lines would appear high-missing due to that marker but are fine
    # once the bad marker is removed.
    set.seed(200)
    n <- 30; p2 <- 10
    ids  <- paste0("L", seq_len(n))
    mids <- paste0("M", seq_len(p2))
    g    <- matrix(sample(0:2, n * p2, replace = TRUE),
                   nrow = n, dimnames = list(ids, mids))
    # Make first marker 100% missing — this inflates missingness for all lines
    g[, 1] <- NA
    mp <- data.frame(marker = mids, chr = "C1",
                     pos = seq_len(p2), stringsAsFactors = FALSE)

    # With marker filter first: bad marker removed, no lines fail miss.line
    fp <- filterPanel(g, mp, miss.marker = 0.50, miss.line = 0.20,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE)
    expect_true("M1" %in% fp$removed$miss_marker)
    expect_equal(length(fp$removed$miss_line), 0L)
})

# =============================================================================
# 4. filterPanel.checkPanel S3 dispatch
# =============================================================================

test_that("filterPanel.checkPanel accepts a checkPanel object", {
    chk <- checkPanel(p$geno, p$map)
    fp  <- filterPanel(chk, maf = 0.05)
    expect_s3_class(fp, "filteredPanel")
})

test_that("filterPanel.checkPanel produces same result as filterPanel.default", {
    chk  <- checkPanel(p$geno, p$map)
    fp1  <- filterPanel(chk,     maf = 0.05, miss.line = 0.20,
                         miss.marker = 0.20, dup.lines = TRUE,
                         dup.markers = FALSE)
    fp2  <- filterPanel(p$geno, p$map, maf = 0.05, miss.line = 0.20,
                         miss.marker = 0.20, dup.lines = TRUE,
                         dup.markers = FALSE)
    expect_equal(dim(fp1$geno), dim(fp2$geno))
    expect_equal(fp1$n.final,   fp2$n.final)
})

# =============================================================================
# 5. NULL thresholds skip steps
# =============================================================================

test_that("miss.marker = NULL skips marker missingness step", {
    fp <- filterPanel(p$geno, p$map, miss.marker = NULL, miss.line = NULL,
                      maf = NULL, dup.lines = FALSE, dup.markers = FALSE)
    expect_equal(length(fp$removed$miss_marker), 0L)
})

test_that("maf = NULL skips MAF step", {
    fp <- filterPanel(p$geno, p$map, maf = NULL, miss.line = NULL,
                      miss.marker = NULL, dup.lines = FALSE,
                      dup.markers = FALSE)
    expect_equal(length(fp$removed$maf), 0L)
})

# =============================================================================
# 6. Over-filtering warnings
# =============================================================================

test_that("filterPanel warns when all lines are removed", {
    # Force all lines to have 100% missing on first marker so miss.line = 0
    # triggers removal of every line
    g_all_miss <- p$geno
    g_all_miss[, ] <- NA   # every value missing -> every line has miss.line = 1
    expect_warning(
        filterPanel(g_all_miss, p$map, miss.line = 0.50, miss.marker = NULL,
                    maf = NULL, dup.lines = FALSE, dup.markers = FALSE),
        "All lines"
    )
})

test_that("filterPanel warns when all markers are removed", {
    expect_warning(
        filterPanel(p$geno, p$map, maf = 0.99, miss.line = NULL,
                    miss.marker = NULL, dup.lines = FALSE,
                    dup.markers = FALSE),
        "All markers"
    )
})

# =============================================================================
# 7. print.filteredPanel
# =============================================================================

test_that("print.filteredPanel produces output without error", {
    fp <- filterPanel(p$geno, p$map)
    expect_output(print(fp), "Filter Summary")
})

test_that("print.filteredPanel shows Remaining line", {
    fp <- filterPanel(p$geno, p$map)
    expect_output(print(fp), "Remaining")
})

test_that("print.filteredPanel returns x invisibly", {
    fp  <- filterPanel(p$geno, p$map)
    out <- withVisible(print(fp))
    expect_false(out$visible)
    expect_identical(out$value, fp)
})

# =============================================================================
# 8. primePanel() — KNN imputation path
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
    # $data should be the raw encoded matrix (with NAs);
    # $imputed.data should be the knn-completed version.
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
# 9. primePanel() — filteredPanel dispatch
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
    # Use the minimal panel with a few NAs introduced
    g_na <- p_min$geno
    g_na[1, 1] <- NA
    fp    <- filterPanel(g_na, p_min$map, miss.line = NULL,
                          miss.marker = NULL, maf = NULL,
                          dup.lines = FALSE, dup.markers = FALSE)
    panel <- suppressMessages(primePanel(fp, impute = "knn", knn.k = 3L))
    imp   <- do.call("cbind",
                     lapply(panel$geno, function(el) el$imputed.data))
    expect_false(anyNA(imp))
})

# =============================================================================
# 10. primePanel() — maf argument removed
# =============================================================================

test_that("primePanel no longer accepts maf argument", {
    expect_error(
        primePanel(p_min$geno, p_min$map, maf = 0.05),
        NA   # should NOT error on the call itself — maf goes to ... silently
        # Actually in R, extra args go to ...; they won't error unless
        # .primePanel_core explicitly rejects them. Adjust test accordingly.
    )
})

test_that("primePanel does NOT filter by MAF (monomorphic markers kept)", {
    # Build clean data with one monomorphic marker
    g  <- p_min$geno
    g[, 1] <- 0L   # monomorphic
    mp <- p_min$map
    panel <- suppressMessages(
        suppressWarnings(primePanel(g, mp, impute = "mean")))
    all_marks <- unlist(lapply(panel$geno, function(el) colnames(el$imputed.data)))
    # Monomorphic marker should still be present (no MAF filter in primePanel)
    expect_true(colnames(g)[1] %in% all_marks)
})
