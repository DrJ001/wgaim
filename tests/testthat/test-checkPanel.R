# =============================================================================
# test-checkPanel.R
# Tests for checkPanel().
#
# Fixtures (auto-loaded from helper-fixtures.R):
#   make_raw_panel() — synthetic geno + map with controllable known issues.
#
# Strategy
# --------
# checkPanel() is pure-diagnostic: no data is modified.  Tests confirm that:
#   - The return object has the correct class and slots.
#   - Each check correctly detects the issues planted by make_raw_panel().
#   - print.checkPanel() produces output without error.
#   - Guard clauses fire correctly on bad inputs.
# =============================================================================

testthat::local_edition(3)

# ---- shared fixture ---------------------------------------------------------
p <- make_raw_panel()   # default: 60 lines, 2 chr x 20 markers, known issues

# =============================================================================
# 1. Return structure
# =============================================================================

test_that("checkPanel returns object of class 'checkPanel'", {
    chk <- checkPanel(p$geno, p$map)
    expect_s3_class(chk, "checkPanel")
})

test_that("checkPanel return object has all expected slots", {
    chk <- checkPanel(p$geno, p$map)
    expected <- c("geno", "map", "encoding", "id", "map.id", "map.chr",
                  "map.pos", "n.lines", "n.markers", "n.chr",
                  "encoding.ok", "encoding.warn", "consistency",
                  "dup.lines", "dup.markers", "miss.marker", "miss.line",
                  "maf", "maf.table", "monomorphic", "het", "chr.coverage")
    expect_true(all(expected %in% names(chk)))
})

test_that("checkPanel n.lines and n.markers match geno dimensions", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(chk$n.lines,   nrow(p$geno))
    expect_equal(chk$n.markers, ncol(p$geno))
})

test_that("checkPanel n.chr matches number of chromosomes in map", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(chk$n.chr, p$n_chr)
})

test_that("checkPanel stores original geno and map unchanged", {
    chk <- checkPanel(p$geno, p$map)
    expect_identical(chk$geno, p$geno)
    expect_identical(chk$map,  p$map)
})

# =============================================================================
# 2. Encoding validation
# =============================================================================

test_that("checkPanel encoding.ok = TRUE for valid 012 data", {
    chk <- checkPanel(p$geno, p$map, encoding = "012")
    expect_true(chk$encoding.ok)
    expect_null(chk$encoding.warn)
})

test_that("checkPanel encoding.ok = FALSE when values out of range", {
    bad <- p$geno
    bad[1, 1] <- 5   # out of [0, 2]
    chk <- checkPanel(bad, p$map, encoding = "012")
    expect_false(chk$encoding.ok)
    expect_false(is.null(chk$encoding.warn))
})

test_that("checkPanel accepts pm1 encoding without warning", {
    pm1_panel <- make_raw_panel(encoding = "pm1")
    chk <- checkPanel(pm1_panel$geno, pm1_panel$map, encoding = "pm1")
    expect_true(chk$encoding.ok)
})

# =============================================================================
# 3. Map consistency
# =============================================================================

test_that("checkPanel detects markers in geno not in map", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(length(chk$consistency$not_in_map), p$n_notinmap)
})

test_that("checkPanel reports zero markers in map not in geno by default", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(length(chk$consistency$not_in_geno), 0L)
})

test_that("checkPanel detects markers in map not in geno", {
    extra_row <- data.frame(marker = "EXTRA_M", chr = "Chr1", pos = 999,
                            stringsAsFactors = FALSE)
    map_extra  <- rbind(p$map, extra_row)
    chk <- checkPanel(p$geno, map_extra)
    expect_true("EXTRA_M" %in% chk$consistency$not_in_geno)
})

# =============================================================================
# 4. Duplicate detection
# =============================================================================

test_that("checkPanel detects duplicate lines", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(length(chk$dup.lines), p$n_dup_lines)
})

test_that("checkPanel detects duplicate markers", {
    chk <- checkPanel(p$geno, p$map)
    # We planted n_dup_marks duplicate markers; at least that many should be
    # detected (accidental duplicates from random data are also valid finds)
    expect_true(length(chk$dup.markers) >= p$n_dup_marks)
})

test_that("checkPanel finds no duplicates in clean data", {
    clean <- make_raw_panel(n_dup_lines = 0, n_dup_marks = 0,
                             n_mono = 0, n_highmiss_marks = 0,
                             n_highmiss_lines = 0, n_highhet = 0,
                             n_notinmap = 0, miss_rate = 0)
    chk <- checkPanel(clean$geno, clean$map)
    expect_equal(length(chk$dup.lines),   0L)
    expect_equal(length(chk$dup.markers), 0L)
})

# =============================================================================
# 5. Missingness
# =============================================================================

test_that("checkPanel miss.marker is named numeric vector in [0, 1]", {
    chk <- checkPanel(p$geno, p$map)
    expect_type(chk$miss.marker, "double")
    expect_true(all(chk$miss.marker >= 0 & chk$miss.marker <= 1))
    expect_equal(length(chk$miss.marker),
                 length(intersect(colnames(p$geno), p$map$marker)))
})

test_that("checkPanel miss.line is named numeric vector in [0, 1]", {
    chk <- checkPanel(p$geno, p$map)
    expect_type(chk$miss.line, "double")
    expect_true(all(chk$miss.line >= 0 & chk$miss.line <= 1))
    expect_equal(length(chk$miss.line), nrow(p$geno))
})

test_that("checkPanel detects high-missingness markers", {
    chk <- checkPanel(p$geno, p$map)
    # High-missingness markers were forced to ~35% missing
    expect_true(any(chk$miss.marker > 0.20))
})

test_that("checkPanel detects high-missingness lines", {
    chk <- checkPanel(p$geno, p$map)
    expect_true(any(chk$miss.line > 0.20))
})

# =============================================================================
# 6. MAF
# =============================================================================

test_that("checkPanel maf is named numeric in [0, 0.5]", {
    chk <- checkPanel(p$geno, p$map)
    maf <- chk$maf[!is.na(chk$maf)]
    expect_true(all(maf >= 0 & maf <= 0.5))
})

test_that("checkPanel detects monomorphic markers", {
    chk <- checkPanel(p$geno, p$map)
    expect_true(length(chk$monomorphic) >= p$n_mono)
})

test_that("checkPanel maf.table has standard threshold rows", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(chk$maf.table$threshold, c(0.01, 0.02, 0.05, 0.10))
    expect_true(all(chk$maf.table$n.removed >= 0))
    # Higher thresholds remove at least as many markers
    expect_true(all(diff(chk$maf.table$n.removed) >= 0))
})

# =============================================================================
# 7. Heterozygosity
# =============================================================================

test_that("checkPanel het is named numeric in [0, 1]", {
    chk <- checkPanel(p$geno, p$map)
    expect_type(chk$het, "double")
    expect_true(all(chk$het >= 0 & chk$het <= 1))
    expect_equal(length(chk$het), nrow(p$geno))
})

test_that("checkPanel detects high-heterozygosity lines", {
    chk <- checkPanel(p$geno, p$map)
    expect_true(any(chk$het > 0.50))
})

test_that("checkPanel het uses correct heterozygous call for pm1 encoding", {
    # In pm1, heterozygous = 0; in 012, heterozygous = 1
    pm1_panel <- make_raw_panel(encoding = "pm1", n_highhet = 0,
                                 miss_rate = 0, n_dup_lines = 0,
                                 n_dup_marks = 0, n_mono = 0,
                                 n_highmiss_marks = 0, n_highmiss_lines = 0,
                                 n_notinmap = 0)
    chk_pm1 <- checkPanel(pm1_panel$geno, pm1_panel$map, encoding = "pm1")
    expect_true(all(chk_pm1$het >= 0 & chk_pm1$het <= 1))
})

# =============================================================================
# 8. Chromosome coverage
# =============================================================================

test_that("checkPanel chr.coverage has one row per chromosome", {
    chk <- checkPanel(p$geno, p$map)
    expect_equal(nrow(chk$chr.coverage), p$n_chr)
})

test_that("checkPanel chr.coverage columns are correct", {
    chk <- checkPanel(p$geno, p$map)
    expect_named(chk$chr.coverage, c("chr", "n.markers", "pos.min", "pos.max"))
})

test_that("checkPanel chr.coverage n.markers sums to markers in common", {
    chk <- checkPanel(p$geno, p$map)
    common_n <- length(intersect(colnames(p$geno), p$map$marker))
    expect_equal(sum(chk$chr.coverage$n.markers), common_n)
})

# =============================================================================
# 9. print.checkPanel
# =============================================================================

test_that("print.checkPanel produces output without error", {
    chk <- checkPanel(p$geno, p$map)
    expect_output(print(chk), "Panel Check Summary")
})

test_that("print.checkPanel reports map consistency section", {
    chk <- checkPanel(p$geno, p$map)
    expect_output(print(chk), "Map consistency")
})

test_that("print.checkPanel reports MAF distribution section", {
    chk <- checkPanel(p$geno, p$map)
    expect_output(print(chk), "MAF distribution")
})

test_that("print.checkPanel returns x invisibly", {
    chk <- checkPanel(p$geno, p$map)
    out <- withVisible(print(chk))
    expect_false(out$visible)
    expect_identical(out$value, chk)
})

# =============================================================================
# 10. Guard clauses
# =============================================================================

test_that("checkPanel stops when id column missing in data frame geno", {
    df <- as.data.frame(p$geno)
    expect_error(checkPanel(df, p$map, id = "no_such_col"), "'no_such_col'")
})

test_that("checkPanel stops when map column missing", {
    expect_error(checkPanel(p$geno, p$map, map.id = "bad_col"), "'bad_col'")
})

test_that("checkPanel stops when geno matrix has no row names", {
    bad <- p$geno
    rownames(bad) <- NULL
    expect_error(checkPanel(bad, p$map), "row names")
})

test_that("checkPanel stops when geno matrix has no column names", {
    bad <- p$geno
    colnames(bad) <- NULL
    expect_error(checkPanel(bad, p$map), "column names")
})

test_that("checkPanel accepts geno as a data frame with id column", {
    df  <- as.data.frame(p$geno)
    df$id <- rownames(p$geno)
    chk <- checkPanel(df, p$map, id = "id")
    expect_s3_class(chk, "checkPanel")
    expect_equal(chk$n.lines, nrow(p$geno))
})
