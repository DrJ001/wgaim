# =============================================================================
# test-primeCross.R
# Tests for primeCross(), primePanel(), and the fixMap() helper.
#
# Fixtures from helper-fixtures.R (loaded automatically by testthat):
#   make_bc_cross(), make_dh_cross(), make_f2_cross()
#   make_wgCross_interval(), make_wgCross_marker()
# =============================================================================

# qtl is a Depends of wgAim and will always be present, but guard explicitly
# so the file is skippable in isolation if needed.
skip_if_not_installed("qtl")

# =============================================================================
# primeCross() — interval type
# =============================================================================

test_that("primeCross bc cross with type='interval' returns wgCross", {
  cr  <- make_bc_cross()
  wgc <- primeCross(cr, type = "interval", impute = "MartinezCurnow")

  expect_s3_class(wgc, "wgCross")
  expect_true(inherits(wgc, "cross"))
  expect_equal(wgc$type, "interval")
})

test_that("primeCross interval: every chromosome gains $imputed.data and $interval.data", {
  cr  <- make_bc_cross()
  wgc <- primeCross(cr, type = "interval", impute = "MartinezCurnow")

  for (chr in names(wgc$geno)) {
    expect_true(!is.null(wgc$geno[[chr]]$imputed.data),
                info = paste("imputed.data missing for chr", chr))
    expect_true(!is.null(wgc$geno[[chr]]$interval.data),
                info = paste("interval.data missing for chr", chr))
  }
})

# =============================================================================
# primeCross() — marker type
# =============================================================================

test_that("primeCross with type='marker' sets $type and omits $interval.data", {
  cr  <- make_bc_cross()
  wgc <- primeCross(cr, type = "marker", impute = "MartinezCurnow")

  expect_equal(wgc$type, "marker")
  for (chr in names(wgc$geno)) {
    expect_true(!is.null(wgc$geno[[chr]]$imputed.data),
                info = paste("imputed.data missing for chr", chr))
    expect_null(wgc$geno[[chr]]$interval.data,
                label = paste("interval.data should be absent in marker mode, chr", chr))
  }
})

# =============================================================================
# primeCross() — input validation
# =============================================================================

test_that("primeCross stops when object has unsupported class", {
  bad <- list(geno = list(), pheno = data.frame(id = "L1"))
  class(bad) <- "notacross"

  expect_error(
    primeCross(bad),
    regexp = "bc.*dh.*f2.*riself",
    ignore.case = TRUE
  )
})

test_that("primeCross stops when id column is missing from pheno", {
  cr <- make_bc_cross()
  cr$pheno$id <- NULL          # remove the 'id' column

  expect_error(
    primeCross(cr, id = "id"),
    regexp = "cannot be found"
  )
})

# =============================================================================
# primeCross() — subset argument
# =============================================================================

test_that("primeCross subset retains only the requested chromosomes", {
  cr      <- make_bc_cross(n_chr = 3)
  keep    <- names(cr$geno)[1:2]           # keep first two chromosomes
  wgc     <- primeCross(cr, type = "marker", subset = keep)

  expect_equal(sort(names(wgc$geno)), sort(keep))
})

# =============================================================================
# primeCross() — consensus.mark argument
# =============================================================================

test_that("primeCross consensus.mark=FALSE skips fixMap (no $colocated.markers)", {
  cr  <- make_bc_cross()
  wgc <- primeCross(cr, type = "marker", consensus.mark = FALSE)

  # fixMap adds $colocated.markers to the cross object; skipping it means
  # the slot will NOT exist (or will be NULL if never set).
  expect_null(wgc$colocated.markers)
})

# =============================================================================
# primeCross() — imputed.data values are in {-1, 0, 1} with no NA
# =============================================================================

test_that("primeCross imputed.data contains no NAs and values in {-1,0,1}", {
  # Use a cross with NO missing data to make the check unambiguous
  set.seed(1)
  cr <- make_bc_cross()
  # Remove all NAs from marker data
  for (i in seq_along(cr$geno)) {
    na_idx <- is.na(cr$geno[[i]]$data)
    cr$geno[[i]]$data[na_idx] <- sample(c(1L, 2L), sum(na_idx), replace = TRUE)
  }

  wgc <- primeCross(cr, type = "marker", impute = "MartinezCurnow")

  for (chr in names(wgc$geno)) {
    imp <- wgc$geno[[chr]]$imputed.data
    expect_false(anyNA(imp),
                 info = paste("NAs found in imputed.data for chr", chr))
    expect_true(all(imp %in% c(-1, 0, 1)),
                info = paste("Values outside {-1,0,1} for chr", chr))
  }
})

# =============================================================================
# primeCross() — DH cross
# =============================================================================

test_that("primeCross works with a dh cross", {
  cr  <- make_dh_cross()
  wgc <- primeCross(cr, type = "interval", impute = "MartinezCurnow")

  expect_s3_class(wgc, "wgCross")
  expect_true(inherits(wgc, "dh"))
  expect_equal(wgc$type, "interval")
})

# =============================================================================
# primeCross() — F2 cross (3-level genotype with dominance)
# =============================================================================

test_that("primeCross works with an f2 cross, returning wgCross", {
  cr  <- make_f2_cross()
  wgc <- primeCross(cr, type = "interval", impute = "MartinezCurnow")

  expect_s3_class(wgc, "wgCross")
  expect_true(inherits(wgc, "f2"))
  expect_equal(wgc$type, "interval")

  # f2 codes 1=AA (+1), 2=AB (0), 3=BB (-1); after imputation additive values
  # should be in [-1, 1] with no NA
  for (chr in names(wgc$geno)) {
    imp <- wgc$geno[[chr]]$imputed.data
    expect_false(anyNA(imp), info = paste("NAs in f2 imputed.data, chr", chr))
    expect_true(all(imp >= -1 & imp <= 1),
                info = paste("f2 imputed.data outside [-1,1], chr", chr))
  }
})

# =============================================================================
# fixMap() — reduces co-located markers; appends $colocated.markers
# =============================================================================

test_that("fixMap collapses co-located markers and adds $colocated.markers", {
  # Build a cross that has deliberate duplicated map positions on one chr
  cr <- make_bc_cross(n_lines = 20, n_chr = 1, n_mar = 4)

  # Force markers 2 and 3 on chr 1 to the same position
  cr$geno[[1]]$map[2] <- cr$geno[[1]]$map[3]

  fixed <- fixMap(cr, rd = 3)

  # After fixMap, the resulting cross should have fewer markers on that chr
  # and a $colocated.markers data.frame
  expect_true(!is.null(fixed$colocated.markers),
              label = "colocated.markers slot should exist after fixMap")
  expect_s3_class(fixed$colocated.markers, "data.frame")
  # The merged chr should have one fewer marker than the original
  expect_lt(length(fixed$geno[[1]]$map), 4L)
})

test_that("fixMap returns a cross with no extra markers when no co-location", {
  # A clean cross with all distinct positions
  cr  <- make_bc_cross(n_lines = 20, n_chr = 1, n_mar = 4)
  # Ensure positions are clearly distinct (they already are from make_bc_cross)
  original_nmar <- length(cr$geno[[1]]$map)

  fixed <- fixMap(cr, rd = 3)

  # No co-location: marker count unchanged on each chr
  expect_equal(length(fixed$geno[[1]]$map), original_nmar)
})

# =============================================================================
# primeCross() — class vector structure
# =============================================================================

test_that("primeCross prepends 'wgCross' to the existing class vector", {
  cr  <- make_bc_cross()
  wgc <- primeCross(cr, type = "marker")

  cls <- class(wgc)
  expect_equal(cls[1], "wgCross")
  expect_true("bc"    %in% cls)
  expect_true("cross" %in% cls)
})

# =============================================================================
# primePanel() — basic matrix input with encoding="012"
# =============================================================================

test_that("primePanel with matrix input returns wgPanel with correct class", {
  set.seed(21)
  n_lines <- 20; n_mar <- 10
  ids     <- paste0("L", seq_len(n_lines))
  gmat    <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                    nrow = n_lines,
                    dimnames = list(ids, paste0("M", seq_len(n_mar))))
  map_df  <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = rep(c("Chr1", "Chr2"), each = n_mar / 2),
    pos    = rep(seq(10, 100, length.out = n_mar / 2), times = 2),
    stringsAsFactors = FALSE
  )

  expect_message(
    panel <- primePanel(gmat, map_df, encoding = "012", maf = 0),
    regexp = "wgPanel object created"
  )

  expect_s3_class(panel, "wgPanel")
  expect_true(inherits(panel, "wgCross"))
  # primePanel always operates in marker mode; it does not set a $type field
  # (that is set only by primeCross).  Verify the class instead.
  expect_true(inherits(panel, "cross"))
})

# =============================================================================
# primePanel() — data.frame input with id column
# =============================================================================

test_that("primePanel accepts data.frame input with id column", {
  set.seed(22)
  n_lines <- 15; n_mar <- 8
  ids     <- paste0("S", seq_len(n_lines))
  gdf     <- as.data.frame(
    matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
           nrow = n_lines,
           dimnames = list(NULL, paste0("M", seq_len(n_mar))))
  )
  gdf$id  <- ids

  map_df  <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = rep("Chr1", n_mar),
    pos    = seq(5, 80, length.out = n_mar),
    stringsAsFactors = FALSE
  )

  expect_message(
    panel <- primePanel(gdf, map_df, id = "id", encoding = "012", maf = 0),
    regexp = "wgPanel object created"
  )

  expect_s3_class(panel, "wgPanel")
  expect_equal(panel$pheno$id, ids)
})

# =============================================================================
# primePanel() — encoding="pm1" passes through without shift
# =============================================================================

test_that("primePanel encoding='pm1' does not alter values", {
  set.seed(23)
  n_lines <- 10; n_mar <- 6
  ids     <- paste0("L", seq_len(n_lines))
  # Already in ±1 coding
  gmat    <- matrix(sample(c(-1L, 0L, 1L), n_lines * n_mar, replace = TRUE),
                    nrow = n_lines,
                    dimnames = list(ids, paste0("M", seq_len(n_mar))))
  map_df  <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "ChrA",
    pos    = seq(0, 50, length.out = n_mar),
    stringsAsFactors = FALSE
  )

  expect_message(
    panel <- primePanel(gmat, map_df, encoding = "pm1", maf = 0),
    regexp = "wgPanel object created"
  )

  # Values should not have been shifted — range must remain in [-1, 1]
  combined <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
  expect_true(all(combined >= -1 & combined <= 1))
})

# =============================================================================
# primePanel() — bad encoding stops
# =============================================================================

test_that("primePanel stops on unrecognised encoding", {
  set.seed(24)
  ids  <- paste0("L", 1:5)
  gmat <- matrix(sample(0:2, 5 * 4, replace = TRUE),
                 nrow = 5,
                 dimnames = list(ids, paste0("M", 1:4)))
  map_df <- data.frame(
    marker = paste0("M", 1:4), chr = "C1",
    pos = 1:4, stringsAsFactors = FALSE
  )

  expect_error(
    suppressMessages(primePanel(gmat, map_df, encoding = "binary")),
    regexp = "encoding must be"
  )
})

# =============================================================================
# primePanel() — MAF filter removes low-frequency markers
# =============================================================================

test_that("filterPanel MAF filter drops monomorphic markers", {
  # MAF filtering is now handled by filterPanel(), not primePanel().
  set.seed(25)
  n_lines <- 30; n_mar <- 10
  ids     <- paste0("L", seq_len(n_lines))
  gmat    <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                    nrow = n_lines,
                    dimnames = list(ids, paste0("M", seq_len(n_mar))))
  # Force the last marker to be monomorphic (MAF = 0)
  gmat[, n_mar] <- 0L

  map_df <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "Chr1",
    pos    = seq(1, n_mar),
    stringsAsFactors = FALSE
  )

  fp <- filterPanel(gmat, map_df, encoding = "012", maf = 0.05,
                    miss.line = NULL, miss.marker = NULL,
                    dup.lines = FALSE, dup.markers = FALSE)

  # Monomorphic marker must be removed
  expect_false(paste0("M", n_mar) %in% colnames(fp$geno))
  expect_true(length(fp$removed$maf) >= 1L)
})

# =============================================================================
# filterPanel() — maf=NULL skips MAF filter
# =============================================================================

test_that("filterPanel with maf=NULL keeps all markers", {
  set.seed(26)
  n_lines <- 20; n_mar <- 6
  ids  <- paste0("L", seq_len(n_lines))
  gmat <- matrix(0L, nrow = n_lines, ncol = n_mar,
                 dimnames = list(ids, paste0("M", seq_len(n_mar))))
  # All monomorphic — with maf=NULL they should survive
  map_df <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "ChrX",
    pos    = seq_len(n_mar),
    stringsAsFactors = FALSE
  )

  fp <- filterPanel(gmat, map_df, encoding = "012", maf = NULL,
                    miss.line = NULL, miss.marker = NULL,
                    dup.lines = FALSE, dup.markers = FALSE)

  expect_equal(ncol(fp$geno), n_mar)
  expect_equal(length(fp$removed$maf), 0L)
})

# =============================================================================
# primePanel() — missing data with impute="none" stops
# =============================================================================

test_that("primePanel stops with informative error when NAs remain and impute='none'", {
  set.seed(27)
  n_lines <- 20; n_mar <- 8
  ids  <- paste0("L", seq_len(n_lines))
  gmat <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                 nrow = n_lines,
                 dimnames = list(ids, paste0("M", seq_len(n_mar))))
  gmat[1, 1] <- NA   # introduce one missing genotype

  map_df <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "Chr1",
    pos    = seq_len(n_mar),
    stringsAsFactors = FALSE
  )

  expect_error(
    suppressMessages(primePanel(gmat, map_df, encoding = "012",
                                maf = 0, impute = "none")),
    regexp = "missing genotype"
  )
})

# =============================================================================
# primePanel() — missing data with impute="mean" emits warning and fills NAs
# =============================================================================

test_that("primePanel impute='mean' emits a warning and fills NAs", {
  set.seed(28)
  n_lines <- 20; n_mar <- 6
  ids  <- paste0("L", seq_len(n_lines))
  gmat <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                 nrow = n_lines,
                 dimnames = list(ids, paste0("M", seq_len(n_mar))))
  gmat[2, 3] <- NA   # one NA

  map_df <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "ChrB",
    pos    = seq_len(n_mar),
    stringsAsFactors = FALSE
  )

  expect_warning(
    suppressMessages(
      panel <- primePanel(gmat, map_df, encoding = "012",
                          maf = 0, impute = "mean")
    ),
    regexp = "Imputing"
  )

  # No NAs should remain after mean imputation
  all_data <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
  expect_false(anyNA(all_data))
})

# =============================================================================
# primePanel() — markers absent from map are dropped with a message
# =============================================================================

test_that("primePanel drops markers absent from map and emits a message", {
  set.seed(29)
  n_lines <- 20; n_mar <- 8
  ids  <- paste0("L", seq_len(n_lines))
  # Create gmat with n_mar+2 markers; map only covers the first n_mar
  gmat <- matrix(sample(0:2, n_lines * (n_mar + 2), replace = TRUE),
                 nrow = n_lines,
                 dimnames = list(ids, paste0("M", seq_len(n_mar + 2))))

  map_df <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "Chr1",
    pos    = seq_len(n_mar),
    stringsAsFactors = FALSE
  )

  expect_message(
    panel <- primePanel(gmat, map_df, encoding = "012", maf = 0),
    regexp = "not found in map"
  )

  all_markers <- unlist(lapply(panel$geno, function(el) names(el$map)))
  # Extra markers M9 and M10 must not appear
  expect_false(any(c(paste0("M", n_mar + 1), paste0("M", n_mar + 2)) %in%
                     all_markers))
})

# =============================================================================
# primePanel() — chromosome names with dots produce a warning
# =============================================================================

test_that("primePanel warns when chromosome names contain dots", {
  set.seed(30)
  n_lines <- 10; n_mar <- 4
  ids  <- paste0("L", seq_len(n_lines))
  gmat <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                 nrow = n_lines,
                 dimnames = list(ids, paste0("M", seq_len(n_mar))))
  # Use a chromosome name with a dot
  map_df <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = "Chr.1",            # dot in chromosome name
    pos    = seq_len(n_mar),
    stringsAsFactors = FALSE
  )

  expect_warning(
    suppressMessages(primePanel(gmat, map_df, encoding = "012", maf = 0)),
    regexp = "[Dd]ots? found"
  )
})

# =============================================================================
# primePanel() — $geno slot structure: $data, $map, $imputed.data per chr
# =============================================================================

test_that("primePanel $geno each chromosome has $data, $map, $imputed.data", {
  set.seed(31)
  n_lines <- 15; n_mar <- 10
  ids     <- paste0("L", seq_len(n_lines))
  gmat    <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                    nrow = n_lines,
                    dimnames = list(ids, paste0("M", seq_len(n_mar))))
  map_df  <- data.frame(
    marker = paste0("M", seq_len(n_mar)),
    chr    = rep(c("A", "B"), each = n_mar / 2),
    pos    = rep(seq(5, 50, length.out = n_mar / 2), times = 2),
    stringsAsFactors = FALSE
  )

  expect_message(
    panel <- primePanel(gmat, map_df, encoding = "012", maf = 0),
    regexp = "wgPanel object created"
  )

  for (chr in names(panel$geno)) {
    chr_el <- panel$geno[[chr]]
    expect_true(!is.null(chr_el$data),         info = paste("$data missing, chr", chr))
    expect_true(!is.null(chr_el$map),          info = paste("$map missing,  chr", chr))
    expect_true(!is.null(chr_el$imputed.data), info = paste("$imputed.data missing, chr", chr))
    expect_true(is.matrix(chr_el$data),        info = paste("$data not matrix, chr", chr))
    expect_true(is.matrix(chr_el$imputed.data),info = paste("$imputed.data not matrix, chr", chr))
  }
})

# =============================================================================
# primePanel() — line IDs preserved as rownames in imputed.data
# =============================================================================

test_that("primePanel preserves line IDs as rownames in imputed.data", {
  set.seed(32)
  n_lines <- 12; n_mar <- 6
  ids     <- paste0("Line_", LETTERS[seq_len(n_lines)])
  gmat    <- matrix(sample(0:2, n_lines * n_mar, replace = TRUE),
                    nrow = n_lines,
                    dimnames = list(ids, paste0("SNP", seq_len(n_mar))))
  map_df  <- data.frame(
    marker = paste0("SNP", seq_len(n_mar)),
    chr    = "ChrAlpha",
    pos    = seq(1, 60, length.out = n_mar),
    stringsAsFactors = FALSE
  )

  expect_message(
    panel <- primePanel(gmat, map_df, encoding = "012", maf = 0),
    regexp = "wgPanel object created"
  )

  for (chr in names(panel$geno)) {
    expect_equal(rownames(panel$geno[[chr]]$imputed.data), ids,
                 info = paste("rownames mismatch in chr", chr))
  }
})
