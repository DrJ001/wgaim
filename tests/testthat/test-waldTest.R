# =============================================================================
# test-waldTest.R
# Additional tests for waldTest() targeting the branches not yet covered.
#
# Existing coverage in test-multivariate.R covers the zero-equality path
# (type = "zero") in waldTest.asreml via .packResults mock tests.
#
# This file adds:
#   - waldTest.default guard
#   - waldTest.asreml type = "con" (treatment contrast) paths:
#       - single contrast vector
#       - contrast matrix (multi-row)
#       - default Helmert contrast (no comp supplied)
#       - group label overrides
#       - keep.fac = FALSE
#   - waldTest.asreml type = "zero" path (direct call)
#   - waldTest.asreml mixed cc list (both "con" and "zero")
#   - Error: invalid type in cc
#   - Error: subscript out of bounds
#   - Error: coefficient names not matched
#   - Error: duplicate zero-test group names
# =============================================================================

# ---- Build a minimal mock asreml model that satisfies waldTest.asreml -------
# waldTest.asreml needs:
#   object$Cfixed      - coefficient variance matrix  (nc x nc)
#   object$coefficients$fixed - named (nc x 1) matrix of fixed effects
#   object$sigma2      - residual variance scalar
.make_wald_mock <- function(nc = 4, seed = 42) {
    set.seed(seed)
    # Random positive-definite Cfixed
    A       <- matrix(rnorm(nc * nc), nc, nc)
    Cfixed  <- crossprod(A) / nc + diag(nc) * 0.1
    nms     <- c("(Intercept)", paste0("X", seq_len(nc - 1)))
    dimnames(Cfixed) <- list(nms, nms)
    coef_fixed <- matrix(rnorm(nc), nc, 1, dimnames = list(nms, "effect"))
    obj <- list(
        Cfixed          = Cfixed,
        sigma2          = 1.0,
        coefficients    = list(fixed = coef_fixed),
        converge        = TRUE
    )
    class(obj) <- "asreml"
    obj
}

# =============================================================================
# 1.  waldTest.default guard
# =============================================================================

test_that("waldTest.default: non-asreml object stops with error", {
    expect_error(waldTest(list()), regexp = "\"asreml\"")
})

test_that("waldTest.default: data.frame stops with error", {
    expect_error(waldTest(data.frame(x = 1)), regexp = "\"asreml\"")
})

# =============================================================================
# 2.  waldTest.asreml -- type = "zero" (zero-equality test, direct call)
# =============================================================================

test_that("waldTest.asreml zero: returns list with $Zero data frame", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(
        grp1 = list(coef = 2:3, type = "zero")
    )
    result <- waldTest(obj, cc)
    expect_type(result, "list")
    expect_s3_class(result$Zero, "data.frame")
})

test_that("waldTest.asreml zero: $Contrasts is NULL when only zero tests", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(grp1 = list(coef = 2:3, type = "zero"))
    result <- waldTest(obj, cc)
    expect_null(result$Contrasts)
})

test_that("waldTest.asreml zero: Wald statistic is non-negative", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(grp1 = list(coef = 2:4, type = "zero"))
    result <- waldTest(obj, cc)
    expect_gte(result$Zero[["Wald Statistic"]], 0)
})

test_that("waldTest.asreml zero: p-value in [0, 1]", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(grp1 = list(coef = 2:3, type = "zero"))
    result <- waldTest(obj, cc)
    pval <- result$Zero[["P-Value"]]
    expect_gte(pval, 0)
    expect_lte(pval, 1)
})

test_that("waldTest.asreml zero: group label used as row name", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(mygroup = list(coef = 2:3, type = "zero", group = "MyQTL"))
    result <- waldTest(obj, cc)
    expect_equal(rownames(result$Zero), "MyQTL")
})

test_that("waldTest.asreml zero: coefficient names resolved by name string", {
    obj <- .make_wald_mock(nc = 4)
    nms <- rownames(obj$coefficients$fixed)
    cc  <- list(grp1 = list(coef = nms[2:3], type = "zero"))
    result <- waldTest(obj, cc)
    expect_s3_class(result$Zero, "data.frame")
})

test_that("waldTest.asreml zero: multiple zero groups produce multiple rows", {
    obj <- .make_wald_mock(nc = 5)
    cc  <- list(
        g1 = list(coef = 2L,   type = "zero", group = "Q1"),
        g2 = list(coef = 3:4L, type = "zero", group = "Q2")
    )
    result <- waldTest(obj, cc)
    expect_equal(nrow(result$Zero), 2L)
})

# =============================================================================
# 3.  waldTest.asreml -- type = "con" (treatment contrast)
# =============================================================================

test_that("waldTest.asreml con: single contrast vector returns $Contrasts data frame", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(
        c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1))
    )
    result <- waldTest(obj, cc)
    expect_s3_class(result$Contrasts, "data.frame")
})

test_that("waldTest.asreml con: $Zero is NULL when only contrast tests", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1)))
    result <- waldTest(obj, cc)
    expect_null(result$Zero)
})

test_that("waldTest.asreml con: Wald statistic is non-negative", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1)))
    result <- waldTest(obj, cc)
    expect_gte(result$Contrasts[["Wald Statistic"]], 0)
})

test_that("waldTest.asreml con: p-value in [0, 1]", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1)))
    result <- waldTest(obj, cc)
    pval <- result$Contrasts[["P-Value"]]
    expect_gte(pval, 0)
    expect_lte(pval, 1)
})

test_that("waldTest.asreml con: contrast matrix (multi-row) produces multiple rows", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(
        c1 = list(
            coef = c(2L, 3L, 4L),
            type = "con",
            comp = matrix(c(1, -1, 0,
                            0,  1, -1), nrow = 2, byrow = TRUE)
        )
    )
    result <- waldTest(obj, cc)
    expect_equal(nrow(result$Contrasts), 2L)
})

test_that("waldTest.asreml con: group labels override row names", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(
        c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1),
                  group = list(left = "TrtA", right = "TrtB"))
    )
    result <- waldTest(obj, cc)
    expect_true(grepl("TrtA", rownames(result$Contrasts)[1]))
})

test_that("waldTest.asreml con: keep.fac = FALSE strips factor prefixes", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1)))
    # Should not error; keep.fac just affects label formatting
    result <- waldTest(obj, cc, keep.fac = FALSE)
    expect_s3_class(result$Contrasts, "data.frame")
})

test_that("waldTest.asreml con: default Helmert contrast (no comp) runs with message", {
    obj <- .make_wald_mock(nc = 4)
    # No comp supplied => Helmert default applied; function calls cat() with warning
    cc  <- list(c1 = list(coef = c(2L, 3L), type = "con"))
    expect_no_error({
        out <- capture.output(result <- waldTest(obj, cc))
    })
    expect_s3_class(result$Contrasts, "data.frame")
})

# =============================================================================
# 4.  waldTest.asreml -- mixed cc (both "con" and "zero")
# =============================================================================

test_that("waldTest.asreml mixed: both $Contrasts and $Zero populated", {
    obj <- .make_wald_mock(nc = 5)
    cc  <- list(
        con1  = list(coef = c(2L, 3L), type = "con", comp = c(1, -1)),
        zero1 = list(coef = 4:5L,      type = "zero", group = "Grp1")
    )
    result <- waldTest(obj, cc)
    expect_s3_class(result$Contrasts, "data.frame")
    expect_s3_class(result$Zero,      "data.frame")
})

# =============================================================================
# 5.  waldTest.asreml -- error guards
# =============================================================================

test_that("waldTest.asreml: invalid cc element names stops with error", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = 2L, type = "con", BADFIELD = 1))
    expect_error(waldTest(obj, cc), regexp = "names")
})

test_that("waldTest.asreml: subscript out of bounds stops with error", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = 99L, type = "zero"))
    expect_error(waldTest(obj, cc), regexp = "subscript|bounds")
})

test_that("waldTest.asreml: unmatched coefficient name stops with error", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = "DOES_NOT_EXIST", type = "zero"))
    expect_error(waldTest(obj, cc), regexp = "Names in contrast|match")
})

test_that("waldTest.asreml: invalid type string stops with error", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = 2L, type = "bad_type"))
    expect_error(waldTest(obj, cc), regexp = "\"con\".*\"zero\"|\"zero\".*\"con\"")
})

test_that("waldTest.asreml: duplicate zero-group names stops with error", {
    obj <- .make_wald_mock(nc = 5)
    cc  <- list(
        g1 = list(coef = 2L, type = "zero", group = "SameName"),
        g2 = list(coef = 3L, type = "zero", group = "SameName")
    )
    expect_error(waldTest(obj, cc), regexp = "[Dd]uplicate")
})

test_that("waldTest.asreml: contrast vector length mismatch stops with error", {
    obj <- .make_wald_mock(nc = 4)
    cc  <- list(c1 = list(coef = c(2L, 3L), type = "con", comp = c(1, -1, 0)))
    expect_error(waldTest(obj, cc), regexp = "[Ll]ength")
})

test_that("waldTest.asreml: contrast matrix column mismatch stops with error", {
    obj <- .make_wald_mock(nc = 4)
    # 2 coef but 3-column matrix
    cc  <- list(c1 = list(
        coef = c(2L, 3L), type = "con",
        comp = matrix(c(1, 0, -1), nrow = 1)
    ))
    expect_error(waldTest(obj, cc), regexp = "[Cc]olumn|mismatch|match")
})
