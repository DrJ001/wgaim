# =============================================================================
# primePanel.R
# Convert association panel marker data into a 'panel'/'interval' class object
# for use with gwasAim() and gpAim().
#
# The resulting object has the same internal structure as an 'interval' object
# from primeCross(), so it feeds directly into the shared wgAim engine pieces
# without modification.
# =============================================================================

#' Prepare Marker Panel Data for Whole-Genome Analysis for gwasAim and gpAim
#'
#' @description
#' Converts a matrix or data frame of marker genotypes and an associated
#' genetic map into a \code{"wgPanel"} object suitable for use with
#' \code{\link{gwasAim}} and \code{\link{gpAim}}.
#'
#' @param geno A numeric matrix (\code{lines x markers}) with row names
#'   identifying lines, or a \code{data.frame} with a column named by
#'   \code{id} identifying lines and remaining columns containing marker
#'   genotypes.
#'
#'   \strong{Genotype values:} The function accepts both integer-coded
#'   genotypes (0, 1, 2) and fractional \emph{dosage} values in the range
#'   \eqn{[0, 2]} output by imputation software such as Beagle, IMPUTE2,
#'   or Minimac. Under \code{encoding = "012"}, the transformation
#'   \eqn{x \mapsto x - 1} is applied to all values whether or not they
#'   are integers, yielding dosage values in \eqn{[-1, 1]}. A pre-imputed
#'   dosage matrix therefore requires no further imputation inside
#'   \code{primePanel}.
#'
#' @param map A \code{data.frame} containing the genetic map with at minimum
#'   three columns: marker names (\code{map.id}), chromosome identifiers
#'   (\code{map.chr}), and cM positions (\code{map.pos}). Markers present
#'   in \code{geno} but absent from \code{map} are dropped with a message.
#'   Markers are sorted by position within each chromosome.
#' @param id Character string naming the line identifier column in \code{geno}
#'   when \code{geno} is a \code{data.frame}. Ignored when \code{geno} is a
#'   matrix (row names are used). Default \code{"id"}.
#' @param map.id Character string naming the marker column in \code{map}.
#'   Default \code{"marker"}.
#' @param map.chr Character string naming the chromosome column in \code{map}.
#'   Default \code{"chr"}. Chromosome names containing dots are replaced with
#'   underscores with a warning, as dots are used as field separators in
#'   internal marker naming.
#' @param map.pos Character string naming the cM position column in \code{map}.
#'   Default \code{"pos"}.
#' @param encoding Character string specifying the genotype encoding in
#'   \code{geno}:
#'   \describe{
#'     \item{\code{"012"} (default)}{Allele-count coding where 0 = AA
#'       homozygote, 1 = AB heterozygote, 2 = BB homozygote. Values are
#'       shifted by \eqn{-1} to produce additive \eqn{\pm 1} coding. This
#'       also accepts fractional dosage values in \eqn{[0, 2]} output by
#'       imputation software -- no additional handling is required.}
#'     \item{\code{"pm1"}}{Genotypes already in additive \eqn{\pm 1}
#'       coding. No transformation is applied. Fractional values in
#'       \eqn{[-1, 1]} are accepted.}
#'   }
#' @param impute Character string controlling how missing values (\code{NA})
#'   are handled \emph{after} encoding:
#'   \describe{
#'     \item{\code{"none"} (default)}{No imputation is performed. If any
#'       \code{NA}s remain, \code{primePanel} stops with an informative error
#'       recommending pre-imputation with dedicated software. This is the
#'       recommended path for large panels.}
#'     \item{\code{"knn"}}{Chromosome-wise \emph{k}-nearest neighbour
#'       imputation.  For each chromosome the pairwise distance matrix is
#'       computed via \code{tcrossprod()} (BLAS), the neighbour ranking is
#'       pre-sorted once, and each missing value is imputed from the mean of
#'       the \code{knn.k} nearest lines that carry a non-missing call at that
#'       marker.  Imputation respects linkage disequilibrium within each
#'       chromosome and is the recommended in-package option when pre-imputation
#'       with dedicated software is not feasible.}
#'     \item{\code{"mean"}}{Missing values are replaced with the column mean
#'       (mean genotype across lines at each marker). This is a crude
#'       approximation that ignores linkage disequilibrium and haplotype
#'       structure; it is provided only as a convenience for small panels or
#'       exploratory analyses. A warning is always issued when this option
#'       is used.}
#'   }
#' @param knn.k Positive integer. Number of nearest neighbours used when
#'   \code{impute = "knn"}.  Default \code{5}.  Larger values produce
#'   smoother imputations but are slower; values between 3 and 10 are
#'   typical.  Ignored when \code{impute != "knn"}.
#'
#' @details
#' \strong{Data quality and filtering:}
#'
#' \code{primePanel} assumes the data are already clean and ready for
#' analysis.  No MAF filtering, missingness filtering, or duplicate removal
#' is performed.  Use \code{\link{checkPanel}} to diagnose data quality
#' issues and \code{\link{filterPanel}} to apply filters before calling
#' \code{primePanel}.  Alternatively, the \code{filteredPanel} S3 method
#' accepts the output of \code{\link{filterPanel}} directly.
#'
#' \strong{Missing data and imputation:}
#'
#' For real association panels, missing genotypes should be imputed using
#' dedicated software \emph{before} calling \code{primePanel}. Tools such as
#' Beagle (\url{https://faculty.washington.edu/browning/beagle/beagle.html}),
#' IMPUTE2, or Minimac use haplotype reference panels and linkage
#' disequilibrium information to produce well-calibrated dosage values.
#' Their output can be supplied directly to \code{primePanel} as a matrix
#' of dosage values in \eqn{[0, 2]} with \code{encoding = "012"} and
#' \code{impute = "none"}.
#'
#' Column-mean imputation (\code{impute = "mean"}) inflates
#' apparent precision and ignores the LD structure that makes imputation
#' informative. For a panel of \eqn{n} lines and \eqn{p} markers it requires
#' \eqn{O(np)} operations and is fast, but the quality of imputed values
#' degrades rapidly as missingness increases. It is not recommended for
#' missing rates above 1--2\%.
#'
#' \strong{Encoding and dosage values:}
#'
#' The \eqn{x \mapsto x - 1} transformation under \code{encoding = "012"}
#' is exact for integers and valid for fractional dosages. A dosage of 0.8
#' (likely heterozygous allele) becomes \eqn{-0.2}; a dosage of 1.6
#' (likely BB) becomes \eqn{+0.6}. Values are not rounded or thresholded;
#' the full dosage information is preserved for use in the relationship
#' matrix or marker effect model.
#'
#' @return An object of class \code{c("wgPanel", "wgCross")} -- a list with:
#' \describe{
#'   \item{\code{$pheno}}{A \code{data.frame} with the line identifier column.}
#'   \item{\code{$geno}}{A named list (one element per chromosome) each
#'     containing:
#'     \describe{
#'       \item{\code{$data}}{Encoded genotype matrix (lines x markers).}
#'       \item{\code{$map}}{Named numeric vector of cM positions.}
#'       \item{\code{$imputed.data}}{The analysis-ready genotype matrix
#'         (same as \code{$data} after any imputation) in additive
#'         \eqn{\pm 1} or dosage coding.}
#'     }
#'   }
#' }
#'
#' @seealso \code{\link{gwasAim}}, \code{\link{gpAim}},
#'   \code{\link{primeCross}}
#'
#' @examples
#' \dontrun{
#' # --- Recommended path: check, filter, then prime with KNN imputation ---
#' chk   <- checkPanel(geno.mat, map.df)
#' print(chk)
#' clean <- filterPanel(chk, maf = 0.05, miss.line = 0.20, miss.marker = 0.20)
#' panel <- primePanel(clean, impute = "knn", knn.k = 5)
#'
#' # --- Direct path (data already clean, no imputation needed) ---
#' panel <- primePanel(geno.mat, map.df, encoding = "012", impute = "none")
#'
#' # --- KNN imputation on raw data ---
#' panel <- primePanel(geno.mat, map.df, encoding = "012",
#'                     impute = "knn", knn.k = 5)
#'
#' # --- Mean imputation for quick exploratory analysis ---
#' panel <- primePanel(geno.mat, map.df, encoding = "012", impute = "mean")
#' }
#'
#' @export
primePanel <- function(geno, map,
                       id      = "id",
                       map.id  = "marker",
                       map.chr = "chr",
                       map.pos = "pos",
                       encoding = "012",
                       impute  = "none",
                       knn.k   = 5L) {
    # Accept a filteredPanel object as the first argument; extract stored fields.
    if (inherits(geno, "filteredPanel")) {
        fp   <- geno
        mat  <- fp$geno
        if (is.null(rownames(mat)))
            stop("filteredPanel object has no row names on $geno.")
        return(.primePanel_core(geno     = mat,
                                map      = fp$map,
                                id       = fp$id,
                                map.id   = fp$map.id,
                                map.chr  = fp$map.chr,
                                map.pos  = fp$map.pos,
                                encoding = fp$encoding,
                                impute   = impute,
                                knn.k    = knn.k))
    }
    .primePanel_core(geno = geno, map = map, id = id, map.id = map.id,
                     map.chr = map.chr, map.pos = map.pos,
                     encoding = encoding, impute = impute, knn.k = knn.k)
}

# =============================================================================
# Internal workhorse (formerly the body of primePanel())
# =============================================================================
.primePanel_core <- function(geno, map, id, map.id, map.chr, map.pos,
                              encoding, impute, knn.k = 5L) {

    impute <- match.arg(impute, c("none", "mean", "knn"))
    knn.k  <- as.integer(knn.k)
    if (impute == "knn" && (length(knn.k) != 1L || knn.k < 1L))
        stop("knn.k must be a positive integer.")

    # -------------------------------------------------------------------------
    # 1. Extract line IDs and raw genotype matrix
    # -------------------------------------------------------------------------
    if (is.data.frame(geno)) {
        if (!id %in% names(geno))
            stop("Column '", id, "' not found in geno data frame.")
        ids      <- as.character(geno[[id]])
        geno.mat <- as.matrix(geno[, !names(geno) %in% id, drop = FALSE])
    } else {
        geno.mat <- as.matrix(geno)
        ids      <- rownames(geno.mat)
        if (is.null(ids))
            stop("geno matrix must have row names identifying lines.")
    }
    storage.mode(geno.mat) <- "numeric"

    # -------------------------------------------------------------------------
    # 2. Validate map columns
    # -------------------------------------------------------------------------
    for (col in c(map.id, map.chr, map.pos)) {
        if (!col %in% names(map))
            stop("Column '", col, "' not found in map.")
    }

    # -------------------------------------------------------------------------
    # 3. Align markers between geno and map
    # -------------------------------------------------------------------------
    map.markers  <- as.character(map[[map.id]])
    geno.markers <- colnames(geno.mat)
    if (is.null(geno.markers))
        stop("geno must have column names matching marker names in map.")
    common <- intersect(map.markers, geno.markers)
    if (length(common) == 0)
        stop("No common markers found between geno columns and map '",
             map.id, "' column.")
    n.drop <- length(geno.markers) - length(common)
    if (n.drop > 0)
        message(n.drop, " marker(s) in geno not found in map and will be dropped.")
    geno.mat <- geno.mat[, common, drop = FALSE]
    map      <- map[map[[map.id]] %in% common, , drop = FALSE]

    # -------------------------------------------------------------------------
    # 4. Encode genotypes to additive +/-1 / dosage coding
    # -------------------------------------------------------------------------
    if (encoding == "012") {
        # Shift by -1: integers 0/1/2 -> -1/0/+1
        # Fractional dosage values in [0,2] -> [-1,+1]; no rounding applied.
        geno.mat <- geno.mat - 1
        # Soft range check -- values outside [-1.05, 1.05] suggest wrong encoding
        rng <- range(geno.mat, na.rm = TRUE)
        if (rng[1] < -1.05 || rng[2] > 1.05)
            warning("After encoding, some values lie outside [-1, 1] ",
                    sprintf("(range: [%.3f, %.3f]). ", rng[1], rng[2]),
                    "Check that encoding = '012' is correct for your data.")
    } else if (encoding == "pm1") {
        rng <- range(geno.mat, na.rm = TRUE)
        if (rng[1] < -1.05 || rng[2] > 1.05)
            warning("With encoding = 'pm1', some values lie outside [-1, 1] ",
                    sprintf("(range: [%.3f, %.3f]). ", rng[1], rng[2]),
                    "Verify that genotypes are correctly pre-coded.")
    } else {
        stop("encoding must be '012' (allele counts or dosage) or 'pm1' ",
             "(already in additive \u00b11 coding).")
    }

    # -------------------------------------------------------------------------
    # 5. Handle missing values
    # -------------------------------------------------------------------------
    n.missing <- sum(is.na(geno.mat))
    if (n.missing > 0) {
        pct <- round(100 * n.missing / length(geno.mat), 2)
        if (impute == "none") {
            stop(
                n.missing, " missing genotype(s) remain (", pct, "% of values). ",
                "primePanel does not impute by default.\n\n",
                "Options:\n",
                "  1. Pre-impute using dedicated software (recommended):\n",
                "       Beagle  : https://faculty.washington.edu/browning/beagle/\n",
                "       IMPUTE2 : https://mathgen.stats.ox.ac.uk/impute/\n",
                "       Minimac : https://genome.sph.umich.edu/wiki/Minimac\n",
                "     Supply the resulting dosage matrix with encoding = '012'.\n\n",
                "  2. Chromosome-wise KNN imputation (recommended in-package option):\n",
                "       primePanel(..., impute = 'knn', knn.k = 5)\n\n",
                "  3. Column-mean imputation for small panels / exploration:\n",
                "       primePanel(..., impute = 'mean')\n",
                "     Warning: mean imputation ignores LD and haplotype structure."
            )
        } else if (impute == "mean") {
            warning(
                "Imputing ", n.missing, " missing value(s) with column means.\n",
                "Column-mean imputation ignores LD and haplotype structure. ",
                "Results may be biased, especially at high missing rates. ",
                "For publication-quality analyses, use dedicated imputation software ",
                "(Beagle, IMPUTE2, Minimac) prior to calling primePanel()."
            )
            col.means <- colMeans(geno.mat, na.rm = TRUE)
            for (j in seq_len(ncol(geno.mat))) {
                mis <- is.na(geno.mat[, j])
                if (any(mis))
                    geno.mat[mis, j] <- col.means[j]
            }
        } else {
            # impute = "knn" -- chromosome-wise; deferred to step 7
            message(sprintf(
                "KNN imputation (k = %d): %d missing value(s) (%.2f%%) ",
                knn.k, n.missing, pct),
                "will be imputed within each chromosome.")
        }
    }
    rownames(geno.mat) <- ids

    # -------------------------------------------------------------------------
    # 6. Sanitise chromosome names (dots are field separators internally)
    # -------------------------------------------------------------------------
    chr.vals <- as.character(map[[map.chr]])
    if (any(grepl("\\.", chr.vals))) {
        warning("Dots found in chromosome names: replacing with underscores ",
                "to avoid conflicts with internal naming conventions.")
        chr.vals       <- gsub("\\.", "_", chr.vals)
        map[[map.chr]] <- chr.vals
    }

    # -------------------------------------------------------------------------
    # 7. Split by chromosome, sort by position
    # -------------------------------------------------------------------------
    chrs      <- unique(chr.vals)
    geno.list <- lapply(chrs, function(ch) {
        chr.map <- map[map[[map.chr]] == ch, , drop = FALSE]
        chr.map <- chr.map[order(as.numeric(chr.map[[map.pos]])), ]
        mnames  <- as.character(chr.map[[map.id]])
        chr.g   <- geno.mat[, mnames, drop = FALSE]
        chr.pos <- setNames(as.numeric(chr.map[[map.pos]]), mnames)
        chr.imp <- if (impute == "knn" && anyNA(chr.g))
            .knn_impute_chr(chr.g, k = knn.k)
        else
            chr.g
        list(
            data         = chr.g,
            map          = chr.pos,
            imputed.data = chr.imp
        )
    })
    names(geno.list) <- chrs

    # -------------------------------------------------------------------------
    # 8. Assemble and return
    # -------------------------------------------------------------------------
    panel <- list(
        pheno = setNames(data.frame(ids, stringsAsFactors = FALSE), id),
        geno  = geno.list
    )
    class(panel) <- c("wgPanel", "wgCross", "cross")

    n.final <- sum(sapply(geno.list, function(el) ncol(el$imputed.data)))
    message("wgPanel object created: ", length(ids), " lines, ",
            n.final, " markers across ", length(chrs), " chromosomes.")
    panel
}

# =============================================================================
# .knn_impute_chr  -- vectorised KNN imputation for a single chromosome matrix
#
# Arguments:
#   geno_chr : numeric matrix (lines x markers), encoded to +/-1, may have NAs
#   k        : number of nearest neighbours
#
# Algorithm:
#   1. Mean-substitute NAs to build a complete matrix for distance calculation.
#   2. Compute the n x n pairwise distance matrix via tcrossprod() (BLAS).
#   3. Pre-sort the neighbour matrix once (n x n; column i = line i's
#      neighbours ranked by distance).
#   4. For each marker with missing values, use vectorised logical indexing
#      on the pre-sorted matrix to identify valid (non-missing) neighbours,
#      then impute as their mean.
#
# Complexity:
#   O(n^2 * p) for the BLAS tcrossprod  (fast in practice)
#   O(n^2 log n) for the one-time presort
#   O(n_miss * n) for the imputation sweep  (n_miss = total missing values)
# =============================================================================
.knn_impute_chr <- function(geno_chr, k = 5L) {
    na_cols <- which(colSums(is.na(geno_chr)) > 0L)
    if (!length(na_cols)) return(geno_chr)

    n         <- nrow(geno_chr)
    col_means <- colMeans(geno_chr, na.rm = TRUE)

    # Step 1: mean-substitute for distance calculation
    geno_sub <- geno_chr
    for (j in na_cols)
        geno_sub[is.na(geno_sub[, j]), j] <- col_means[j]

    # Step 2: pairwise distance via tcrossprod (BLAS inner-product formula)
    #   dist^2(i,j) = ||x_i||^2 + ||x_j||^2 - 2 <x_i, x_j>
    row_ss  <- rowSums(geno_sub^2)
    cp      <- tcrossprod(geno_sub)
    dist_sq <- outer(row_ss, row_ss, "+") - 2 * cp
    d_mat   <- sqrt(pmax(dist_sq, 0))
    diag(d_mat) <- 0

    # Step 3: pre-sort neighbour matrix once
    #   nn_rank[r, i] = index of r-th nearest neighbour of line i
    nn_rank <- apply(d_mat, 1L, order)

    # Step 4: impute each marker with missing values
    result <- geno_chr
    for (j in na_cols) {
        miss_j    <- which( is.na(geno_chr[, j]))
        nonmiss_j <- which(!is.na(geno_chr[, j]))
        if (!length(nonmiss_j)) {
            result[miss_j, j] <- col_means[j]   # all missing: fall back to mean
            next
        }

        # Logical flag: which neighbour slots are valid at this marker
        nonmiss_flag          <- logical(n)
        nonmiss_flag[nonmiss_j] <- TRUE

        # Sub-matrix of pre-sorted neighbours for missing lines only
        nn_sub    <- nn_rank[, miss_j, drop = FALSE]         # n x n_miss
        valid_mat <- matrix(nonmiss_flag[nn_sub],
                            nrow = nrow(nn_sub))             # n x n_miss logical

        # Impute each missing line from its k nearest valid neighbours
        result[miss_j, j] <- vapply(seq_along(miss_j), function(idx) {
            vr <- which(valid_mat[, idx])
            if (!length(vr)) return(col_means[j])
            nn_idx <- nn_sub[vr[seq_len(min(k, length(vr)))], idx]
            mean(geno_chr[nn_idx, j])
        }, numeric(1L))
    }
    result
}


