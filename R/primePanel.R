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
#'       imputation software â€” no additional handling is required.}
#'     \item{\code{"pm1"}}{Genotypes already in additive \eqn{\pm 1}
#'       coding. No transformation is applied. Fractional values in
#'       \eqn{[-1, 1]} are accepted.}
#'   }
#' @param maf Numeric scalar. Markers with minor allele frequency below this
#'   threshold are removed. Computed from column means after encoding.
#'   Default \code{0.05}. Set to \code{0} or \code{NULL} to disable.
#' @param impute Character string controlling how missing values (\code{NA})
#'   are handled \emph{after} encoding and MAF filtering:
#'   \describe{
#'     \item{\code{"none"} (default)}{No imputation is performed. If any
#'       \code{NA}s remain, \code{primePanel} stops with an informative error
#'       recommending pre-imputation with dedicated software. This is the
#'       recommended path for large panels.}
#'     \item{\code{"mean"}}{Missing values are replaced with the column mean
#'       (mean genotype across lines at each marker). This is a crude
#'       approximation that ignores linkage disequilibrium and haplotype
#'       structure; it is provided only as a convenience for small panels or
#'       exploratory analyses. A warning is always issued when this option
#'       is used.}
#'   }
#'
#' @details
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
#'   \code{\link[wgAim]{primeCross}}
#'
#' @examples
#' \dontrun{
#' # --- Integer 0/1/2 genotypes, no missing data ---
#' geno.mat <- matrix(sample(0:2, 200 * 500, replace = TRUE),
#'                    nrow = 200, ncol = 500)
#' rownames(geno.mat) <- paste0("Line_", 1:200)
#' colnames(geno.mat) <- paste0("M",     1:500)
#' map.df <- data.frame(
#'   marker = paste0("M", 1:500),
#'   chr    = rep(paste0("Chr", 1:5), each = 100),
#'   pos    = rep(seq(1, 100), times = 5))
#'
#' panel <- primePanel(geno.mat, map.df, encoding = "012", maf = 0.05)
#'
#' # --- Dosage values from Beagle (fractional, no NAs) ---
#' # dosage.mat is a matrix of values in [0, 2] from Beagle output
#' panel <- primePanel(dosage.mat, map.df, encoding = "012", impute = "none")
#'
#' # --- Small panel with some missingness (mean imputation, with warning) ---
#' panel <- primePanel(geno.mat, map.df, encoding = "012", impute = "mean")
#' }
#'
#' @export
primePanel <- function(geno, map, id = "id",
                      map.id   = "marker",
                      map.chr  = "chr",
                      map.pos  = "pos",
                      encoding = "012",
                      maf      = 0.05,
                      impute   = "none") {

    impute <- match.arg(impute, c("none", "mean"))

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
    # 4. Encode genotypes to additive Â±1 / dosage coding
    # -------------------------------------------------------------------------
    if (encoding == "012") {
        # Shift by -1: integers 0/1/2 -> -1/0/+1
        # Fractional dosage values in [0,2] -> [-1,+1]; no rounding applied.
        geno.mat <- geno.mat - 1
        # Soft range check â€” values outside [-1.05, 1.05] suggest wrong encoding
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
    # 5. MAF filter
    # -------------------------------------------------------------------------
    if (!is.null(maf) && maf > 0) {
        # In Â±1 coding: mean = 2p - 1, so p = (mean + 1) / 2
        allele.freq <- colMeans(geno.mat, na.rm = TRUE)
        p           <- (allele.freq + 1) / 2
        minor.freq  <- pmin(p, 1 - p)
        keep        <- minor.freq >= maf
        n.filt      <- sum(!keep)
        if (n.filt > 0)
            message(n.filt, " marker(s) removed with MAF < ", maf, ".")
        geno.mat <- geno.mat[, keep, drop = FALSE]
        map      <- map[map[[map.id]] %in% colnames(geno.mat), , drop = FALSE]
    }

    # -------------------------------------------------------------------------
    # 6. Handle missing values
    # -------------------------------------------------------------------------
    n.missing <- sum(is.na(geno.mat))
    if (n.missing > 0) {
        if (impute == "none") {
            pct <- round(100 * n.missing / length(geno.mat), 2)
            stop(
                n.missing, " missing genotype(s) remain (", pct, "% of values). ",
                "primePanel does not impute by default.\n\n",
                "Options:\n",
                "  1. Pre-impute using dedicated software (recommended):\n",
                "       Beagle  : https://faculty.washington.edu/browning/beagle/\n",
                "       IMPUTE2 : https://mathgen.stats.ox.ac.uk/impute/\n",
                "       Minimac : https://genome.sph.umich.edu/wiki/Minimac\n",
                "     Supply the resulting dosage matrix with encoding = '012'.\n\n",
                "  2. Use column-mean imputation for small panels / exploration:\n",
                "       primePanel(..., impute = 'mean')\n",
                "     Warning: mean imputation ignores LD and haplotype structure.\n",
                "     Not recommended for missing rates above 1-2%."
            )
        } else {
            # impute = "mean"
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
        }
    }
    rownames(geno.mat) <- ids

    # -------------------------------------------------------------------------
    # 7. Sanitise chromosome names (dots are field separators internally)
    # -------------------------------------------------------------------------
    chr.vals <- as.character(map[[map.chr]])
    if (any(grepl("\\.", chr.vals))) {
        warning("Dots found in chromosome names: replacing with underscores ",
                "to avoid conflicts with internal naming conventions.")
        chr.vals       <- gsub("\\.", "_", chr.vals)
        map[[map.chr]] <- chr.vals
    }

    # -------------------------------------------------------------------------
    # 8. Split by chromosome, sort by position
    # -------------------------------------------------------------------------
    chrs      <- unique(chr.vals)
    geno.list <- lapply(chrs, function(ch) {
        chr.map <- map[map[[map.chr]] == ch, , drop = FALSE]
        chr.map <- chr.map[order(as.numeric(chr.map[[map.pos]])), ]
        mnames  <- as.character(chr.map[[map.id]])
        chr.g   <- geno.mat[, mnames, drop = FALSE]
        chr.pos <- setNames(as.numeric(chr.map[[map.pos]]), mnames)
        list(
            data         = chr.g,
            map          = chr.pos,
            imputed.data = chr.g   # analysis-ready; marker mode uses this
        )
    })
    names(geno.list) <- chrs

    # -------------------------------------------------------------------------
    # 9. Assemble and return
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


