# =============================================================================
# makePanel.R
# Convert association panel marker data into a 'panel'/'interval' class object
# for use with GWASAim().
#
# The resulting object has the same internal structure as an 'interval' object
# from cross2int(), so it feeds directly into the shared wgAim engine pieces
# without modification.
# =============================================================================

#' Convert Marker Panel Data to a Panel Object for GWASAim and GPAim
#'
#' @description
#' Converts a matrix or data frame of marker genotypes and an associated
#' genetic map into a \code{c("panel","interval")} object suitable for use with
#' \code{\link{GWASAim}} and \code{\link{GPAim}}. The function handles genotype
#' encoding, minor allele frequency filtering, mean imputation of missing values,
#' and chromosome name sanitisation.
#'
#' @param geno A numeric matrix (\code{lines x markers}) with row names
#'   identifying lines, or a \code{data.frame} with a column named by the
#'   \code{id} argument identifying lines and remaining columns containing
#'   marker genotypes. Marker genotypes should be encoded as specified by
#'   \code{encoding}.
#' @param map A \code{data.frame} containing the genetic map. Must contain
#'   at minimum three columns: marker names (see \code{map.id}), chromosome
#'   identifiers (see \code{map.chr}), and cM positions (see \code{map.pos}).
#'   Markers present in \code{geno} but absent from \code{map} are dropped
#'   with a message. The map is sorted by position within each chromosome.
#' @param id Character string naming the line identifier column in \code{geno}
#'   when \code{geno} is a \code{data.frame}. Ignored when \code{geno} is a
#'   matrix (row names are used directly). Default is \code{"id"}.
#' @param map.id Character string naming the marker name column in \code{map}.
#'   Default is \code{"marker"}.
#' @param map.chr Character string naming the chromosome identifier column in
#'   \code{map}. Default is \code{"chr"}. Chromosome names containing dots are
#'   replaced with underscores with a warning, as dots are used internally as
#'   field separators in marker naming conventions.
#' @param map.pos Character string naming the cM position column in \code{map}.
#'   Default is \code{"pos"}.
#' @param encoding Character string specifying the genotype encoding in
#'   \code{geno}. Use \code{"012"} (default) when genotypes are coded as
#'   allele counts (0 = AA homozygote, 1 = AB heterozygote, 2 = BB
#'   homozygote); these are recoded to additive \eqn{\pm 1} coding
#'   (\eqn{-1, 0, +1}) internally. Use \code{"pm1"} when genotypes are
#'   already in additive \eqn{\pm 1} coding.
#' @param maf Numeric scalar giving the minor allele frequency threshold below
#'   which markers are removed. Markers with minor allele frequency
#'   \eqn{< maf} are dropped prior to analysis. Default is \code{0.05}.
#'   Set to \code{0} or \code{NULL} to disable filtering.
#'
#' @details
#' The function performs the following steps in order:
#' \enumerate{
#'   \item Extracts line IDs and the genotype matrix from \code{geno}.
#'   \item Aligns markers between \code{geno} and \code{map}, retaining only
#'     the intersection.
#'   \item Re-encodes genotypes to additive \eqn{\pm 1} coding if
#'     \code{encoding = "012"}.
#'   \item Removes markers with minor allele frequency below \code{maf}.
#'   \item Imputes any remaining missing values (\code{NA}) with the
#'     column mean (i.e. the mean genotype at each marker across lines).
#'   \item Splits the marker matrix by chromosome and sorts markers by
#'     position within each chromosome.
#' }
#'
#' The resulting object has the same internal structure as an \code{"interval"}
#' object produced by \code{\link[wgaim]{cross2int}}, meaning it is
#' compatible with the shared wgAim engine without modification.
#' The \code{"panel"} class designation distinguishes it from biparental
#' cross objects in the argument validation of \code{\link{GWASAim}} and
#' \code{\link{GPAim}}.
#'
#' @return An object of class \code{c("panel","interval")} — a list with:
#' \describe{
#'   \item{\code{$pheno}}{A \code{data.frame} with one column (the line
#'     identifier) for all genotyped lines.}
#'   \item{\code{$geno}}{A named list with one element per chromosome. Each
#'     element contains:
#'     \describe{
#'       \item{\code{$data}}{Raw encoded genotype matrix (lines x markers).}
#'       \item{\code{$map}}{Named numeric vector of cM positions.}
#'       \item{\code{$imputed.data}}{Mean-imputed \eqn{\pm 1} coded genotype
#'         matrix used by \code{GWASAim} and \code{GPAim}.}
#'     }
#'   }
#' }
#'
#' @seealso \code{\link{GWASAim}}, \code{\link{GPAim}},
#'   \code{\link[wgaim]{cross2int}}
#'
#' @examples
#' \dontrun{
#' # A simple panel from a 0/1/2 encoded genotype matrix
#' n.lines   <- 200
#' n.markers <- 500
#' geno.mat  <- matrix(sample(0:2, n.lines * n.markers, replace = TRUE),
#'                     nrow = n.lines, ncol = n.markers)
#' rownames(geno.mat) <- paste0("Line_", 1:n.lines)
#' colnames(geno.mat) <- paste0("M",     1:n.markers)
#'
#' map.df <- data.frame(
#'   marker = paste0("M", 1:n.markers),
#'   chr    = rep(paste0("Chr", 1:5), each = 100),
#'   pos    = rep(seq(1, 100), times = 5)
#' )
#'
#' panel <- makePanel(geno = geno.mat, map = map.df,
#'                    encoding = "012", maf = 0.05)
#' class(panel)  # "panel" "interval"
#' names(panel$geno)
#' }
#'
#' @export
makePanel <- function(geno, map, id = "id",
                      map.id  = "marker",
                      map.chr = "chr",
                      map.pos = "pos",
                      encoding = "012",
                      maf = 0.05) {

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
        stop("No common markers found between geno columns and map '", map.id, "' column.")
    n.drop <- length(geno.markers) - length(common)
    if (n.drop > 0)
        message(n.drop, " marker(s) in geno not found in map and will be dropped.")
    geno.mat <- geno.mat[, common, drop = FALSE]
    map      <- map[map[[map.id]] %in% common, , drop = FALSE]

    # -------------------------------------------------------------------------
    # 4. Encode genotypes to additive ±1 coding
    # -------------------------------------------------------------------------
    if (encoding == "012") {
        # 0 (AA) -> -1,  1 (AB) -> 0,  2 (BB) -> +1
        geno.mat <- geno.mat - 1
    } else if (encoding == "pm1") {
        # Already ±1 additive coding — nothing to do
    } else {
        stop("encoding must be '012' (0/1/2 counts) or 'pm1' (already \u00b11).")
    }

    # -------------------------------------------------------------------------
    # 5. MAF filter
    # -------------------------------------------------------------------------
    if (!is.null(maf) && maf > 0) {
        allele.freq <- colMeans(geno.mat, na.rm = TRUE)   # in ±1 coding, mean = 2p - 1
        p           <- (allele.freq + 1) / 2              # convert back to [0,1] freq
        minor.freq  <- pmin(p, 1 - p)
        keep        <- minor.freq >= maf
        n.filt      <- sum(!keep)
        if (n.filt > 0)
            message(n.filt, " marker(s) removed with MAF < ", maf, ".")
        geno.mat <- geno.mat[, keep, drop = FALSE]
        map      <- map[map[[map.id]] %in% colnames(geno.mat), , drop = FALSE]
    }

    # -------------------------------------------------------------------------
    # 6. Mean-impute missing genotypes
    # -------------------------------------------------------------------------
    n.missing <- sum(is.na(geno.mat))
    if (n.missing > 0) {
        message("Imputing ", n.missing, " missing genotype(s) with column means.")
        col.means <- colMeans(geno.mat, na.rm = TRUE)
        for (j in seq_len(ncol(geno.mat))) {
            mis <- is.na(geno.mat[, j])
            if (any(mis))
                geno.mat[mis, j] <- col.means[j]
        }
    }
    rownames(geno.mat) <- ids

    # -------------------------------------------------------------------------
    # 7. Check and sanitise chromosome names (dots break QTL naming)
    # -------------------------------------------------------------------------
    chr.vals <- as.character(map[[map.chr]])
    if (any(grepl("\\.", chr.vals))) {
        warning("Dots found in chromosome names: replacing with underscores ",
                "to avoid conflicts with internal naming conventions.")
        chr.vals        <- gsub("\\.", "_", chr.vals)
        map[[map.chr]]  <- chr.vals
    }

    # -------------------------------------------------------------------------
    # 8. Split by chromosome and build per-chromosome structures
    # -------------------------------------------------------------------------
    chrs     <- unique(chr.vals)
    geno.list <- lapply(chrs, function(ch) {
        chr.map <- map[map[[map.chr]] == ch, , drop = FALSE]
        chr.map <- chr.map[order(as.numeric(chr.map[[map.pos]])), ]
        mnames  <- as.character(chr.map[[map.id]])
        chr.g   <- geno.mat[, mnames, drop = FALSE]
        chr.pos <- setNames(as.numeric(chr.map[[map.pos]]), mnames)
        list(
            data         = chr.g,
            map          = chr.pos,
            imputed.data = chr.g    # already imputed above; marker mode uses this
        )
    })
    names(geno.list) <- chrs

    # -------------------------------------------------------------------------
    # 9. Assemble and return panel object
    # -------------------------------------------------------------------------
    panel <- list(
        pheno = setNames(data.frame(ids, stringsAsFactors = FALSE), id),
        geno  = geno.list
    )
    class(panel) <- c("panel", "interval")

    n.markers.final <- sum(sapply(geno.list, function(el) ncol(el$imputed.data)))
    message("Panel object created: ", length(ids), " lines, ",
            n.markers.final, " markers across ", length(chrs), " chromosomes.")
    panel
}
