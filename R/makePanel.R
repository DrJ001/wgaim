# =============================================================================
# makePanel.R
# Convert association panel marker data into a 'panel'/'interval' class object
# for use with GWASAim().
#
# The resulting object has the same internal structure as an 'interval' object
# from cross2int(), so it feeds directly into the shared wgAim engine pieces
# without modification.
# =============================================================================

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
