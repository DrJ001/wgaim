# =============================================================================
# fineMap.R
# Fine-mapping around selected QTL / markers from qtlAim or gwasAim models.
#
# For qtlAim models a fresh inferred-marker grid (cM step = 'step') is
# computed inside the window around each nominated QTL, regardless of the
# gen.type or infer setting used in primeCross().  This always uses the
# lambdaf() weight-matrix approach.
#
# For gwasAim models the window is scanned using the actual markers in the
# panel (no midpoint imputation).
#
# In both cases the scan loop refits the ASReml model with each window
# candidate substituted for the existing fixed QTL / marker term, and
# reports z-ratio derived p-values and LOD scores.
# =============================================================================

#' Fine-map QTL or GWAS markers
#'
#' @description
#' Performs single-QTL fine mapping in a window around each significant QTL
#' (for \code{\link{qtlAim}} models) or marker (for \code{\link{gwasAim}}
#' models) found in a fitted model.  Both univariate and multivariate
#' (\code{Trait}) models are supported.
#'
#' For \code{qtlAim} models a dense inferred-marker grid is constructed inside
#' the window using the Haldane-map weight approach of \code{\link{primeCross}},
#' irrespective of the \code{type} or \code{infer} argument originally passed
#' to \code{primeCross()}.  Each grid position is tested as a fixed effect in
#' turn, with all other detected QTL effects held in the model.
#'
#' For \code{gwasAim} models the window is defined around the significant
#' marker (in cM) and all panel markers within the window are tested
#' individually; no inferred midpoints are added.
#'
#' For multivariate models, QTL that were retained as main effects after
#' Wald pruning are scanned using the same single-df z-ratio as the
#' univariate case.  QTL retained as \code{Trait}-interaction effects are
#' scanned by substituting \code{Trait:candidate} for \code{Trait:QTL} and
#' computing a joint Wald statistic across all non-aliased interaction
#' coefficients via an internal Wald test; the LOD score is derived as
#' \eqn{\chi^2 / (2 \log 10)}, consistent with the univariate formula.
#'
#' @param object A fitted object of class \code{"qtlAim"} or \code{"gwasAim"}.
#' @param genObj The \code{wgCross} or \code{wgPanel} object used in the
#'   analysis, produced by \code{\link{primeCross}} or \code{\link{primePanel}}.
#' @param qtl Character vector of QTL / marker keys to fine-map, using the
#'   internal \code{"Chr.<chr>.<idx>"} format stored in \code{object$QTL$qtl}.
#'   If \code{NULL} (default) all detected QTL / markers are fine-mapped.
#' @param window Numeric scalar.  Half-width of the fine-mapping window in cM.
#'   Default \code{50}.
#' @param step Numeric scalar.  cM step size for the inferred-marker grid.
#'   Only used for \code{qtlAim} models.  Smaller values give finer resolution
#'   at the cost of more model refits.  Default \code{2}.
#' @param exclusion.window Numeric scalar.  cM exclusion zone around other
#'   (non-target) detected QTL during each scan.  Default \code{20}.
#' @param \dots Additional arguments passed to \code{update.asreml}, e.g.\
#'   \code{na.action = na.method(x = "include")}.
#'
#' @return An object of class \code{"fineMap"}: a named list with one element
#'   per fine-mapped QTL.  Each element is a data frame with columns:
#'   \describe{
#'     \item{\code{mark}}{Name of the inferred grid point or panel marker.}
#'     \item{\code{dist}}{cM position on the chromosome.}
#'     \item{\code{pvalue}}{Two-tailed p-value from the z-ratio test.}
#'     \item{\code{LOD}}{LOD score derived from the z-ratio.}
#'   }
#'   The list has attribute \code{"qtl.key"} recording the original QTL key
#'   names, and attribute \code{"type"} recording \code{"qtlAim"} or
#'   \code{"gwasAim"}.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{gwasAim}},
#'   \code{\link{primeCross}}, \code{\link{primePanel}},
#'   \code{\link{print.fineMap}}, \code{\link{plot.fineMap}}
#'
#' @export
fineMap <- function(object, genObj, qtl = NULL, window = 50, step = 2,
                    exclusion.window = 20, ...)
{
    # -------------------------------------------------------------------------
    # 0.  Class checks and setup
    # -------------------------------------------------------------------------
    if (!inherits(object, c("qtlAim", "gwasAim")))
        stop("object must be of class \"qtlAim\" or \"gwasAim\".")

    is_gwas <- inherits(object, "gwasAim")

    if (is_gwas) {
        if (!inherits(genObj, "wgPanel"))
            stop("genObj must be of class \"wgPanel\" produced by primePanel().")
    } else {
        if (!inherits(genObj, "wgCross"))
            stop("genObj must be of class \"wgCross\" produced by primeCross().")
    }

    if (is.null(object$QTL$qtl) || length(object$QTL$qtl) == 0L)
        stop("No significant QTL / markers found in object.")

    # Resolve which QTL to fine-map
    all_qtl <- object$QTL$qtl          # e.g. "Chr.1A.3", "Chr.2B.11"
    if (!is.null(qtl)) {
        bad <- setdiff(qtl, all_qtl)
        if (length(bad))
            stop("The following QTL keys were not found in object$QTL$qtl: ",
                 paste(bad, collapse = ", "))
        targets <- qtl
    } else {
        targets <- all_qtl
    }

    # -------------------------------------------------------------------------
    # 1.  Retrieve stored information from the fitted model
    # -------------------------------------------------------------------------
    gterm    <- object$QTL$diag$genetic.term
    state    <- object$QTL$diag$state
    method   <- object$QTL$method     # always "fixed" for gwasAim

    # Multivariate flag and per-QTL interaction flag
    is_mv    <- !is.null(object$QTL$Trait)
    mv_trait <- object$QTL$Trait      # e.g. "Trial" or NULL

    # Retrieve phenoData: qtlAim/gwasAim write it back to the caller's
    # environment under the name  "<response>.data" (e.g. "yld.data").
    # The formula environment still holds a reference to that environment.
    resp      <- deparse(object$call$fixed[[2]])
    data_name <- paste0(resp, ".data")
    fe        <- environment(object$call$fixed)
    phenoData <- tryCatch(
        get(data_name, envir = fe, inherits = TRUE),
        error = function(e)
            stop("Cannot locate phenoData ('", data_name, "') -- ",
                 "ensure fineMap() is called in the same session as ",
                 "qtlAim() / gwasAim().")
    )

    # -------------------------------------------------------------------------
    # 2.  Build the background genoData matrix used for covObj rebuilds.
    #     This must match exactly what the engine used:
    #       - qtlAim interval type  -> interval.data  (ncol = n intervals)
    #       - qtlAim marker type    -> imputed.data   (ncol = n markers)
    #       - gwasAim               -> imputed.data   (ncol = n markers)
    # -------------------------------------------------------------------------
    bg_type <- object$QTL$type   # "interval" or "marker"
    if (!is_gwas && bg_type == "interval") {
        gdat <- lapply(genObj$geno, function(el) el$interval.data)
    } else {
        gdat <- lapply(genObj$geno, function(el) el$imputed.data)
    }
    genoData <- do.call("cbind", gdat)
    # state names are the canonical column names
    colnames(genoData) <- names(state)
    rownames(genoData) <- as.character(genObj$pheno[[gterm]])
    genoData <- genoData[
        rownames(genoData) %in% as.character(phenoData[[gterm]]), ,
        drop = FALSE
    ]

    # Chromosomes index for each state column
    chrs <- sapply(strsplit(names(state), "\\."), "[", 2)

    # -------------------------------------------------------------------------
    # 3.  Fine-map each target QTL
    # -------------------------------------------------------------------------
    results <- vector("list", length(targets))
    names(results) <- targets

    for (tgt in targets) {

        tgt_parts <- strsplit(tgt, "\\.")[[1]]  # ["Chr", chrname, idx]
        mchr      <- tgt_parts[2]
        tgt_idx   <- as.integer(tgt_parts[3])
        tgt_col   <- paste0("X.", mchr, ".", tgt_idx)  # column name in phenoData

        # The existing QTL fixed effect column name in phenoData
        mark_qtl  <- gsub("Chr\\.", "X.", tgt)

        # For MV: did the Wald test keep the interaction for this QTL?
        # object$QTL$is.interaction is a named logical vector (name = mark_qtl or
        # its interaction form); we match by position in object$QTL$qtl.
        tgt_is_interaction <- if (is_mv && !is.null(object$QTL$is.interaction)) {
            tgt_pos <- match(tgt, all_qtl)
            if (!is.na(tgt_pos) && tgt_pos <= length(object$QTL$is.interaction))
                object$QTL$is.interaction[tgt_pos]
            else
                TRUE   # conservative default: assume interaction present
        } else {
            FALSE
        }

        # ----------------------------------------------------------------
        # 3a.  Derive the fine-mapping map and genotype matrix for window
        # ----------------------------------------------------------------
        chr_map  <- genObj$geno[[mchr]]$map   # named cM positions of markers

        if (is_gwas) {
            # GWAS: scan actual markers within +/-window cM of the significant
            # marker's cM position
            sig_pos  <- chr_map[tgt_idx]
            ql       <- max(1L, which(chr_map >= sig_pos - window)[1L])
            qr_cands <- which(chr_map <= sig_pos + window)
            qr       <- if (length(qr_cands)) max(qr_cands) else length(chr_map)

            # Warn and record when the chromosome is shorter than the window
            clipped <- chr_map[ql] > sig_pos - window || chr_map[qr] < sig_pos + window
            if (clipped) {
                message(sprintf(
                    "  Note: %s -- window clipped to chromosome boundaries [%.1f, %.1f] cM",
                    tgt, chr_map[ql], chr_map[qr]
                ))
            }

            fine_map    <- chr_map[ql:qr]
            fine_geno   <- genObj$geno[[mchr]]$imputed.data[
                dimnames(genoData)[[1]], ql:qr, drop = FALSE
            ]
            fine_names  <- names(fine_map)
            grid_type   <- "marker"

        } else {
            # qtlAim: use inferred.map to find the QTL cM position,
            # then build a fresh fine grid via lambdaf()
            if (!is.null(genObj$geno[[mchr]]$inferred.map) &&
                length(genObj$geno[[mchr]]$inferred.map) >= tgt_idx) {
                qtl_pos <- genObj$geno[[mchr]]$inferred.map[tgt_idx]
            } else {
                qtl_pos <- chr_map[tgt_idx]
            }

            win_lo  <- max(0,            qtl_pos - window)
            win_hi  <- min(max(chr_map), qtl_pos + window)

            # Warn and record when the chromosome is shorter than the window
            clipped <- win_lo > qtl_pos - window || win_hi < qtl_pos + window
            if (clipped) {
                message(sprintf(
                    "  Note: %s -- window clipped to chromosome boundaries [%.1f, %.1f] cM",
                    tgt, win_lo, win_hi
                ))
            }

            # Subset markers flanking the window (one extra each side for weights)
            lo_idx  <- max(1L,              findInterval(win_lo, chr_map))
            hi_idx  <- min(length(chr_map), findInterval(win_hi, chr_map) + 1L)
            sub_map <- chr_map[lo_idx:hi_idx]

            # Build the fine grid by stepping outward from qtl_pos so that
            # qtl_pos is guaranteed to be an exact grid point.
            left_seq  <- seq(qtl_pos,        win_lo, by = -step)
            right_seq <- seq(qtl_pos + step,  win_hi, by =  step)
            imark     <- sort(c(rev(left_seq), right_seq))
            imark     <- imark[imark >= win_lo & imark <= win_hi]

            # Compute Haldane interpolation weights for each grid position.
            # For position p between markers at lm and rm (absolute cM):
            #   theta_p  = 0.5*(1 - exp(-2*(p  - lm)/100))
            #   theta_LR = 0.5*(1 - exp(-2*(rm - lm)/100))
            #   w_L = (1 - theta_LR - theta_p)*(theta_LR - theta_p) /
            #              (theta_LR*(1 - theta_LR)*(1 - 2*theta_p))
            #   w_R = theta_p*(1 - theta_p)*(1 - 2*theta_LR) /
            #              (theta_LR*(1 - theta_LR)*(1 - 2*theta_p))
            n_marks <- length(sub_map)
            n_grid  <- length(imark)
            lambda  <- matrix(0, nrow = n_marks, ncol = n_grid)

            for (k in seq_len(n_grid)) {
                p   <- imark[k]
                j   <- findInterval(p, sub_map)
                j   <- max(1L, min(n_marks - 1L, j))
                lm_k  <- sub_map[j]
                rm_k  <- sub_map[j + 1L]
                theta_p  <- 0.5 * (1 - exp(-2 * (p   - lm_k) / 100))
                theta_LR <- 0.5 * (1 - exp(-2 * (rm_k - lm_k) / 100))
                den <- theta_LR * (1 - theta_LR) * (1 - 2 * theta_p)
                if (abs(den) < 1e-12) {
                    # p coincides with a marker -- full weight on that marker
                    lambda[j, k] <- 1
                } else {
                    lambda[j,     k] <- (1 - theta_LR - theta_p) *
                                        (theta_LR - theta_p) / den
                    lambda[j + 1L, k] <- theta_p * (1 - theta_p) *
                                         (1 - 2 * theta_LR) / den
                }
            }

            # Subsetted imputed marker matrix for this chromosome / window
            sub_imp <- genObj$geno[[mchr]]$imputed.data[
                dimnames(genoData)[[1]], lo_idx:hi_idx, drop = FALSE
            ]

            # Fine interval scores:  lines x grid_positions
            fine_geno  <- sub_imp %*% lambda
            fine_map   <- imark
            fine_names <- paste0("fm", seq_along(imark))
            colnames(fine_geno) <- fine_names
            grid_type  <- "interval"
        }

        n_pos <- length(fine_map)
        if (n_pos == 0L) {
            message("No positions in window for QTL ", tgt, " -- skipping.")
            results[[tgt]] <- data.frame(mark = character(), dist = numeric(),
                                         pvalue = numeric(), LOD = numeric())
            next
        }

        # ----------------------------------------------------------------
        # 3b.  Build the exclusion-adjusted state for other QTL on this chr
        # ----------------------------------------------------------------
        other_qtl  <- setdiff(all_qtl, tgt)
        state_work <- state

        if (length(other_qtl)) {
            for (oq in other_qtl) {
                oq_parts <- strsplit(oq, "\\.")[[1]]
                ochr     <- oq_parts[2]
                oidx     <- as.integer(oq_parts[3])
                if (ochr == mchr) {
                    if (!is.null(genObj$geno[[mchr]]$inferred.map) &&
                        length(genObj$geno[[mchr]]$inferred.map) >= oidx) {
                        opos <- genObj$geno[[mchr]]$inferred.map[oidx]
                    } else {
                        opos <- chr_map[oidx]
                    }
                    # Exclude state columns within exclusion.window of this other QTL
                    col_pos <- if (!is_gwas && !is.null(genObj$geno[[mchr]]$inferred.map)) {
                        genObj$geno[[mchr]]$inferred.map
                    } else {
                        chr_map
                    }
                    col_names <- paste0("Chr.", mchr, ".", seq_along(col_pos))
                    wind_idx  <- abs(col_pos - opos) <= exclusion.window
                    state_work[col_names[wind_idx]] <- 0L
                }
            }
        }

        # ----------------------------------------------------------------
        # 3c.  Remove target QTL column from phenoData for scan
        # ----------------------------------------------------------------
        phenoScan <- phenoData
        if (mark_qtl %in% names(phenoScan))
            phenoScan <- phenoScan[, !(names(phenoScan) %in% mark_qtl), drop = FALSE]

        # Add ordering handle so merge doesn't scramble rows
        phenoScan$.ord <- seq_len(nrow(phenoScan))

        # ----------------------------------------------------------------
        # 3d.  Scan loop
        # ----------------------------------------------------------------
        pvalue <- numeric(n_pos)
        lod    <- numeric(n_pos)

        for (k in seq_len(n_pos)) {

            cand_name  <- fine_names[k]
            cand_x     <- paste0("X.", mchr, ".", k, "_fm")  # unique col name
            cand_col   <- fine_geno[, k, drop = FALSE]
            colnames(cand_col) <- cand_x

            # Merge candidate genotype column into scan phenoData
            tmp        <- cbind.data.frame(
                rownames(cand_col), as.vector(cand_col),
                stringsAsFactors = FALSE
            )
            names(tmp) <- c(gterm, cand_x)
            pScan_k    <- merge(phenoScan, tmp, by = gterm, all.x = TRUE)
            pScan_k    <- pScan_k[order(pScan_k$.ord), ]

            # Rebuild covObj excluding the target QTL's state column
            mout    <- which(!as.logical(state_work))
            genoSub <- if (length(mout)) genoData[, -mout, drop = FALSE] else genoData

            if (ncol(genoSub) > nrow(genoSub)) {
                cov_env_k <- .constructCM(genoSub)
                covObj    <- cov_env_k$relm
            } else {
                covObj    <- cbind.data.frame(rownames(genoSub), genoSub)
                names(covObj)[1] <- gterm
            }
            assign("covObj", covObj, envir = parent.frame())

            # ------------------------------------------------------------------
            # Update formula: swap old QTL term for candidate.
            #
            # MV interaction QTL (is.interaction = TRUE):
            #   old term = "Trait:mark_qtl"  →  new term = "Trait:cand_x"
            #
            # UV or MV main-effect QTL (is.interaction = FALSE):
            #   old term = "mark_qtl"         →  new term = "cand_x"
            # ------------------------------------------------------------------
            fixed_labels <- labels(terms(object$call$fixed))
            if (is_mv && isTRUE(tgt_is_interaction)) {
                old_term  <- paste0(mv_trait, ":", mark_qtl)
                new_term  <- paste0(mv_trait, ":", cand_x)
            } else {
                old_term  <- mark_qtl
                new_term  <- cand_x
            }
            remove_str <- if (old_term %in% fixed_labels) paste("-", old_term) else ""

            if (method == "fixed") {
                fix_form  <- as.formula(paste(". ~ . +", new_term, remove_str))
                tempmodel <- update(object, fixed. = fix_form, data = pScan_k, ...)
            } else {
                ran_form  <- update(
                    formula(object$call$random),
                    as.formula(paste("~ . -", mark_qtl, "+", cand_x))
                )
                tempmodel <- .vModify(object, gterm)
                tempmodel <- update(tempmodel, random. = ran_form, data = pScan_k, ...)
            }

            # ------------------------------------------------------------------
            # Extract test statistic for the candidate effect.
            #
            # UV or MV main-effect QTL:
            #   Single-df z-ratio from the candidate's main-effect coefficient.
            #
            # MV interaction QTL:
            #   Joint Wald chi-square via .waldTest() across all non-aliased
            #   Trait:cand_x coefficients.  LOD = stat / (2 * log(10)),
            #   consistent with the UV formula since zrat^2 = 1-df Wald stat.
            # ------------------------------------------------------------------
            cf  <- tempmodel$coefficients[[method]]
            whr <- grep(cand_x, rownames(cf), fixed = TRUE)
            if (length(whr) == 0L) {
                pvalue[k] <- NA_real_
                lod[k]    <- NA_real_
                next
            }

            if (is_mv && isTRUE(tgt_is_interaction)) {
                # Drop aliased (zero-variance) coefficients before Wald test
                vcf_all <- tempmodel$vcoeff[[method]][whr]
                whr_nz  <- whr[vcf_all > 0]
                if (length(whr_nz) == 0L) {
                    pvalue[k] <- NA_real_
                    lod[k]    <- NA_real_
                    next
                }
                wt        <- .waldTest(tempmodel, cc = list(list(coef = whr_nz, type = "zero")))
                stat      <- wt[1L, "Wald Statistic"]
                pvalue[k] <- round(wt[1L, "P-Value"], 4L)
                lod[k]    <- round(stat / (2 * log(10)), 4L)
            } else {
                mcf  <- cf[whr[1L], 1L]
                vcf  <- tempmodel$vcoeff[[method]][whr[1L]]
                zrat <- mcf / sqrt(vcf * tempmodel$sigma2)
                pvalue[k] <- round(1 - pchisq(zrat^2, df = 1L), 4L)
                lod[k]    <- round(0.5 * log(exp(zrat^2), base = 10), 4L)
            }
        }

        df <- data.frame(
            mark   = fine_names,
            dist   = as.numeric(fine_map),
            pvalue = pvalue,
            LOD    = lod,
            stringsAsFactors = FALSE
        )
        # Store original QTL position and clipping flag for plot.fineMap()
        attr(df, "qtl_pos") <- if (is_gwas) as.numeric(sig_pos) else as.numeric(qtl_pos)
        attr(df, "clipped")  <- clipped
        results[[tgt]] <- df
    }

    attr(results, "qtl.key") <- targets
    attr(results, "type")    <- if (is_gwas) "gwasAim" else "qtlAim"
    class(results) <- "fineMap"
    results
}
