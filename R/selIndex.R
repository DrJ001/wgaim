# =============================================================================
# selIndex.R
# Selection index for gpAim objects.
#
# Combines per-environment (or per-trait) GEBVs from a multivariate gpAim fit
# into a single selection index, enabling breeding decisions that simultaneously
# account for performance across multiple environments and economic priorities.
#
# Three index types are supported:
#
#   "weighted"      -- Direct economic index: I = b'g where b = user-supplied
#                      weights.  The simplest and most transparent choice when
#                      the breeder can directly specify the relative value of
#                      performance in each environment.
#
#   "smith-hazel"   -- Smith-Hazel optimal index (Smith 1936; Hazel 1943):
#                      b = P^{-1} Ga w
#                      where P = empirical variance-covariance matrix of GEBVs
#                      across lines (estimated from data), Ga = genetic
#                      covariance matrix from the fitted model, and w =
#                      user-supplied economic weights.  Maximises the
#                      correlation between the index I and the aggregate
#                      genotype H = w'g.  Unlike the "weighted" type, this
#                      accounts for genetic correlations when deriving the
#                      optimal combining weights.
#
#   "desired-gains" -- Pesek-Baker desired gains index (Pesek & Baker 1969):
#                      b proportional to Ga^{-1} d
#                      where d = user-supplied vector of desired genetic gains
#                      per environment.  The breeder specifies HOW MUCH gain
#                      is wanted in each environment rather than the economic
#                      value; the index weights are computed to achieve those
#                      proportional gains.
#
# For univariate gpAim objects the index reduces to the GEBV ranking
# regardless of type; a message is issued and the weighted index is returned.
#
# References:
#   Smith H.F. (1936). A discriminant function for plant selection.
#     Ann. Eugenics 7, 240-250.
#   Hazel L.N. (1943). The genetic basis for constructing selection indexes.
#     Genetics 28, 476-490.
#   Pesek J. & Baker R.J. (1969). Desired improvement in relation to selection
#     indices. Can. J. Plant Sci. 49, 803-804.
# =============================================================================

#' Selection Index for Genomic Prediction Objects
#'
#' @description
#' Combines per-environment GEBVs from a multivariate \code{\link{gpAim}} fit
#' into a single selection index, allowing simultaneous optimisation across
#' multiple environments. Three index types are supported: a direct
#' \strong{weighted} index, the \strong{Smith-Hazel} optimal index, and the
#' \strong{Pesek-Baker} desired-gains index.
#'
#' @param object A fitted \code{"gpAim"} object, as returned by
#'   \code{\link{gpAim}}.  For univariate objects the index reduces to the
#'   GEBV ranking regardless of \code{type}.
#' @param weights Named numeric vector of economic weights, one per
#'   environment/trait level.  Names must match the levels of the \code{Trait}
#'   factor used in \code{gpAim} (i.e. \code{object$GP$trait.levels}).
#'   Required for \code{type = "weighted"} and \code{type = "smith-hazel"}.
#'   If \code{NULL} and \code{type = "weighted"}, equal weights are used.
#' @param desired Named numeric vector of desired genetic gains, one per
#'   environment/trait level.  Names must match \code{object$GP$trait.levels}.
#'   Required for \code{type = "desired-gains"}.
#' @param type Character string specifying the index type.  One of:
#'   \describe{
#'     \item{\code{"weighted"} (default)}{Direct economic index
#'       \eqn{I = \mathbf{b}'\mathbf{g}}, where \eqn{\mathbf{b}} equals the
#'       supplied \code{weights} (or equal weights if \code{NULL}).}
#'     \item{\code{"smith-hazel"}}{Smith-Hazel optimal index
#'       \eqn{\mathbf{b} = \mathbf{P}^{-1}\mathbf{G}_a\mathbf{w}}, where
#'       \eqn{\mathbf{P}} is the empirical variance-covariance matrix of GEBVs
#'       across lines and \eqn{\mathbf{G}_a} is the genetic covariance matrix
#'       from the fitted model.  Maximises the correlation between the index
#'       and the aggregate genotype \eqn{H = \mathbf{w}'\mathbf{g}}.}
#'     \item{\code{"desired-gains"}}{Pesek-Baker desired-gains index
#'       \eqn{\mathbf{b} \propto \mathbf{G}_a^{-1}\mathbf{d}}, where
#'       \eqn{\mathbf{d}} is the user-supplied vector of desired gains.}
#'   }
#' @param standardise Logical.  If \code{TRUE} (default), GEBVs are
#'   standardised to unit variance within each environment before combining,
#'   so that weights are applied on a common scale irrespective of
#'   environment-specific variance differences.
#' @param prop.select Numeric proportion of lines to select.  Used to compute
#'   the selection threshold when \code{threshold = NULL} and \code{selected}
#'   is not supplied.  Default \code{0.10} (top 10\%).
#' @param threshold Numeric.  Explicit selection threshold applied to the index.
#'   If between 0 and 1 exclusive, treated as a quantile of the index
#'   distribution.  If outside (0, 1), treated as a raw index value.  If
#'   \code{NULL} (default), \code{prop.select} is used.  Ignored when
#'   \code{selected} is supplied.
#' @param selected Optional character vector of line identifiers chosen by an
#'   external method (e.g. pedigree-based selection, phenotypic ranking, expert
#'   judgement).  When supplied, \code{prop.select} and \code{threshold} are
#'   ignored for gain computation: \code{$gain} is computed from exactly these
#'   lines, and the selection index is still fitted with equal weights (or
#'   \code{weights}/\code{desired} if supplied) purely for context and display.
#'   All line IDs must be present in the \code{gpAim} object.  Default
#'   \code{NULL} (index-based selection).
#' @param \dots Currently unused.
#'
#' @details
#' \strong{Smith-Hazel index:}
#' The empirical variance-covariance matrix \eqn{\mathbf{P}} is estimated
#' directly from the \eqn{n \times t} matrix of GEBVs (lines \eqn{\times}
#' environments) after any standardisation.  When \eqn{\mathbf{P}} is
#' numerically singular (e.g. when the number of environments approaches
#' the number of lines) a small ridge is added to the diagonal.
#'
#' \strong{Desired-gains index:}
#' \eqn{\mathbf{b}} is scaled so that the expected gain in each environment
#' is proportional to the corresponding entry of \code{desired}.  The
#' magnitude of \eqn{\mathbf{b}} is arbitrary (only the direction matters
#' for ranking), so the weights are re-scaled to have unit L2 norm before
#' indexing.
#'
#' \strong{Univariate objects:}
#' When \code{object$GP$Trait} is \code{NULL}, all three index types reduce
#' to ranking by GEBV.  A message is issued and the \code{"weighted"} index
#' is returned.
#'
#' @return An object of class \code{c("selIndex", "list")} containing:
#' \describe{
#'   \item{\code{$index}}{A \code{data.frame} with one row per line, sorted
#'     descending by \code{Index}.  Columns: line identifier, one GEBV column
#'     per environment (MV only), \code{Index}, \code{Rank}.}
#'   \item{\code{$weights}}{The index weight vector \eqn{\mathbf{b}} actually
#'     used (named, one per environment for MV; length-one scalar for univariate).}
#'   \item{\code{$type}}{The index type string.}
#'   \item{\code{$standardise}}{Logical: whether GEBVs were standardised.}
#'   \item{\code{$trait.levels}}{Environment/trait level names (NULL for
#'     univariate).}
#'   \item{\code{$genetic.term}}{Line identifier column name.}
#'   \item{\code{$Ga}}{Genetic covariance matrix (NULL for univariate).}
#'   \item{\code{$Gcor}}{Genetic correlation matrix (NULL for univariate).}
#'   \item{\code{$gain}}{Named numeric vector (one entry per environment for
#'     multivariate; scalar for univariate) of expected genetic gains at the
#'     chosen selection threshold:
#'     \eqn{\Delta G_j = \bar{g}_j^{\text{selected}} - \bar{g}_j^{\text{all}}}.
#'     Gains are on the original GEBV scale (not standardised).}
#'   \item{\code{$selected}}{Character vector of the selected line IDs when
#'     \code{selected} was supplied; \code{NULL} otherwise.}
#'   \item{\code{$thr}}{The index threshold value used to define the selected
#'     group; \code{NA} when \code{selected} was supplied externally.}
#'   \item{\code{$n.selected}}{Number of lines in the selected group.}
#'   \item{\code{$prop.select}}{The \code{prop.select} value passed at
#'     construction time.}
#'   \item{\code{$n.lines}}{Number of lines.}
#'   \item{\code{$n.environments}}{Number of environments (1 for univariate).}
#' }
#'
#' @seealso \code{\link{print.selIndex}}, \code{\link{summary.selIndex}},
#'   \code{\link{plot.selIndex}}, \code{\link{gpAim}}
#'
#' @references
#' Smith H.F. (1936). A discriminant function for plant selection.
#' \emph{Ann. Eugenics} \strong{7}, 240--250.
#'
#' Hazel L.N. (1943). The genetic basis for constructing selection indexes.
#' \emph{Genetics} \strong{28}, 476--490.
#'
#' Pesek J. & Baker R.J. (1969). Desired improvement in relation to
#' selection indices. \emph{Can. J. Plant Sci.} \strong{49}, 803--804.
#'
#' @examples
#' \dontrun{
#' # After running multivariate gpAim with Trait = "Trial":
#' # Equal weights across all trials (default)
#' si <- selIndex(gp.mv)
#'
#' # Economic weights: Trial 1 worth twice Trial 2
#' si <- selIndex(gp.mv, weights = c(Trial1 = 2, Trial2 = 1))
#'
#' # Smith-Hazel optimal index
#' si <- selIndex(gp.mv, weights = c(Trial1 = 1, Trial2 = 1),
#'                type = "smith-hazel")
#'
#' # Desired-gains index: want 0.5 units gain in Trial1, 0.3 in Trial2
#' si <- selIndex(gp.mv, desired = c(Trial1 = 0.5, Trial2 = 0.3),
#'                type = "desired-gains")
#'
#' print(si)
#' summary(si)
#' plot(si)
#' plot(si, type = "biplot")
#' plot(si, type = "weights")
#' }
#'
#' @name selIndex
#' @export
selIndex <- function(object, ...)
    UseMethod("selIndex")

#' @rdname selIndex
#' @exportS3Method
selIndex.default <- function(object, ...)
    stop("selIndex() currently only supports objects of class \"gpAim\".")

#' @rdname selIndex
#' @exportS3Method
selIndex.gpAim <- function(object,
                            weights     = NULL,
                            desired     = NULL,
                            type        = c("weighted", "smith-hazel", "desired-gains"),
                            standardise = TRUE,
                            prop.select = 0.10,
                            threshold   = NULL,
                            selected    = NULL,
                            ...) {

    type <- match.arg(type)
    gp   <- object$GP
    is.mv <- !is.null(gp$Trait)

    # -------------------------------------------------------------------------
    # Univariate path -- index reduces to GEBV ranking
    # -------------------------------------------------------------------------
    if (!is.mv) {
        message("selIndex: univariate gpAim -- index is equivalent to GEBV ranking.")
        gebv <- gp$gebv
        gebv <- gebv[order(gebv$GEBV, decreasing = TRUE), ]
        rownames(gebv) <- NULL
        idx.df <- data.frame(
            Line  = gebv[[gp$genetic.term]],
            GEBV  = round(gebv$GEBV, 4),
            Index = round(gebv$GEBV, 4),
            Rank  = seq_len(nrow(gebv)),
            stringsAsFactors = FALSE
        )
        names(idx.df)[1] <- gp$genetic.term
        # Selection: externally supplied lines OR threshold-based
        all.ids <- gebv[[gp$genetic.term]]
        if (!is.null(selected)) {
            selected <- as.character(selected)
            miss <- setdiff(selected, as.character(all.ids))
            if (length(miss))
                stop("selected: line(s) not found in genotypic data: ",
                     paste(miss, collapse = ", "))
            sel.uni    <- as.character(all.ids) %in% selected
            thr.uni    <- NA_real_
            ps.uni     <- length(selected) / nrow(gebv)
        } else {
            thr.uni    <- .si_threshold(gebv$GEBV, threshold, prop.select)
            sel.uni    <- gebv$GEBV >= thr.uni
            ps.uni     <- prop.select
            selected   <- NULL
        }
        gain.uni <- mean(gebv$GEBV[sel.uni]) - mean(gebv$GEBV)
        result <- list(
            index          = idx.df,
            weights        = setNames(1, gp$genetic.term),
            type           = "weighted",
            standardise    = FALSE,
            prop.select    = ps.uni,
            thr            = thr.uni,
            n.selected     = sum(sel.uni),
            gain           = gain.uni,
            selected       = selected,
            trait.levels   = NULL,
            genetic.term   = gp$genetic.term,
            Ga             = NULL,
            Gcor           = NULL,
            n.lines        = nrow(idx.df),
            n.environments = 1L
        )
        class(result) <- c("selIndex", "list")
        return(invisible(result))
    }

    # -------------------------------------------------------------------------
    # Multivariate path
    # -------------------------------------------------------------------------
    tl   <- gp$trait.levels
    nt   <- length(tl)
    Ga   <- gp$Ga
    Gcor <- gp$Gcor
    gterm <- gp$genetic.term

    # Build the (lines x environments) GEBV matrix in trait-level order
    gebv.long <- gp$gebv
    line.ids  <- unique(gebv.long[[gterm]])
    n.lines   <- length(line.ids)

    G.mat <- matrix(NA_real_, nrow = n.lines, ncol = nt,
                    dimnames = list(line.ids, tl))
    for (j in seq_len(nt)) {
        sub.j  <- gebv.long[gebv.long[[gp$Trait]] == tl[j], ]
        # Match by line ID to handle any row ordering in gebv
        G.mat[as.character(sub.j[[gterm]]), j] <- sub.j$GEBV
    }

    # -------------------------------------------------------------------------
    # Validate weights / desired arguments
    # -------------------------------------------------------------------------
    if (type %in% c("weighted", "smith-hazel")) {
        if (is.null(weights)) {
            message("selIndex: weights not supplied -- using equal weights.")
            weights <- setNames(rep(1 / nt, nt), tl)
        } else {
            weights <- .check_index_vector(weights, tl, "weights")
        }
    }
    if (type == "desired-gains") {
        if (is.null(desired))
            stop("type = \"desired-gains\" requires a named 'desired' vector ",
                 "with one entry per environment.")
        desired <- .check_index_vector(desired, tl, "desired")
    }

    # -------------------------------------------------------------------------
    # Standardise GEBVs per environment (unit variance) if requested
    # -------------------------------------------------------------------------
    sd.vec <- apply(G.mat, 2, sd, na.rm = TRUE)
    sd.vec[sd.vec < .Machine$double.eps^0.5] <- 1  # guard against zero-variance

    G.std <- if (standardise) {
        sweep(G.mat, 2, sd.vec, "/")
    } else {
        G.mat
    }

    # -------------------------------------------------------------------------
    # Compute index weights b
    # -------------------------------------------------------------------------
    b <- switch(type,

        "weighted" = {
            setNames(as.numeric(weights), tl)
        },

        "smith-hazel" = {
            # Empirical variance-covariance matrix of GEBVs across lines
            P <- stats::cov(G.std)
            # Add a small ridge if P is near-singular
            ev.min <- min(eigen(P, only.values = TRUE)$values)
            if (ev.min < 1e-8) {
                ridge <- max(1e-6, abs(ev.min) + 1e-8)
                P     <- P + diag(ridge, nt)
            }
            # Ga on same scale as G.std (divide rows & cols by sd.vec)
            Ga.std <- sweep(sweep(Ga, 1, sd.vec, "/"), 2, sd.vec, "/")
            # b = P^{-1} Ga w  (w on standardised scale too)
            w.std <- weights / sd.vec
            b.raw <- solve(P, Ga.std %*% w.std)
            setNames(as.numeric(b.raw), tl)
        },

        "desired-gains" = {
            # Ga on original scale; b proportional to Ga^{-1} d
            ev.min <- min(eigen(Ga, only.values = TRUE)$values)
            Ga.reg <- if (ev.min < 1e-8) {
                Ga + diag(max(1e-6, abs(ev.min) + 1e-8), nt)
            } else Ga
            b.raw <- solve(Ga.reg, desired)
            # Rescale to unit L2 norm (direction only matters for ranking)
            b.raw <- b.raw / sqrt(sum(b.raw^2))
            setNames(as.numeric(b.raw), tl)
        }
    )

    # -------------------------------------------------------------------------
    # Compute index scores
    # -------------------------------------------------------------------------
    index.vals <- as.numeric(G.std %*% b)
    names(index.vals) <- line.ids

    # -------------------------------------------------------------------------
    # Assemble output data frame: line ID + per-env GEBVs + index + rank
    # -------------------------------------------------------------------------
    idx.df <- as.data.frame(G.mat, stringsAsFactors = FALSE)
    idx.df <- round(idx.df, 4)
    idx.df <- cbind(
        data.frame(Line = line.ids, stringsAsFactors = FALSE),
        idx.df,
        data.frame(
            Index = round(index.vals, 4),
            stringsAsFactors = FALSE
        )
    )
    names(idx.df)[1] <- gterm
    idx.df <- idx.df[order(idx.df$Index, decreasing = TRUE), ]
    idx.df$Rank <- seq_len(nrow(idx.df))
    rownames(idx.df) <- NULL

    # -------------------------------------------------------------------------
    # Selection: externally supplied lines OR threshold from index
    # ΔG_j = mean(GEBV_j[selected]) − mean(GEBV_j[all])  on original GEBV scale
    # -------------------------------------------------------------------------
    if (!is.null(selected)) {
        selected  <- as.character(selected)
        miss      <- setdiff(selected, line.ids)
        if (length(miss))
            stop("selected: line(s) not found in genotypic data: ",
                 paste(miss, collapse = ", "))
        sel.lines  <- selected
        thr        <- NA_real_
        ps.stored  <- length(selected) / n.lines
    } else {
        thr        <- .si_threshold(index.vals, threshold, prop.select)
        sel.lines  <- line.ids[index.vals >= thr]
        ps.stored  <- prop.select
        selected   <- NULL
    }

    gain <- setNames(
        vapply(tl, function(j)
            mean(G.mat[sel.lines, j], na.rm = TRUE) -
            mean(G.mat[,           j], na.rm = TRUE),
            numeric(1L)),
        tl)

    # -------------------------------------------------------------------------
    # Return
    # -------------------------------------------------------------------------
    result <- list(
        index          = idx.df,
        weights        = b,
        type           = type,
        standardise    = standardise,
        prop.select    = ps.stored,
        thr            = thr,
        n.selected     = length(sel.lines),
        gain           = gain,
        selected       = selected,
        trait.levels   = tl,
        genetic.term   = gterm,
        Ga             = Ga,
        Gcor           = Gcor,
        n.lines        = n.lines,
        n.environments = nt
    )
    class(result) <- c("selIndex", "list")
    invisible(result)
}

# =============================================================================
# Internal helper: resolve a selection threshold from a raw value or a
# proportion of lines to select.  Mirrors .gp_threshold() in plot.gpAim.R
# but lives here so selIndex.R has no dependency on a plot helper.
# =============================================================================
#' @keywords internal
.si_threshold <- function(vals, threshold, prop.select) {
    if (!is.null(threshold)) {
        if (threshold > 0 && threshold < 1)
            stats::quantile(vals, threshold, names = FALSE)
        else
            as.numeric(threshold)
    } else {
        stats::quantile(vals, 1 - prop.select, names = FALSE)
    }
}

# =============================================================================
# Internal helper: validate and reorder a user-supplied named weight/desired
# vector to match trait level ordering, with informative errors.
# =============================================================================
#' @keywords internal
.check_index_vector <- function(v, tl, arg.name) {
    if (is.null(names(v))) {
        if (length(v) != length(tl))
            stop(arg.name, ": must be a named vector with one entry per environment ",
                 "(", paste(tl, collapse = ", "), "), or an unnamed vector of ",
                 "length ", length(tl), ".")
        names(v) <- tl
        return(as.numeric(v))
    }
    miss <- setdiff(tl, names(v))
    if (length(miss))
        stop(arg.name, ": missing entries for environment(s): ",
             paste(miss, collapse = ", "), ".\n",
             "  Required: ", paste(tl, collapse = ", "))
    extra <- setdiff(names(v), tl)
    if (length(extra))
        warning(arg.name, ": ignoring extra entries not in trait levels: ",
                paste(extra, collapse = ", "))
    setNames(as.numeric(v[tl]), tl)
}
