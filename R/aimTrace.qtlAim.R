# =============================================================================
# aimTrace.qtlAim.R
# aimTrace() generic + S3 method for qtlAim objects.
# Shared internal helpers (.plot_lrt, .build_stability_df, .plot_stability)
# are also used by aimTrace.gwasAim (defined in aimTrace.gwasAim.R).
# =============================================================================

#' Trace the forward-selection algorithm
#'
#' @description
#' Prints a p-value matrix and likelihood ratio test (LRT) table showing the
#' incremental state of the forward-selection algorithm across iterations,
#' and optionally returns diagnostic plots.
#'
#' The \strong{p-value matrix} shows the p-value of each detected QTL (or
#' marker) at every iteration in which it was present in the model, making it
#' easy to assess whether effects remain stable as further QTL enter.
#'
#' The \strong{LRT table} records the base and full model log-likelihoods, the
#' LRT statistic, and its p-value at every iteration including the final
#' non-significant one.
#'
#' For multivariate analyses (\code{Trait} non-\code{NULL}), the p-value
#' matrix shows the main-effect p-value per QTL per iteration, and the
#' stability plot shows per-trial effect estimates coloured by trial level.
#'
#' @param object A fitted object of class \code{"qtlAim"} or \code{"gwasAim"}.
#' @param iter Integer vector of iterations to include in the p-value matrix.
#'   Default is all iterations: \code{1:length(object$QTL$effects)}.
#' @param lik.out Logical. If \code{TRUE} (default), print the LRT table.
#' @param plot Controls optional diagnostic plot output (console output always
#'   printed regardless). One of:
#'   \describe{
#'     \item{\code{FALSE} (default)}{No plot is returned.}
#'     \item{\code{"lrt"}}{Returns a \code{ggplot} of the LRT statistic across
#'       iterations with the significance threshold marked and each detected
#'       QTL labelled.}
#'     \item{\code{"stability"}}{Returns a \code{ggplot} of QTL effect
#'       estimates \eqn{\pm}1 SE across every iteration in which the QTL was
#'       in the model, one facet per QTL. For multivariate analyses, one line
#'       per trial is shown within each facet. Large jumps suggest confounding.}
#'     \item{\code{"both"}}{Returns a named list
#'       \code{list(lrt = ..., stability = ...)} invisibly.}
#'   }
#' @param sig.col Colour for highlighted points and ribbons in diagnostic
#'   plots. Default \code{"firebrick"}.
#' @param \dots Further arguments (e.g.\ \code{digits}) passed to methods.
#'
#' @return \code{NULL} invisibly when \code{plot = FALSE}; a \code{ggplot}
#'   when \code{plot} is \code{"lrt"} or \code{"stability"}; a named list
#'   \code{list(lrt, stability)} when \code{plot = "both"}.
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{gwasAim}},
#'   \code{\link{summary.qtlAim}}, \code{\link{summary.gwasAim}}
#'
#' @examples
#' \dontrun{
#' # After running qtlAim():
#' aimTrace(qtl.fit)
#' aimTrace(qtl.fit, plot = "lrt")
#' aimTrace(qtl.fit, plot = "stability")
#' }
#'
#' @export
aimTrace <- function(object, ...) UseMethod("aimTrace")

#' @rdname aimTrace
#' @exportS3Method
aimTrace.qtlAim <- function(object,
                             iter    = seq_along(object$QTL$diag$coef.list),
                             lik.out = TRUE,
                             plot    = FALSE,
                             sig.col = "firebrick",
                             ...) {
    dots <- list(...)
    dig  <- if (!is.na(pmatch("digits", names(dots)))) dots$digits else options()$digits

    cl <- object$QTL$diag$coef.list
    vl <- object$QTL$diag$vcoef.list
    sigma2 <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1

    is.mv  <- !is.null(object$QTL$Trait)
    n.qtl  <- length(object$QTL$qtl)
    # Column labels in forward detection order
    qnams  <- sub("^Chr\\.", "", object$QTL$qtl)

    # Build z-ratio list: one element per iteration, one value per unique QTL.
    # For multivariate, extract only the main-effect coefficient (no ":" in name)
    # so the matrix has one column per QTL regardless of ntrait.
    if (is.mv) {
        zrl <- lapply(seq_along(cl), function(i) {
            cof  <- cl[[i]]; vcof <- vl[[i]]
            main <- !grepl(":", names(cof))
            cof[main] / sqrt(vcof[main] * sigma2)
        })
    } else {
        zrl <- lapply(seq_along(cl), function(i)
            cl[[i]] / sqrt(vl[[i]] * sigma2))
    }

    if (any(ret <- is.na(pmatch(iter, seq_along(zrl))))) {
        warning("\"iter\" values outside expected range. Using values within range.")
        iter <- iter[!ret]
    }

    # p-values -- coef.list stores QTL in reversed detection order (most recent
    # first).  Reverse each z-vector so column k = QTL k (detection order).
    if (object$QTL$method == "random") {
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- round((1 - pchisq(rev(el)^2, df = 1)) / 2, dig)
            c(pv, rep(NA, len - length(pv)))
        }, len = n.qtl, dig = dig)
    } else {
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- round(2 * (1 - pnorm(abs(rev(el)))), dig)
            c(pv, rep(NA, len - length(pv)))
        }, len = n.qtl, dig = dig)
    }

    qtlmat <- do.call("rbind", pvals)
    dimnames(qtlmat) <- list(paste0("Iter.", seq_along(zrl)), qnams)

    cat("\nIncremental QTL P-value Matrix")
    if (is.mv) cat("  (main effect per QTL)")
    cat("\n===============================\n")
    qtlmat[qtlmat < 0.001] <- "<0.001"
    qtlmat[is.na(qtlmat)]  <- ""
    print.default(qtlmat[iter, seq_len(iter[length(iter)]), drop = FALSE],
                  quote = FALSE, right = TRUE, ...)

    if (lik.out) {
        cat("\nLikelihood Ratio Test of Additive Variance Parameter.\n")
        cat("======================================================\n")
        dmat <- round(as.matrix(object$QTL$diag$lik.mat), dig)
        dimnames(dmat)[[1]] <- paste0("Iter.", seq_len(nrow(dmat)))
        dmat[, 4][dmat[, 4] < 0.001] <- "<0.001"
        print.default(dmat, quote = FALSE, right = TRUE, ...)
    }

    if (identical(plot, FALSE)) return(invisible(NULL))
    plot <- match.arg(as.character(plot), c("lrt", "stability", "both"))
    if (plot == "lrt")       return(.plot_lrt(object, sig.col))
    if (plot == "stability") return(.plot_stability(object, sig.col))
    invisible(list(lrt       = .plot_lrt(object, sig.col),
                   stability = .plot_stability(object, sig.col)))
}

# =============================================================================
# Shared engine: .plot_lrt
#
# LRT statistic trace across all forward-selection iterations.
#
#   x-axis    : iteration number (including the final non-significant one)
#   y-axis    : LRT statistic
#   Line      : grey connector between all points
#   Points    : filled sig.col for passing, open grey for the failed iteration
#   Threshold : dashed horizontal line; for ntrait == 1 uses qchisq boundary
#               test threshold; for ntrait > 1 uses qchisq.mixture().
#   Labels    : short QTL name above each passing point via geom_text_repel
# =============================================================================
.plot_lrt <- function(object, sig.col = "firebrick") {

    lik    <- object$QTL$diag$lik.mat
    n_rows <- nrow(lik)
    n_qtl  <- length(object$QTL$qtl)
    TypeI  <- if (!is.null(object$QTL$TypeI)) object$QTL$TypeI else 0.05

    is.mv  <- !is.null(object$QTL$Trait)
    ntrait <- if (is.mv) length(object$QTL$trait.levels) else 1L
    thresh <- if (ntrait == 1L)
        qchisq(1 - 2 * TypeI, 1)
    else
        qchisq.mixture(1 - TypeI, ntrait)

    qtl_labels <- sub("^Chr\\.", "", object$QTL$qtl)

    df <- data.frame(
        iter   = seq_len(n_rows),
        stat   = lik[, "Statistic"],
        passed = seq_len(n_rows) <= n_qtl,
        label  = c(qtl_labels, rep("", n_rows - n_qtl)),
        stringsAsFactors = FALSE
    )

    gp <- ggplot2::ggplot(df,
               ggplot2::aes(x = .data$iter, y = .data$stat)) +
        ggplot2::geom_hline(yintercept = thresh,
                            linetype  = "dashed",
                            colour    = "grey40",
                            linewidth = 0.6) +
        ggplot2::annotate("text",
                          x      = max(df$iter),
                          y      = thresh,
                          label  = sprintf("threshold  (TypeI = %.4f)", TypeI),
                          hjust  = 1, vjust = -0.5,
                          size   = 3, colour = "grey40") +
        ggplot2::geom_line(colour = "grey60", linewidth = 0.7) +
        ggplot2::geom_point(
            ggplot2::aes(fill = .data$passed),
            shape = 21, size = 3.5, stroke = 0.7,
            colour = "grey30", show.legend = FALSE) +
        ggplot2::scale_fill_manual(
            values = c("TRUE" = sig.col, "FALSE" = "grey85")) +
        ggrepel::geom_text_repel(
            ggplot2::aes(label = .data$label),
            size               = 3,
            colour             = sig.col,
            nudge_y            = diff(range(df$stat)) * 0.06,
            min.segment.length = 0,
            show.legend        = FALSE) +
        ggplot2::scale_x_continuous(breaks = seq_len(n_rows)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
        ggplot2::xlab("Iteration") +
        ggplot2::ylab("LRT Statistic") +
        ggplot2::coord_cartesian(clip = "off") +
        theme_scatter()
    gp
}

# =============================================================================
# Shared engine: .build_stability_df
#
# Builds a long-format data frame for the effect stability plot.
#
# Univariate (ntrait == 1):
#   For each detected QTL k, extracts its main effect estimate and +/-1 SE at
#   every iteration j >= k.  coef.list[[j]] is in reversed detection order so
#   QTL k's coefficient is at position (j - k + 1) from the end, but we match
#   by name ("X.Chr1.20") to be safe.
#
# Multivariate (ntrait > 1):
#   For each QTL k, at each iteration j >= k:
#     - Main effect (Trial1 reference): coefficient named "X.Chr1.20"
#     - Per-trial effects for Trial_t (t > 1):
#         beta_t = main + (Trial_t interaction coefficient)
#         SE_t   = sqrt(se_main^2 + se_int_t^2)  [approx; ignores covariance]
#   Returns one row per (QTL, iteration, trial).
# =============================================================================
.build_stability_df <- function(object) {

    cl     <- object$QTL$diag$coef.list
    vl     <- object$QTL$diag$vcoef.list
    n_qtl  <- length(object$QTL$qtl)
    is.mv  <- !is.null(object$QTL$Trait)

    sigma2 <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1

    qtl.x.keys <- sub("Chr\\.", "X.", object$QTL$qtl)  # "X.Chr1.20" format

    if (!is.mv) {
        # Univariate: match main coefficient by name, one row per (QTL, iter)
        rows <- lapply(seq_len(n_qtl), function(k) {
            label <- sub("^Chr\\.", "", object$QTL$qtl[k])
            xk    <- qtl.x.keys[k]
            do.call(rbind, lapply(k:n_qtl, function(j) {
                idx <- which(names(cl[[j]]) == xk)
                if (length(idx) == 0L) return(NULL)
                eff <- cl[[j]][idx]
                se  <- sqrt(vl[[j]][idx] * sigma2)
                data.frame(qtl_label = label, qtl_idx = k, iter = j,
                           effect = eff, se = se,
                           lo = eff - se, hi = eff + se,
                           trial = "univariate",
                           stringsAsFactors = FALSE)
            }))
        })

    } else {
        # Multivariate: reconstruct per-trial effects at each iteration
        trait.levels <- object$QTL$trait.levels

        rows <- lapply(seq_len(n_qtl), function(k) {
            label <- sub("^Chr\\.", "", object$QTL$qtl[k])
            xk    <- qtl.x.keys[k]

            do.call(rbind, lapply(k:n_qtl, function(j) {
                cof  <- cl[[j]];  vcof <- vl[[j]]

                # Main effect (Trial1 reference level)
                main.idx <- which(names(cof) == xk)
                if (length(main.idx) == 0L) return(NULL)
                main.eff <- cof[main.idx]
                main.var <- vcof[main.idx] * sigma2

                # One row per trial level
                do.call(rbind, lapply(seq_along(trait.levels), function(t) {
                    tname <- trait.levels[t]
                    if (t == 1L) {
                        # Trial 1 is the reference: effect = main coefficient
                        eff <- main.eff
                        se  <- sqrt(main.var)
                    } else {
                        # Trial t: effect = main + interaction deviation
                        int.nm  <- paste0(tname, ":", xk)
                        int.idx <- which(names(cof) == int.nm)
                        if (length(int.idx) == 0L) {
                            eff <- main.eff
                            se  <- sqrt(main.var)
                        } else {
                            int.var <- vcof[int.idx] * sigma2
                            eff <- main.eff + cof[int.idx]
                            # Approximate SE: ignores main/interaction covariance
                            se  <- sqrt(main.var + int.var)
                        }
                    }
                    data.frame(qtl_label = label, qtl_idx = k, iter = j,
                               effect = eff, se = se,
                               lo = eff - se, hi = eff + se,
                               trial = tname,
                               stringsAsFactors = FALSE)
                }))
            }))
        })
    }

    rows <- rows[!sapply(rows, is.null)]
    do.call(rbind, rows)
}

# =============================================================================
# Shared engine: .plot_stability
#
# QTL effect stability plot.
#
# Univariate: one facet per QTL (3 columns), one line, lollipop + error bars.
#
# Multivariate: same facet layout but one line per trial level, coloured with
#   hcl.colors("Dark 2") palette to match the blups plot.  Facet labels carry
#   a [MAIN] or [INT] suffix (when is.interaction is available).  The SE is
#   approximate (ignores main/interaction coefficient covariance) so error bars
#   are labelled accordingly in the y-axis title.
# =============================================================================
.plot_stability <- function(object, sig.col = "firebrick") {

    sdf   <- .build_stability_df(object)
    n_qtl <- length(object$QTL$qtl)
    is.mv <- !is.null(object$QTL$Trait)

    # Build facet labels (with optional MAIN/INT tag for multivariate)
    base_labels <- sub("^Chr\\.", "", object$QTL$qtl)
    if (is.mv && !is.null(object$QTL$is.interaction)) {
        facet_labels <- paste0(
            base_labels, "  [",
            ifelse(object$QTL$is.interaction, "INT", "MAIN"), "]")
    } else {
        facet_labels <- base_labels
    }
    label_map         <- stats::setNames(facet_labels, base_labels)
    sdf$qtl_label_fac <- factor(label_map[sdf$qtl_label],
                                levels = facet_labels)

    all_iters <- seq_len(max(sdf$iter))

    if (is.mv) {
        trial.levels <- object$QTL$trait.levels
        trial.cols   <- grDevices::hcl.colors(length(trial.levels), "Dark 2")
        names(trial.cols) <- trial.levels
        sdf$trial <- factor(sdf$trial, levels = trial.levels)

        gp <- ggplot2::ggplot(sdf,
                   ggplot2::aes(x      = .data$iter,
                                y      = .data$effect,
                                colour = .data$trial,
                                group  = .data$trial)) +
            ggplot2::facet_wrap(~ qtl_label_fac, scales = "free_y", ncol = 3) +

            ggplot2::geom_hline(yintercept = 0,
                                linetype = "dashed", colour = "grey60",
                                linewidth = 0.5) +

            ggplot2::geom_errorbar(
                ggplot2::aes(ymin = .data$lo, ymax = .data$hi,
                             colour = .data$trial),
                width = 0.15, linewidth = 0.5) +

            ggplot2::geom_line(linewidth = 0.5, linetype = "dashed") +

            ggplot2::geom_point(
                ggplot2::aes(fill = .data$trial),
                shape = 21, size = 2.5, stroke = 0.5,
                colour = "white", show.legend = TRUE) +

            ggplot2::scale_colour_manual(values = trial.cols,
                                         name   = object$QTL$Trait) +
            ggplot2::scale_fill_manual(values   = trial.cols,
                                       name     = object$QTL$Trait) +
            ggplot2::scale_x_continuous(breaks = all_iters) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.18, 0.18))) +
            ggplot2::xlab("Iteration") +
            ggplot2::ylab("Per-trial Effect \u00b1 approx. 1 SE") +
            theme_scatter() +
            ggplot2::theme(
                panel.spacing = ggplot2::unit(0.8, "lines"),
                strip.text    = ggplot2::element_text(size = 8, face = "bold"),
                legend.position = "right")

    } else {
        gp <- ggplot2::ggplot(sdf,
                   ggplot2::aes(x = .data$iter, y = .data$effect)) +
            ggplot2::facet_wrap(~ qtl_label_fac, scales = "free_y", ncol = 3) +

            ggplot2::geom_hline(yintercept = 0,
                                linetype = "dashed", colour = "grey60",
                                linewidth = 0.5) +

            ggplot2::geom_errorbar(
                ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                colour = sig.col, width = 0.15, linewidth = 0.7) +

            ggplot2::geom_line(colour = "grey50", linewidth = 0.5,
                               linetype = "dashed") +

            ggplot2::geom_point(fill = sig.col, colour = "white",
                                shape = 21, size = 3, stroke = 0.6) +

            ggplot2::scale_x_continuous(breaks = all_iters) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.18, 0.18))) +
            ggplot2::xlab("Iteration") +
            ggplot2::ylab("Effect \u00b1 1 SE") +
            theme_scatter() +
            ggplot2::theme(
                panel.spacing = ggplot2::unit(0.8, "lines"),
                strip.text    = ggplot2::element_text(size = 8, face = "bold"))
    }
    gp
}
