# =============================================================================
# aimTrace.qtlAim.R
# aimTrace() generic + S3 method for qtlAim objects.
# Shared internal helpers (.plot_lrt, .build_stability_df, .plot_stability)
# are also used by aimTrace.gwasAim (defined in aimTrace.gwasAim.R).
# =============================================================================

#' Trace the forward-selection algorithm
#'
#' Prints a p-value matrix and likelihood ratio test (LRT) table showing the
#' incremental state of the forward-selection algorithm across iterations,
#' and optionally returns diagnostic plots.
#'
#' @param object A fitted object of class \code{"qtlAim"} or
#'   \code{"gwasAim"}.
#' @param ... Further arguments passed to methods.
#' @export
aimTrace <- function(object, ...) UseMethod("aimTrace")

#' @describeIn aimTrace Trace method for \code{qtlAim} objects.
#'
#'   The \strong{p-value matrix} shows the p-value of each detected QTL at
#'   every iteration in which it was present in the model, making it easy to
#'   assess whether QTL effects remain stable as further QTL enter.
#'
#'   The \strong{LRT table} records the base and full model log-likelihoods,
#'   the LRT statistic, and its p-value at every iteration including the final
#'   non-significant one.
#'
#' @param iter Integer vector of iterations to include in the p-value matrix.
#'   Default is all iterations: \code{1:length(object\$QTL\$effects)}.
#' @param lik.out Logical. If \code{TRUE} (default), print the LRT table.
#' @param plot Controls optional diagnostic plot output (console output always
#'   printed regardless). One of:
#'   \describe{
#'     \item{\code{FALSE} (default)}{No plot returned.}
#'     \item{\code{"lrt"}}{Returns a \code{\link[ggplot2]{ggplot}} of the LRT
#'       statistic across iterations with the significance threshold marked and
#'       each detected QTL labelled.}
#'     \item{\code{"stability"}}{Returns a \code{\link[ggplot2]{ggplot}} of
#'       QTL effect estimates \eqn{\pm} 1 SE across every iteration in which
#'       the QTL was in the model, one facet per QTL. Instability (large jumps)
#'       flags possible confounding between QTL.}
#'     \item{\code{"both"}}{Returns \code{list(lrt = \ldots, stability =
#'       \ldots)} invisibly.}
#'   }
#' @param sig.col Colour for highlighted points and ribbons in the diagnostic
#'   plots. Default \code{"firebrick"}.
#' @return \code{NULL} invisibly when \code{plot = FALSE}; otherwise a
#'   \code{ggplot} object or named list of two \code{ggplot} objects.
#' @exportS3Method
aimTrace.qtlAim <- function(object,
                             iter    = 1:length(object$QTL$effects),
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

    zrl <- lapply(seq_along(cl), function(i)
        cl[[i]] / sqrt(vl[[i]] * sigma2))

    if (any(ret <- is.na(pmatch(iter, seq_along(zrl))))) {
        warning("\"iter\" values outside expected range. Using values within range.")
        iter <- iter[!ret]
    }

    if (object$QTL$method == "random") {
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- round((1 - pchisq(el^2, df = 1)) / 2, dig)
            c(pv, rep(NA, len - length(pv)))
        }, len = length(zrl), dig = dig)
    } else {
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- round(2 * (1 - pnorm(abs(el))), dig)
            c(pv, rep(NA, len - length(pv)))
        }, len = length(zrl), dig = dig)
    }

    qtlmat <- do.call("rbind", pvals)
    qnams  <- gsub("X\\.", "", names(cl))
    dimnames(qtlmat) <- list(paste0("Iter.", seq_along(zrl)), qnams)

    cat("\nIncremental QTL P-value Matrix.\n")
    cat("===============================\n")
    qtlmat[qtlmat < 0.001] <- "<0.001"
    qtlmat[is.na(qtlmat)]  <- ""
    print.default(qtlmat[iter, 1:iter[length(iter)], drop = FALSE],
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
#   Threshold : dashed horizontal line at qchisq(1 - 2*TypeI, 1), annotated
#   Labels    : short QTL name above each passing point via geom_text_repel
# =============================================================================
.plot_lrt <- function(object, sig.col = "firebrick") {

    lik    <- object$QTL$diag$lik.mat
    n_rows <- nrow(lik)
    n_qtl  <- length(object$QTL$qtl)
    TypeI  <- if (!is.null(object$QTL$TypeI)) object$QTL$TypeI else 0.05
    thresh <- qchisq(1 - 2 * TypeI, 1)

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
# For each detected QTL k, extracts its effect estimate and ±1 SE at every
# iteration j >= k (all iterations in which it was in the model).
# .addEffect() stores coefs with rev(), so within coef.list[[j]] (which has j
# elements) the detection order is reversed: element 1 = newest QTL (j),
# element j = first QTL (1).  To get QTL k's value at iteration j, use
# index (j - k + 1).
# =============================================================================
.build_stability_df <- function(object) {

    cl    <- object$QTL$diag$coef.list
    vl    <- object$QTL$diag$vcoef.list
    n_qtl <- length(object$QTL$qtl)

    sigma2 <- object$sigma2
    if (object$vparameters.con[length(object$vparameters.con)] == 4)
        sigma2 <- 1

    rows <- lapply(seq_len(n_qtl), function(k) {
        label <- sub("^Chr\\.", "", object$QTL$qtl[k])
        do.call(rbind, lapply(k:n_qtl, function(j) {
            idx <- j - k + 1L          # rev() offset: QTL k is at position j-k+1
            eff <- cl[[j]][idx]
            se  <- sqrt(vl[[j]][idx] * sigma2)
            data.frame(qtl_label = label,
                       qtl_idx   = k,
                       iter      = j,
                       effect    = eff,
                       se        = se,
                       lo        = eff - se,
                       hi        = eff + se,
                       stringsAsFactors = FALSE)
        }))
    })
    do.call(rbind, rows)
}

# =============================================================================
# Shared engine: .plot_stability
#
# QTL effect stability plot — lollipop + line style.
#
# For each detected QTL (one facet, free y-scale, 3 columns):
#   — dashed zero reference line
#   — vertical lollipop stem from 0 to the effect estimate at each iteration
#   — error bar (±1 SE) at the tip of each stem
#   — line connecting the effect estimates across iterations
#   — filled circle at each tip
# =============================================================================
.plot_stability <- function(object, sig.col = "firebrick") {

    sdf   <- .build_stability_df(object)
    n_qtl <- length(object$QTL$qtl)

    lev           <- unique(sdf$qtl_label[order(sdf$qtl_idx)])
    sdf$qtl_label <- factor(sdf$qtl_label, levels = lev)
    all_iters     <- seq_len(max(sdf$iter))

    gp <- ggplot2::ggplot(sdf,
               ggplot2::aes(x = .data$iter, y = .data$effect)) +
        ggplot2::facet_wrap(~ qtl_label, scales = "free_y", ncol = 3) +

        # Zero reference
        ggplot2::geom_hline(yintercept = 0,
                            linetype = "dashed", colour = "grey60",
                            linewidth = 0.5) +

        # ±1 SE error bars at the tip
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
            colour = sig.col, width = 0.15, linewidth = 0.7) +

        # Line connecting the tips across iterations
        ggplot2::geom_line(colour = "grey50", linewidth = 0.5,
                           linetype = "dashed") +

        # Filled circle at each tip
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
    gp
}
