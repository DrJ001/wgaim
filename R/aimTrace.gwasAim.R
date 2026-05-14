# =============================================================================
# aimTrace.gwasAim.R
# S3 aimTrace() method for gwasAim objects.
# Shared plot helpers (.plot_lrt, .plot_stability) live in aimTrace.qtlAim.R.
# =============================================================================

#' @rdname aimTrace
#' @exportS3Method
aimTrace.gwasAim <- function(object,
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

    # gwasAim always uses method = "fixed"
    zrl <- lapply(seq_along(cl), function(i)
        cl[[i]] / sqrt(vl[[i]] * sigma2))

    if (any(ret <- is.na(pmatch(iter, seq_along(zrl))))) {
        warning("\"iter\" values outside expected range. Using values within range.")
        iter <- iter[!ret]
    }

    pvals <- lapply(zrl, function(el, len, dig) {
        pv <- round(2 * (1 - pnorm(abs(el))), dig)
        c(pv, rep(NA, len - length(pv)))
    }, len = length(zrl), dig = dig)

    qtlmat <- do.call("rbind", pvals)
    qnams  <- gsub("X\\.", "", names(cl))
    dimnames(qtlmat) <- list(paste0("Iter.", seq_along(zrl)), qnams)

    cat(sprintf("\nGWAS  TypeI = %.4f  (%d markers in panel)\n",
                object$QTL$TypeI, object$QTL$n.markers))
    cat("\nIncremental Marker P-value Matrix.\n")
    cat("==================================\n")
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
