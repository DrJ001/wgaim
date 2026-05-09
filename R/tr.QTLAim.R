# =============================================================================
# tr.QTLAim.R
# S3 tr() method for QTLAim objects.
# Prints incremental QTL p-value matrix and LRT table.
# Note: the tr() generic is defined in wgaim16.R and is already available.
# =============================================================================

tr.QTLAim <- function(object, iter = 1:length(object$QTL$effects),
                      lik.out = TRUE, ...) {
    dots <- list(...)
    dig <- if (!is.na(pmatch("digits", names(dots)))) dots$digits else options()$digits
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
        cat("====================================================\n")
        dmat <- round(as.matrix(object$QTL$diag$lik.mat), dig)
        dimnames(dmat)[[1]] <- paste0("Iter.", 1:(length(zrl) + 1))
        dmat[, 4][dmat[, 4] < 0.001] <- "<0.001"
        print.default(dmat, quote = FALSE, right = TRUE, ...)
    }
}
