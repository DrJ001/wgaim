# =============================================================================
# waldTest.R
# Wald test generic and asreml method for testing linear hypotheses on fixed
# effects from a fitted ASReml model.
#
# Used internally by qtlAim.asreml() and gwasAim.asreml() after the
# forward-selection loop to determine, for each detected multivariate QTL,
# whether the Trait x QTL interaction term is significant (interaction QTL)
# or whether it collapses to a main-effect-only QTL.
#
# Also exported for direct user access.
# =============================================================================

#' Wald Test for Linear Hypotheses on Fixed Effects
#'
#' @description
#' Computes Wald test statistics for linear hypotheses on the fixed effects of
#' a fitted \code{asreml} model. Supports two test types: treatment contrasts
#' (\code{type = "con"}) and zero-equality tests (\code{type = "zero"}).
#'
#' The zero-equality test is used internally by \code{\link{qtlAim}} and
#' \code{\link{gwasAim}} when \code{Trait} is non-\code{NULL} to determine
#' whether each detected QTL has a significant Trait-by-QTL interaction or
#' whether the interaction coefficients are jointly zero (main effect only).
#'
#' @param object A fitted model object.
#' @param \dots Further arguments passed to methods.
#'
#' @return A list with components \code{$Contrasts} (data frame, or
#'   \code{NULL} if no contrast tests requested) and \code{$Zero} (data frame,
#'   or \code{NULL} if no zero-equality tests requested).
#'
#' @seealso \code{\link{qtlAim}}, \code{\link{gwasAim}}
#'
#' @name waldTest
#' @export
waldTest <- function(object, ...)
    UseMethod("waldTest")

#' @rdname waldTest
#' @exportS3Method
waldTest.default <- function(object, ...)
    stop("waldTest is currently only implemented for \"asreml\" objects.")

#' @rdname waldTest
#'
#' @param cc A named list of comparison objects. Each element must be a list
#'   with at minimum a \code{coef} component (numeric indices or coefficient
#'   names) and a \code{type} component (\code{"con"} for treatment contrast
#'   or \code{"zero"} for zero-equality test). Contrast tests also require a
#'   \code{comp} component (numeric contrast weights). An optional
#'   \code{group} component (list with \code{left} and \code{right} elements)
#'   provides custom row labels.
#' @param keep.fac Logical. If \code{TRUE} (default), factor-level prefixes
#'   are retained in contrast row labels.
#'
#' @exportS3Method
waldTest.asreml <- function(object, cc, keep.fac = TRUE, ...) {
    if (!inherits(object, "asreml"))
        stop("waldTest.asreml requires an object of class \"asreml\".")
    if (is.null(object$Cfixed)) {
        warning("Cfixed matrix not found. Refitting with asreml.options(Cfixed = TRUE).\n")
        asreml.options(Cfixed = TRUE)
        object <- update(object)
    }

    vrb    <- as.matrix(object$Cfixed)
    tau    <- c(object$coefficients$fixed)
    names(tau) <- rownames(object$coefficients$fixed)
    nc     <- length(tau)
    sigma2 <- object$sigma2
    vrb    <- vrb / sigma2
    ccnams <- names(tau)

    # ------------------------------------------------------------------
    # Resolve coefficient identifiers in cc: numeric indices or names
    # ------------------------------------------------------------------
    cc <- lapply(cc, function(el, ccnams) {
        if (!all(names(el) %in% c("coef", "type", "comp", "group")))
            stop("Each element of cc must have names from: \"coef\", \"type\", \"comp\", \"group\".")
        if (is.numeric(el$coef)) {
            if (max(el$coef) > length(ccnams))
                stop("Coefficient subscript out of bounds.")
            names(el$coef) <- ccnams[el$coef]
        } else {
            idx <- pmatch(el$coef, ccnams)
            if (any(is.na(idx)))
                stop("Names in contrast do not match coefficient names of object.")
            names(idx) <- el$coef
            el$coef    <- idx
        }
        el
    }, ccnams)

    ctype <- vapply(cc, function(el) el$type, character(1L))
    if (!all(ctype %in% c("con", "zero")))
        stop("Contrast types must be \"con\" (treatment contrast) or \"zero\" (zero-equality test).")

    cons <- cc[ctype == "con"]
    zero <- cc[ctype == "zero"]

    cdf <- zdf <- NULL

    # ------------------------------------------------------------------
    # Treatment contrasts
    # ------------------------------------------------------------------
    if (length(cons)) {
        CRows <- lapply(cons, function(el, nc) {
            if (length(el) < 3L) {
                con       <- contr.helmert(length(el$coef))[, length(el$coef) - 1L]
                names(con) <- names(el$coef)
                cat("Warning: default Helmert contrast taken for",
                    paste(names(el$coef), collapse = ", "), ":", con, "\n")
                row         <- rep(0, nc)
                row[el$coef] <- con
                row
            } else if (is.matrix(el$comp)) {
                if (length(el$coef) != ncol(el$comp))
                    stop("Contrast matrix columns must match length of coef.")
                do.call(rbind, lapply(seq_len(nrow(el$comp)), function(r) {
                    row          <- rep(0, nc)
                    row[el$coef] <- el$comp[r, ]
                    row
                }))
            } else {
                if (length(el$coef) != length(el$comp))
                    stop("Length of contrast vector must match length of coef.")
                row          <- rep(0, nc)
                row[el$coef] <- el$comp
                row
            }
        }, nc)
        Cmat <- do.call(rbind, CRows)

        if (!keep.fac)
            ccnams <- substring(ccnams, regexpr("_", ccnams) + 1L, nchar(ccnams))

        # Build row labels from contrast signs, overriding with group names if given
        cnam <- do.call(rbind, lapply(split(Cmat, seq_len(nrow(Cmat))), function(row, nm) {
            c(paste(nm[row > 0], collapse = ":"), paste(nm[row < 0], collapse = ":"))
        }, nm = ccnams))

        gnams <- do.call(rbind, lapply(cons, function(el) {
            nr <- if (is.matrix(el$comp)) nrow(el$comp) else 1L
            gl <- if (!is.null(el$group$left))  rep_len(el$group$left,  nr) else rep(NA, nr)
            gr <- if (!is.null(el$group$right)) rep_len(el$group$right, nr) else rep(NA, nr)
            cbind(gl, gr)
        }))
        cnam[!is.na(gnams[, 1L]), 1L] <- gnams[!is.na(gnams[, 1L]), 1L]
        cnam[!is.na(gnams[, 2L]), 2L] <- gnams[!is.na(gnams[, 2L]), 2L]

        cse <- ctau <- cwtest <- numeric(nrow(Cmat))
        for (i in seq_len(nrow(Cmat))) {
            varmat   <- sum(Cmat[i, ] * crossprod(vrb, t(Cmat)[, i]))
            cse[i]   <- sqrt(varmat * sigma2)
            ctau[i]  <- sum(Cmat[i, ] * tau)
            cwtest[i] <- (ctau[i] / cse[i])^2
        }
        cdf <- data.frame(
            "Wald Statistic" = round(cwtest, 6),
            "P-Value"        = round(1 - pchisq(cwtest, 1L), 6),
            "Cont. Coef."    = round(ctau, 6),
            "Std. Error"     = round(cse, 6),
            check.names = FALSE,
            row.names   = paste(cnam[, 1L], cnam[, 2L], sep = " vs ")
        )
    }

    # ------------------------------------------------------------------
    # Zero-equality tests
    # ------------------------------------------------------------------
    if (length(zero)) {
        ZRows <- lapply(zero, function(el, nc) {
            nr   <- length(el$coef)
            rows <- rep(rep(0, nc), nr)
            dum  <- seq(0L, (nr - 1L) * nc, by = nc)
            rows[el$coef + dum] <- 1
            matrix(rows, nrow = nr, byrow = TRUE)
        }, nc)

        znam <- vapply(zero, function(el, ccnams) {
            if (is.null(el$group)) paste(ccnams[el$coef], collapse = ":")
            else el$group
        }, character(1L), ccnams = ccnams)

        if (any(table(znam) > 1L))
            stop("Duplicate names found in group structures for zero-equality tests.")

        zwtest <- zpval <- numeric(length(ZRows))
        for (i in seq_along(ZRows)) {
            varmat    <- ZRows[[i]] %*% crossprod(vrb, t(ZRows[[i]]))
            Ctau      <- ZRows[[i]] %*% tau
            zwtest[i] <- sum(Ctau * crossprod(solve(varmat), Ctau)) / sigma2
            zpval[i]  <- 1 - pchisq(zwtest[i], nrow(ZRows[[i]]))
        }
        zdf <- data.frame(
            "Wald Statistic" = round(zwtest, 6),
            "P-Value"        = round(zpval, 6),
            check.names = FALSE,
            row.names   = znam
        )
    }

    invisible(list(Contrasts = cdf, Zero = zdf))
}
