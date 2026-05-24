# =============================================================================
# engine_results.R
# Internal function to package analysis results into a structured list.
# =============================================================================

#' @keywords internal
#'
#' Trait     : character name of the Trait factor (NULL = univariate).
#' qtlModel  : the final fitted asreml model (needed for multivariate waldTest
#'             pruning). Only inspected when Trait is non-NULL and QTL were
#'             found; passed through unchanged for the univariate case.
#'
#' When Trait is non-NULL and QTL were found, a waldTest() zero-equality test
#' is performed per QTL to determine whether the Trait:QTL interaction is
#' significant. If not (p > TypeI), the QTL is retained as a main effect only
#' and the interaction term is dropped from the model. The pruned model is
#' returned as qtlModel.pruned in the output list.
.packResults <- function(qtl, coef.list, vcoef.list, ldiag, state, iter,
                         breakout, cov.env, genetic.term, method, gen.type,
                         selection, TypeI, Trait = NULL, qtlModel = NULL,
                         trait.levels = NULL, phenoData = NULL) {
    qtl.list <- list()
    qtl.list$selection <- selection
    qtl.list$method    <- method
    qtl.list$type      <- gen.type
    qtl.list$TypeI     <- TypeI
    qtl.list$diag      <- ldiag
    qtl.list$iterations <- iter
    qtlModel.pruned    <- qtlModel   # returned even if no pruning needed

    if (length(qtl)) {
        qtl.list$diag$coef.list  <- coef.list
        qtl.list$diag$vcoef.list <- vcoef.list
        qtl.list$diag$lik.mat    <- matrix(unlist(ldiag$lik), ncol = 4, byrow = TRUE)
        dimnames(qtl.list$diag$lik.mat)[[2]] <- c("L0", "L1", "Statistic", "Pvalue")
        qtl.list$diag$state       <- state
        qtl.list$diag$genetic.term <- genetic.term
        qtl.list$diag$rel.scale   <- ifelse(!is.null(cov.env), cov.env$scale, 1)
        qtl.list$breakout <- breakout != -1
        qtl.list$qtl      <- qtl
        qtl.list$effects  <- coef.list[[iter]]
        qtl.list$veffects <- vcoef.list[[iter]]

        # ------------------------------------------------------------------
        # Multivariate post-loop: waldTest pruning of Trait:QTL interactions
        # ------------------------------------------------------------------
        if (!is.null(Trait) && !is.null(qtlModel)) {
            list.coefs <- qtlModel$coefficients$fixed
            trms       <- attr(terms(qtlModel$call$fixed), "term.labels")
            # All fixed terms involving a QTL column (X.chr.idx)
            marks  <- trms[grep("X\\.", trms)]
            imarks <- marks[grep(":", marks)]    # interaction terms: Trait:X.chr.idx
            mmarks <- marks[!grepl(":", marks)]  # main-effect terms: X.chr.idx

            int.test <- vector("list", length(imarks))
            for (i in seq_along(imarks)) {
                # Identify fixed coefficients for this interaction term
                forw <- paste0(Trait, ".*", mmarks[i])
                reve <- paste0(mmarks[i], ".*", Trait)
                zind <- grep(paste(forw, reve, sep = "|"), rownames(list.coefs))
                ci   <- list.coefs[zind, 1L]
                # Test only the non-aliased (non-zero reference) coefficients
                int.test[[i]] <- list(coef = zind[ci != 0], type = "zero")
            }
            wt <- .waldTest(qtlModel, cc = int.test)

            # Decide per QTL: keep interaction or reduce to main effect
            keep.int    <- wt[, "P-Value"] <= TypeI
            final.terms <- ifelse(keep.int, imarks, mmarks)
            other.terms <- trms[!grepl("X\\.", trms)]
            fix.form    <- as.formula(paste(
                ". ~", paste(c(other.terms, final.terms), collapse = " + ")))

            cat("\nWald Test: Pruning Trait x QTL interaction terms:\n")
            cat("==================================================\n")
            print(wt)
            cat("\n")

            qtlModel.pruned <- update(qtlModel, fixed. = fix.form)

            # Store multivariate-specific slots
            qtl.list$Trait        <- Trait
            qtl.list$trait.levels <- trait.levels
            qtl.list$wald.test    <- wt
            qtl.list$final.terms  <- final.terms
            qtl.list$is.interaction <- keep.int
        }
    }
    list(qtl.list = qtl.list, qtlModel.pruned = qtlModel.pruned)
}

# =============================================================================
# .waldTest -- internal zero-equality Wald test for multivariate QTL pruning.
#
# Tests whether a set of fixed-effect coefficients (identified by numeric
# index in cc) are jointly zero, using the Cfixed variance matrix from a
# fitted ASReml model.  Used by .packResults() to decide whether each
# detected Trait:QTL interaction is significant or reduces to a main effect.
#
# Arguments:
#   object : fitted asreml model with $Cfixed, $coefficients$fixed, $sigma2
#   cc     : named list, each element a list with:
#              $coef  -- integer indices of the coefficients to test
#              $group -- optional character label for the output row
# Returns:
#   data.frame with columns "Wald Statistic" and "P-Value", one row per
#   element of cc.
# =============================================================================
#' @keywords internal
.waldTest <- function(object, cc) {
    if (is.null(object$Cfixed)) {
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

    # Resolve any name-based coef references to integer indices
    cc <- lapply(cc, function(el) {
        if (!is.numeric(el$coef)) {
            idx <- pmatch(el$coef, ccnams)
            if (any(is.na(idx)))
                stop("Names in cc$coef do not match coefficient names of model.")
            el$coef <- idx
        }
        if (max(el$coef) > nc)
            stop("Coefficient subscript out of bounds.")
        el
    })

    # Build row labels
    znam <- vapply(cc, function(el) {
        if (!is.null(el$group)) el$group
        else paste(ccnams[el$coef], collapse = ":")
    }, character(1L))

    if (any(table(znam) > 1L))
        stop("Duplicate group names in zero-equality tests.")

    # Zero-equality Wald statistics
    zwtest <- zpval <- numeric(length(cc))
    for (i in seq_along(cc)) {
        coef_idx <- cc[[i]]$coef
        nr       <- length(coef_idx)
        rows     <- rep(rep(0, nc), nr)
        dum      <- seq(0L, (nr - 1L) * nc, by = nc)
        rows[coef_idx + dum] <- 1
        Zmat      <- matrix(rows, nrow = nr, byrow = TRUE)
        varmat    <- Zmat %*% crossprod(vrb, t(Zmat))
        Ctau      <- Zmat %*% tau
        zwtest[i] <- sum(Ctau * crossprod(solve(varmat), Ctau)) / sigma2
        zpval[i]  <- 1 - pchisq(zwtest[i], nr)
    }

    data.frame(
        "Wald Statistic" = round(zwtest, 6),
        "P-Value"        = round(zpval,  6),
        check.names = FALSE,
        row.names   = znam
    )
}
