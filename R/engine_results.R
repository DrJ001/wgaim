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
        qtl.list$effects  <- coef.list[[iter - 1L]]
        qtl.list$veffects <- vcoef.list[[iter - 1L]]

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
            wt <- waldTest(qtlModel, cc = int.test)

            # Decide per QTL: keep interaction or reduce to main effect
            keep.int    <- wt$Zero[, "P-Value"] <= TypeI
            final.terms <- ifelse(keep.int, imarks, mmarks)
            other.terms <- trms[!grepl("X\\.", trms)]
            fix.form    <- as.formula(paste(
                ". ~", paste(c(other.terms, final.terms), collapse = " + ")))

            cat("\nWald Test: Pruning Trait x QTL interaction terms:\n")
            cat("==================================================\n")
            print(wt$Zero)
            cat("\n")

            qtlModel.pruned <- update(qtlModel, fixed. = fix.form)

            # Store multivariate-specific slots
            qtl.list$Trait        <- Trait
            qtl.list$trait.levels <- trait.levels
            qtl.list$wald.test    <- wt$Zero
            qtl.list$final.terms  <- final.terms
            qtl.list$is.interaction <- keep.int
        }
    }
    list(qtl.list = qtl.list, qtlModel.pruned = qtlModel.pruned)
}
