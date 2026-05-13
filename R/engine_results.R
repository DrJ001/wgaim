# =============================================================================
# engine_results.R
# Internal function to package analysis results into a structured list.
# =============================================================================

#' @keywords internal
.packResults <- function(qtl, coef.list, vcoef.list, ldiag, state, iter,
                         breakout, cov.env, genetic.term, method, gen.type,
                         selection, TypeI) {
    qtl.list <- list()
    qtl.list$selection <- selection
    qtl.list$method <- method
    qtl.list$type <- gen.type
    qtl.list$TypeI <- TypeI
    qtl.list$diag <- ldiag
    qtl.list$iterations <- iter
    if (length(qtl)) {
        qtl.list$diag$coef.list <- coef.list
        qtl.list$diag$vcoef.list <- vcoef.list
        qtl.list$diag$lik.mat <- matrix(unlist(ldiag$lik), ncol = 4, byrow = TRUE)
        dimnames(qtl.list$diag$lik.mat)[[2]] <- c("L0", "L1", "Statistic", "Pvalue")
        qtl.list$diag$state <- state
        qtl.list$diag$genetic.term <- genetic.term
        qtl.list$diag$rel.scale <- ifelse(!is.null(cov.env), cov.env$scale, 1)
        qtl.list$breakout <- breakout != -1
        qtl.list$qtl <- qtl
        qtl.list$effects <- coef.list[[iter - 1]]
        qtl.list$veffects <- vcoef.list[[iter - 1]]
    }
    qtl.list
}
