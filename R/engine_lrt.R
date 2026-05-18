# =============================================================================
# engine_lrt.R
# Internal LRT helpers for significance testing in the wgAim engine.
#
# Univariate (ntrait == 1):
#   Standard one-sided boundary LRT on a single variance parameter.
#   Null distribution is a 50:50 mixture of chi^2(0) and chi^2(1),
#   giving a halved p-value: (1 - pchisq(stat, 1)) / 2.
#
# Multivariate (ntrait > 1):
#   Boundary LRT on a covariance matrix parameter. The null distribution
#   is a mixture of chi^2(k) for k = 0, ..., ntrait, with binomial
#   mixing weights (Stram & Lee 1994; Self & Liang 1987).
#   pchisq.mixture() and qchisq.mixture() implement this mixture.
#
# .compute_vqtilde():
#   Pure-R implementation of the multivariate posterior variance diagonal
#   needed by .qtlSelect() when ntrait > 1. For ntrait == 1 the simpler
#   scalar computation in .qtlSelect() is used instead.
# =============================================================================

# -----------------------------------------------------------------------------
# pchisq.mixture
# CDF of the mixture chi-squared distribution for boundary LRT on an
# ntrait-dimensional covariance matrix parameter.
#
# x      : numeric vector of test statistics
# ntrait : number of traits (>= 1). ntrait == 1 reduces to standard
#          one-sided chi^2(1) boundary test.
# -----------------------------------------------------------------------------
#' @keywords internal
pchisq.mixture <- function(x, ntrait = 2L) {
    df       <- 0:ntrait
    mixprobs <- dbinom(df, size = ntrait, prob = 0.5)
    vapply(x, function(xi) {
        sum(mixprobs * pchisq(xi, df = df))
    }, numeric(1L))
}

# -----------------------------------------------------------------------------
# qchisq.mixture
# Quantile (percentile) of the mixture chi-squared distribution.
# Uses Newton-Raphson iteration starting from the plain chi^2(ntrait) quantile.
#
# prob    : probability (scalar, 0 < prob < 1)
# ntrait  : number of traits (>= 1)
# maxit   : maximum Newton-Raphson iterations
# tol     : convergence tolerance on the Newton correction
# -----------------------------------------------------------------------------
#' @keywords internal
qchisq.mixture <- function(prob, ntrait = 2L, tol = .Machine$double.eps^0.5) {
    # Special case: ntrait == 1 reduces to standard one-sided boundary test.
    if (ntrait == 1L)
        return(qchisq(prob, df = 1L))

    # Use uniroot rather than Newton-Raphson.  The N-R approach diverges when
    # ntrait is large (e.g. 8) because the starting value qchisq(prob, ntrait)
    # places the mixture CDF already very close to 1; the derivative (chi-sq
    # density) is near-zero there, so the Newton step overshoots to a large
    # negative value where the density is exactly 0, giving NaN.
    #
    # The root is bracketed in [0, qchisq(prob, ntrait)]:
    #   lower: CDF at 0+ = dbinom(0, ntrait, 0.5) = 2^{-ntrait} << prob
    #   upper: CDF at qchisq(prob, ntrait) >= prob (shown above)
    uniroot(
        function(cv) pchisq.mixture(cv, ntrait = ntrait) - prob,
        lower = 0,
        upper = qchisq(prob, df = ntrait),
        tol   = tol
    )$root
}

# -----------------------------------------------------------------------------
# .compute_vqtilde
# Compute the diagonal of  trans %*% vatilde %*% t(trans)  in a way that
# scales to ntrait > 1 without building the full (nlines*ntrait) x (nlines*ntrait)
# intermediate matrix explicitly.
#
# trans   : nmarkers x nlines transformation matrix (from .constructCM)
# Ginv    : ntrait x ntrait precision matrix (inverse of Ga)
# vatilde : (nlines*ntrait) x (nlines*ntrait) posterior variance of atilde
#           arranged in trait-major order (i.e. kronecker(Ga, relm) - vcov)
# ntrait  : number of traits
#
# Returns : numeric vector of length nmarkers
# -----------------------------------------------------------------------------
#' @keywords internal
.compute_vqtilde <- function(trans, Ginv, vatilde, ntrait) {
    nmarkers <- nrow(trans)
    nlines   <- ncol(trans)

    vapply(seq_len(nmarkers), function(i) {
        # Build the (ntrait * nlines) x ntrait block-row selector for marker i:
        # T_i = kronecker(I_ntrait, trans[i, ])  -> ntrait x (ntrait*nlines)
        Ti    <- kronecker(diag(ntrait), trans[i, , drop = FALSE])
        # tmp2 = Ti %*% vatilde %*% t(Ti)  : ntrait x ntrait
        tmp2  <- Ti %*% vatilde %*% t(Ti)
        # scalar = tr(Ginv %*% tmp2)
        sum(diag(Ginv %*% tmp2))
    }, numeric(1L))
}
