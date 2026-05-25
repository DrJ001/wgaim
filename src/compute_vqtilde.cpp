/**
 * Rcpp/RcppArmadillo implementation of compute_vqtilde for the
 * multivariate outlier statistic denominator in .qtlSelect().
 *
 * For each genetic marker k, computes
 *   vqtilde[k] = tr( Ginv %*% Var_k )
 * where
 *   Var_k[i,j] = trans[k,] %*% vatilde[ i*p:(i+1)*p-1, j*p:(j+1)*p-1 ] %*% trans[k,]'
 * and p = ncol(trans) (number of lines), vatilde is arranged in trait-major
 * block order.
 *
 * Original code: Russell Edson, Biometry Hub (19/11/2024).
 * Integrated into wgAim package: Session 34.
 */
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute vqtilde — denominator of the MV outlier statistic
//'
//' For each marker row of \code{trans}, computes
//' \code{tr(Ginv \%*\% (Ti \%*\% vatilde \%*\% t(Ti)))} where
//' \code{Ti = kronecker(diag(ntrait), trans[k,,drop=FALSE])}.
//' This is the denominator of the multivariate outlier statistic
//' \eqn{t^2_{kl}} from Verbyla & Cullis (2012), equation (11).
//'
//' @param trans   nmarkers x nlines back-transform matrix M^T(MM^T)^{-1}.
//' @param Ginv    ntrait x ntrait generalised inverse of Ga.
//' @param vatilde (ntrait*nlines) x (ntrait*nlines) posterior variance
//'   (Ga x relm - PEV), arranged in trait-major block order.
//' @param ntrait  Number of traits/environments.
//' @return Numeric vector of length nmarkers.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector compute_vqtilde(Rcpp::NumericMatrix trans,
                                     Rcpp::NumericMatrix Ginv,
                                     Rcpp::NumericMatrix vatilde,
                                     int ntrait) {
    arma::uword n = trans.nrow(), p = trans.ncol();
    arma::mat trans_mat(trans.begin(),   n, p,               false);
    arma::mat Ginv_mat (Ginv.begin(),    Ginv.nrow(),  Ginv.ncol(),  false);
    arma::mat vat_mat  (vatilde.begin(), vatilde.nrow(), vatilde.ncol(), false);

    Rcpp::NumericVector vqtilde(n);
    Rcpp::List dimnames = trans.attr("dimnames");
    if (!Rf_isNull(dimnames)) {
        Rcpp::CharacterVector rn = Rcpp::as<Rcpp::List>(dimnames)[0];
        vqtilde.attr("names") = rn;
    }

    arma::mat varq(ntrait, ntrait);
    for (arma::uword k = 0; k < n; k++) {
        for (arma::uword i = 0; i < (arma::uword)ntrait; i++) {
            for (arma::uword j = 0; j < (arma::uword)ntrait; j++) {
                varq(i, j) = arma::as_scalar(
                    trans_mat.row(k) *
                    vat_mat.submat(i * p, j * p, (i + 1) * p - 1, (j + 1) * p - 1) *
                    trans_mat.row(k).t()
                );
            }
        }
        vqtilde[k] = arma::trace(Ginv_mat * varq);
    }
    return vqtilde;
}
