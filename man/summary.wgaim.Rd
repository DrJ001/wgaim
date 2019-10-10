\name{summary.wgaim}
\alias{summary.wgaim}
\alias{print.wgaim}
\alias{print.summary.wgaim}
\title{Summary and print methods for the class \code{"wgaim"}}
\description{
Prints a QTL summary form the \code{"wgaim"} object in a presentable format
}
\usage{
\method{summary}{wgaim}(object, intervalObj, LOD = TRUE, \ldots)
\method{print}{wgaim}(x, intervalObj, \ldots)
}
\arguments{
  \item{object}{an object of class \code{"wgaim"} (see Details)}
  \item{x}{an object of class \code{"wgaim"}}
  \item{intervalObj}{a data structure of class \code{"cross"} or
    \code{"interval"} containing the genotypic data}
  \item{LOD}{logical value. If TRUE LOD scores for QTL are calculated, defaults to \code{TRUE}}
  \item{\ldots}{further arguments passed to or from other methods}
}
\details{
  It is important that the \code{intervalObj} is not missing in
  \code{summary.wgaim()} or \code{print.wgaim()} as it
  contains vital summary information about each of the QTL
  detected.

  The summary of the QTL differs depending on the method chosen
  in the \code{wgaim.asreml} call. If \code{method = "random"}
  then the significance of the QTL are summarized using a probablistic
  argument based on the conditional distribution of the QTL sizes given
  the data (see Verbyla et. al, 2012 in References) Thus, for each
  QTL, a value is calculated that represents the probability that the
  QTL size is greater than zero (or less than zero if the effect is
  negative). If \code{method = "fixed"} then the significance of the QTL is
  summarized using a one degree of freedom Wald statistic.
}
\value{A summary of the QTL component of the \code{"wgaim"} object is
  printed to the screen. For each QTL detected, if an \code{"interval"}
  analysis was performed then \code{summary.wgaim()} prints which
  chromosome, name and distance of each flanking marker, size,
  probability/p-value, contribution of genetic variance and LOD
  score if desired. If a \code{"marker"} analysis was performed then the
  chromosome, name and distance of the associated marker, size,
  probability/p-value, contribution of genetic variance
  and LOD score are printed. \code{print.wgaim()} provides a narrative
  brief of the QTL detected.
}
\references{
  Verbyla, A. P & Taylor, J. D, Verbyla, K. L (2012). RWGAIM: An efficient high
  dimensional random whole genome average (QTL) interval mapping approach.
  \emph{Genetics Research}. \bold{94}, 291-306.

  Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
  Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
  Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.

  Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
  by simultaneous use of the full linkage map. \emph{Theoretical and
    Applied Genetics}, \bold{116}, 195-211.
  }
\author{Julian Taylor and Ari Verbyla}
\seealso{\code{\link{wgaim.asreml}}}
\examples{

\dontrun{
# read in data

data(phenoRxK, package = "wgaim")
data(genoRxK, package = "wgaim")

# subset linkage map and convert to "interval" object

genoRxK <- subset(genoRxK, chr = c("1A", "2D1", "2D2", "3B"))
genoRxK <- cross2int(genoRxK, impute = "Martinez", id = "Genotype")

# base model

rkyld.asf <- asreml(yld ~ lrow, random = ~ Genotype + Range,
                   residual = ~ ar1(Range):ar1(Row), data = phenoRxK)

# find QTL

rkyld.qtl <- wgaim(rkyld.asf, intervalObj = genoRxK, merge.by = "Genotype",
                   trace = "trace.txt", na.action = na.method(x = "include"))

# summarise

print(rkyld.qtl, genoRxK)
summary(rkyld.qtl, genoRxK)

}
}
\keyword{regression}

