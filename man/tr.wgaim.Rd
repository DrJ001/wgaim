\name{tr.wgaim}
\alias{tr}
\alias{tr.wgaim}
\title{Display diagnostic information about the QTL detected.}
\description{Displays diagnostic infomation about QTL detection and
  significance for the sequence of models generated in a \code{wgaim}
  analysis.
}
\usage{
\method{tr}{wgaim}(object, iter = 1:length(object$QTL$effects),
      lik.out = TRUE, \ldots)
}
\arguments{
  \item{object}{an object of class \code{"wgaim"}}
  \item{iter}{a vector of integers specifying what rows of the p-value matrix to display}
  \item{lik.out}{logical value. If \code{TRUE} then diagnostic information
about the testing of the genetic variance is given for all iterations.}
  \item{\ldots}{arguments passed to \code{print.default} for displaying
of information}
}
\details{
By default the printing of the objects occur with arguments \code{quote = FALSE}
and \code{right = TRUE}. Users should avoid using these arguments.
}
\value{For the selected QTL, a probability value matrix is displayed
with rows specified by \code{iter}. If \code{lik.out =
TRUE} then a matrix with rows consisting of the likelihood with
additive genetic variance, the likelihood without additive genetic
variance (NULL model), the test statistic and the p-value for the statistic.
 }
 \references{
Verbyla, A. P & Taylor, J. D, Verbyla, K. L (2012). RWGAIM: An efficient high
dimensional random whole genome average (QTL) interval mapping approach.
\emph{Genetics Research}. \bold{94}, 291-306.

Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{wgaim}}}
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

# diagnostic check

tr(rkyld.qtl, digits = 4)
}
}

\keyword{regression}

