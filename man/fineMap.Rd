\name{fineMap}
\alias{fineMap}
\title{Fine mapping of a QTL region using a sliding window approach}
\description{
  Performs fine mapping of a previously detected QTL by scanning a window of
  markers or intervals flanking a nominated marker. At each position in the
  window, the nominated QTL is replaced in the model and a likelihood ratio
  test statistic and LOD score are computed, allowing the identification of the
  most likely QTL position within the region.
}
\usage{
fineMap(model, intervalObj, mark = NULL, flanking = 50,
        exclusion.window = 50, \ldots)
}
\arguments{
  \item{model}{a fitted model object of class \code{"wgaim"} obtained from a
    call to \code{wgaim()} (see Details).
  }
  \item{intervalObj}{a list object of class \code{"cross"} or
    \code{"interval"} containing the genotypic data, usually obtained
    from a call to \code{cross2int()} (see Details).
  }
  \item{mark}{a character string giving the name of the QTL marker or interval
    midpoint to be fine mapped. This marker must exist in
    \code{intervalObj}. The default is \code{NULL} and a non-\code{NULL}
    value must be supplied.
  }
  \item{flanking}{a numerical value giving the number of markers or intervals
    to the left and right of \code{mark} that define the fine mapping
    window. Defaults to \code{50}.
  }
  \item{exclusion.window}{for each position tested in the window, the
    distance in centimorgans on the left and right side of that position
    that is excluded from the remaining marker set during model fitting.
    Defaults to \code{50}.
  }
  \item{\ldots}{any extra arguments passed to the \code{asreml()} calls
    (see \code{?asreml} and \code{?asreml.options} for more information).
  }
}
\details{
  \code{fineMap()} is designed to be used after a successful \code{wgaim()}
  analysis to refine the location of a single detected QTL. The function
  requires a fitted \code{"wgaim"} model object and the corresponding
  \code{intervalObj} used in the original analysis.

  The marker name supplied to \code{mark} is used to identify the chromosome
  and position of the QTL to be fine mapped. If the original \code{wgaim()}
  analysis used \code{gen.type = "interval"} then \code{mark} should be an
  interval midpoint name from the \code{inferred.map} component of the relevant
  chromosome in \code{intervalObj}. If \code{gen.type = "markers"} was used
  then \code{mark} should be a marker name from the \code{map} component.

  A sliding window of width \code{2 * flanking + 1} positions centred on
  \code{mark} is constructed. At each position in the window, the nominated
  QTL term is replaced in the model by the current window position and the
  model is refitted using \code{update.asreml()}. The method of model fitting
  (i.e. \code{method = "fixed"} or \code{method = "random"}) is inherited
  directly from the \code{wgaim()} model supplied to \code{model}.

  At each window position a Wald-type likelihood ratio test statistic is
  computed from the coefficient and its variance, and a LOD score is derived
  from the test statistic. The p-value is calculated using a chi-squared
  distribution with 1 degree of freedom.

  The phenotypic data used for model updating is recovered from the name of the
  response variable in the fixed formula of \code{model} by looking for an
  object named \code{<response>.data} in the parent environment. This object
  must exist and be accessible at the time \code{fineMap()} is called.

  When the number of active markers in the reduced marker set exceeds the
  number of lines, a genomic relationship matrix is constructed via an internal
  call to \code{constructCM()} and passed to the \code{vm()} term in the
  \code{asreml()} model. Otherwise the reduced marker matrix is used directly.
}
\value{
  A \code{data.frame} with one row per position in the fine mapping window,
  containing the following columns:
  \item{mark}{the name of the marker or interval midpoint at each window
    position.}
  \item{dist}{the genetic map position (in cM) of each window position.}
  \item{pvalue}{the p-value from the chi-squared likelihood ratio test at
    each window position.}
  \item{LOD}{the LOD score at each window position.}
}
\references{
Verbyla, A. P & Taylor, J. D, Verbyla, K. L (2012). RWGAIM: An efficient high
dimensional random whole genome average (QTL) interval mapping approach.
\emph{Genetics Research}. \bold{94}, 291-306.

Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
Statistical Software}, \bold{40}(7), 1-18. URL \url{https://www.jstatsoft.org/v40/i07/}.

Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
by simultaneous use of the full linkage map. \emph{Theoretical and
Applied Genetics}, \bold{116}, 95-111.
}
\author{Julian Taylor and Ari Verbyla}
\seealso{\code{\link{wgaim.asreml}}, \code{\link{cross2int}}, \code{\link{summary.wgaim}}}
\examples{
\dontrun{
data(phenoRxK, package = "wgaim")
data(genoRxK, package = "wgaim")

## convert cross object to interval object
genoRxKi <- cross2int(genoRxK, consensus.mark = TRUE,
                      impute = "MartinezCurnow", id = "Genotype")

## fit base asreml model
rkyld.asf <- asreml(yld ~ Type + lrow, random = ~ Genotype + Range,
                    residual = ~ ar1(Range):ar1(Row), data = phenoRxK)

## run wgaim QTL analysis
rkyld.qtl <- wgaim(rkyld.asf, intervalObj = genoRxKi, merge.by = "Genotype",
                   gen.type = "interval", method = "random",
                   trace = "trace.txt", na.action = na.method(x = "include"))

## summarise detected QTL to identify a marker name for fine mapping
summary(rkyld.qtl, genoRxKi)

## fine map around a detected QTL interval midpoint
rkyld.fm <- fineMap(rkyld.qtl, genoRxKi, mark = "3B.139.2",
                    flanking = 50, exclusion.window = 50,
                    na.action = na.method(x = "include"))

## inspect fine mapping results
print(rkyld.fm)
}
}
\keyword{regression}
