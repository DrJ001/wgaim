\name{wgaim.asreml}
\alias{wgaim.asreml}
\alias{wgaim}
\title{wgaim method for class "asreml"}
\description{Implements the iterative Whole Genome Average Interval Mapping
  (wgaim) algorithm using the functionality of the flexible linear mixed
  modelling R package \pkg{ASReml-R} V4.
}
\usage{
\method{wgaim}{asreml}(baseModel, intervalObj, merge.by = NULL,
         fix.lines = TRUE, gen.type = "interval", method = "fixed",
         selection = "interval", force = FALSE, exclusion.window = 20,
         breakout = -1, TypeI = 0.05, trace = TRUE, verboseLev = 0, \ldots)
}
\arguments{
\item{baseModel}{a linear mixed model object of class \code{"asreml"} usually representing
    a base \code{asreml()} model to be extended.
}
\item{intervalObj}{a list object containing the genotypic data, usually an \code{"interval"}
    object obtained from using \code{cross2int} (see Details).
}
 \item{merge.by}{a character string or name of the column(s) in \code{phenoData} and
   \code{intervalObj} to merge the phenotypic and genotypic data sets.
 }
\item{fix.lines}{a logical value. If \code{TRUE} then lines existing in
  the phenotype data that do not exist in \code{intervalObj} are fixed and placed in
  the \code{fixed} component of the \code{asreml()} models (see
  Details). It is recommended to set this to TRUE. Defaults to \code{TRUE}.
}
\item{gen.type}{a character string determining the type of genetic data to
  be used in the analysis. Possibilities are \code{"interval"} and
  \code{"markers"}. The default is \code{"interval"}. (see Details).
}
\item{method}{a character string determining the type of algorithm to
  be used in the analysis. Possibilities are \code{"random"} and
  \code{"fixed"}. The default is \code{"random"}. (see Details).
}
\item{selection}{a character string determining the type of selection
  method that is used to select QTL in the analysis. Possibilities are
  \code{"interval"} and \code{"chromosome"}. The default is
  \code{"interval"}. (see Details).
}
\item{force}{a logical value. If \code{TRUE} then force the algorithm to
  to envoke the low-dimensional model that is normally used when the
  number of markers is less than the number of lines. Defaults to \code{FALSE}.
}
\item{exclusion.window}{For each QTL, the distance in centimorgans on
  the left and right side of each QTL that is excluded from further analysis.
}
\item{breakout}{a numerical integer equivalent to the iteration
 where the algorithm breaks out. The default is
  -1 which ensures the algorithm finds all QTL before halting. (see
  Details)
}
\item{TypeI}{a numerical value determining the familywise alpha level of significance
    for detecting a QTL. The default is 0.05.
  }
\item{trace}{a automatic tracing facility. If \code{trace = TRUE} then
    all \code{asreml} output is piped to the screen during the analysis.
  If \code{trace = "file.txt"}, then output from all asreml models is
  piped to "\code{file.txt}". Both trace machanisms will display a
  message if a QTL is detected.
}
 \item{verboseLev}{Numerical value, either 0 or 1, determining the level of tracing
   outputted during execution of the algorithm A 0 value will produce
   the standard model fitting output from the fitted ASReml models
   involved in the forward selection. A value of 1 will add a table of
   chromosome and interval outlier statistics for each iteration.
}
\item{\ldots}{Any extra arguments to be passed to the \code{asreml}
  calls. (see \code{?asreml} and \code{?asreml.options} for more
  information).
}
}
\details{
In the initial call to \code{wgaim.asreml()}, the marker or interval
information is collated from \code{intervalObj}. If \code{gen.type =
  "interval"} then midpoints of intervals are collated from the
\code{"interval.data"} components of the chromosomes in
\code{intervalObj}. If \code{gen.type = "markers"} then markers are
collated from the \code{"imputed.data"} components of the chromosomes in
\code{intervalObj}.

It is recommended to set \code{fix.lines = TRUE} to ensure additive and
non-additive genetic variances are estimated from lines in the
\code{merge.by} component of the phenotypic data that have genetic
marker data in \code{intervalObj}. Lines in the phenotype
\code{merge.by} factor not existing in \code{intervalObj} will be placed
as a fixed factor (called \code{Gomit}) in the \code{asreml} model.
Note, if there are others factors in the model that have some
potential confounding with \code{Gomit} then \code{asreml} will
indicate this with a simple message \code{'Terms with zero df listed
in attribute 'zerodf' of the wald table'} at the end of its iterative
maximisation. This confounding will have no effect on the outcome and
can be safely ignored. If \code{fix.lines = FALSE} is set then all
available lines in the \code{merge.by} component of the phenotypic data
will be used to estimate the non-additive genetic variance.
In this instance, users also need to be aware that \code{asreml} will
output a large number of warnings due to an inherent mismatch in the levels
of the lines contained in the phenotype data compared to the lines in
\code{intervalObj}.

The \code{method} argument in \code{wgaim.asreml()} allows the user access
to two algorithms. If \code{method = "fixed"} the algorithm
places selected QTL as an additive set of fixed
effects in the model as the forward selection algorithm proceeds. If \code{method = "random"}
places selected QTL in the random part of the model as
an additive set of random effects. This new formulation is outlined in
Verbyla et. al (2012).

The \code{selection} argument determines the type of selection algorithm
for the analysis. If \code{selection = "chromosome"} then outlier
statistics for each chromosome are calculated and the largest
chromosome or linkage group is chosen. The largest marker/interval outlier
statistic in this linkage group is then selected as the putative QTL. If
\code{selection = "interval"}, only marker/interval statistics are calculated
and the largest marker/interval is chosen as the putative QTL.

Note: If a genetic map has a small number of markers on a linkage group
then using \code{selection = "chromosome"} as the selection algorithm
is known to be flawed (see Verbyla et. al, 2012). For this reason it is suggested
that this option only be used when there are a moderate number of
markers on each linkage group.

Users can break out of the algorithm using the \code{breakout}
argument. If a numerical value greater than zero is given, then the forward
selection algorithm breaks at the iteration equal to that value and
returns the collected information to this point. This includes
fixed/random QTL effects, diaganostic components such as interval/marker
BLUPs and outlier statistics as well as the trace components of the
algorithm. It should be noted that the algorithm breaks out before a QTL
has been moved to the fixed/random effects and estimated. Therefore a
positive integer, say \code{n} will not return an estimate of the nth
QTL but it will return the outlier statistics or BLUPs for the nth iteration.

It is recommended that \code{trace = "file.txt"} be used to pipe the
sometimes invasive tracing of \code{asreml} licensing and fitting
numerics for each model to a file. Errors, warnings and messages will
still appear on screen during this process. Note some warnings that
appear may be passed through from an asreml call and are outputted upon
exit. These may be ignored as they are handled during the execution of
the function.
}
\value{An object of class \code{"wgaim"} which also inherits the class
  \code{"asreml"} by default. The object returned is actually an \code{asreml}
  object (see \code{asreml.object}) with the addition of components from
  the QTL detection listed below.

\item{QTL}{A list of components from the significant QTL detected
    including a character vector of the significant QTL along with a
    vector of the QTL effect sizes. There are also a number of diagnostic meausres that
    can be found in \code{diag} that are used in conjunction with
    \code{tr.wgaim} and \code{outStat}.}
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
Applied Genetics}, \bold{116}, 95-111.
}
\author{Julian Taylor and Ari Verbyla}
\seealso{\code{\link{print.wgaim}}, \code{\link{summary.wgaim}}}
\examples{

\dontrun{
# read in data

data(phenoRxK, package = "wgaim")
data(genoRxK, package = "wgaim")

# subset linkage map and convert to "interval" object

genoRxK <- subset(genoRxK, chr = c("1A", "2D1", "2D2", "3B"))
genoRxK <- cross2int(genoRxK, impute = "Martinez", id = "Genotype")

# base model

rkyld.asf <- asreml(yld ~ Type + lrow, random = ~ Genotype + Range,
                   residual = ~ ar1(Range):ar1(Row), data = phenoRxK)

# detect and estimate QTL

rkyld.qtl <- wgaim(rkyld.asf, intervalObj = genoRxK, merge.by = "Genotype",
                   trace = "trace.txt", na.action = na.method(x = "include"))

}
}
\keyword{regression}
