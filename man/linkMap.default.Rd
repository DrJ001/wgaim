\name{linkMap.default}
\alias{linkMap.default}
\title{Plot a genetic linkage map with QTL for multiple traits}
\description{
Neatly plots the genetic linkage map with marker locations,
marker names and highlights QTL with their associated flanking markers
for multiple traits obtained from a list of \code{wgaim} models.
}
\usage{
\method{linkMap}{default}(object, intervalObj, chr, chr.dist, marker.names
    = "markers", flanking = TRUE, list.col = list(q.col = rainbow(length(object)),
    m.col = "red", t.col = rainbow(length(object))), list.cex =
    list(m.cex = 0.6, t.cex = 0.6), trait.labels = NULL, tick = FALSE, \ldots)
}
\arguments{
  \item{object}{a list object with elements inheriting the class \code{"wgaim"}}
  \item{intervalObj}{object of class \code{"cross"} or \code{"interval"}}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{chr.dist}{a list containing named elements \code{"start"} and
    \code{"end"} containing the start and end distances in cM the
    genetic map should be subsetted by. Each of these may also be
    a vector of distances equal to the length of the number of linkage
    groups to be plotted.}
  \item{marker.names}{a character string naming the type of marker
    information to plot. If \code{"dist"} then distances are
    plotted alongside each chromosome on the left. If \code{"markers"}
    then marker names are plotted instead. Defaults to
    \code{"markers"}.}
  \item{flanking}{logical value. If \code{TRUE} then only plot marker
     names or distances for flanking markers of the QTL. Defaults to
     \code{TRUE}}
  \item{list.col}{named list of colors used to highlight the QTL regions and
    their flanking markers. \code{q.col} is the colors of the QTL
    regions (defaults to \code{rainbow(n)} where \code{n} is the length
    of \code{object}). \code{m.col} is the color the flanking markers. \code{t.col}
    is the color of the trait names used in each model (defaults to the
    same color as the QTL regions). See \code{par} for color options}
  \item{list.cex}{a named list object containing the character expansion
    factors for the marker names \code{m.cex} and the trait labels \code{t.cex}}
  \item{trait.labels}{character string naming the trait used in the
    model object, defaults to the response names extracted from the
    fixed compoenent of the base model.}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names}
  \item{\ldots}{arguments passed to the \code{plot()} function to set up the plot
    region. Arguments may also be passed to the \code{text()} function for the
    manipulation of the marker names}
}
\details{
This plotting procedure is a wrapper for \code{linkMap.wgaim} and displays
QTL for multiple traits obtained from a list of models given by \code{object}.
Alternative labels for the traits can be given, in model order, using
\code{trait.labels}.

Color specific highligting of the QTL is also available using
\code{clist}. This differs slightly from \code{linkMap.wgaim}. Here
the \code{q.col} and \code{t.col} should be given a set of colors equal to
the length of \code{object}. Let \code{n} be the length of \code{object}.
Then if \code{q.col} is \code{NULL} or length of \code{q.col} is not equal
to \code{n} then it defaults to \code{rainbow(n)}. If \cite{t.col} is
\code{NULL} or length of \code{t.col} is not equal to \code{n} or 1 then
it defaults to the colors of \code{q.col}. Examples of different color
combinations are given below in the examples.

The \code{list.cex} argument can be used to manipulate the character expansion of
the marker names using \code{m.cex} or the character expansion of the
\code{trait.labels} using \code{t.cex}. If a set of \code{"marker"} analyses has been
performed then \code{pch} is used to plot a symbol at the
location of the QTL. This character can be changed using the usual
arguments such as \code{pch} or \code{cex} that are passed through the
\ldots argument.
}
\value{For a set of \code{"interval"} analyses, the genetic linkage map is
  plotted with shaded QTL regions and highlighted flanking markers. For
  a set of \code{"marker"} analyses, symbols are
  placed at the QTL locations and the markers are highlighted.
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
\author{Julian Taylor}
\seealso{\code{\link{linkMap.cross}}, \code{\link{linkMap.wgaim}}}
\examples{

\dontrun{
## fit wgaim models

rktgw.qtl <- wgaim(rktgw.asf, intervalObj = genoRxK, merge.by = "Genotype",
                 trace = "trace.txt", na.action = na.method(x = "include"))

rkyld.qtl <- wgaim(rkyld.asf, intervalObj = genoRxK, merge.by = "Genotype",
                 trace = "trace.txt", na.action = na.method(x = "include"))

## plot QTL intervals

# matching rainbow QTL color and trait names, red flanking markers
# (default) and gray background markers.

linkMap(list(rktgw.qtl,rkyld.qtl), genoRxK, col = "gray")

# rainbow QTL color and black trait names, red flanking markers
# (default) and gray background markers.

linkMap(list(rktgw.qtl,rkyld.qtl), genoRxK, list.col = list(t.col =
"black", m.col = "red"), col = "gray")

# monochromatic plot: gray QTLs, black trait names, black flanking
# markers and gray background markers

linkMap(list(rktgw.qtl,rkyld.qtl), genoRxK, list.col = list(q.col =
rep(gray(0.8), 2), t.col = "black", mcol = "black"), col = "gray")

}
}
\keyword{hplot}