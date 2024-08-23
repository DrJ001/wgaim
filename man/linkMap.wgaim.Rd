\name{linkMap.wgaim}
\alias{linkMap.wgaim}
\title{Plot a genetic linkage map with QTL}
\description{
Neatly plots the genetic linkage map with marker locations,
marker names and highlights QTL with their associated flanking markers
obtained from a \code{wgaim} model.
}
\usage{
\method{linkMap}{wgaim}(object, intervalObj, chr, chr.dist,
    marker.names = "markers", flanking = TRUE, list.col = list(q.col = "light blue",
    m.col = "red", t.col = "light blue"), list.cex = list(t.cex = 0.6,
    m.cex = 0.6), trait.labels = NULL, tick = FALSE, \ldots)
}
\arguments{
  \item{object}{object of class \code{"wgaim"}}
  \item{intervalObj}{object of class \code{"cross"} or \code{"interval"}}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{chr.dist}{a list containing named elements \code{"start"} and
    \code{"end"} containing the start and end distances in cM the
    genetic map should be subsetted by. Each of these may also be
    a vector of distances equal to the length of the number of linkage
    groups to be plotted.}
   \item{marker.names}{a character string naming the type of marker
    information to plot. If \code{"dist"} then distances names
    plotted alongside each chromosome on the left. If \code{"markers"}
    then marker names are plotted instead. Defaults to
    \code{"markers"}.}
  \item{flanking}{logical value. If \code{TRUE} then only plot marker
     names or distances for flanking markers of the QTL. Defualts to \code{TRUE}}
  \item{list.col}{named list of colours used to highlight the QTL regions and
    their flanking markers. \code{q.col} is the color of the QTL
    regions. \code{m.col} is the color the flanking markers. \code{t.col}
    is the color of the trait name used in the model object (see
    \code{par} for colour options)}
  \item{list.cex}{a named list object containing the character expansion
    factors for the marker names \code{m.cex} and the trait labels \code{t.cex}}
  \item{trait.labels}{character string naming the trait used in the
    model object}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names}
  \item{\ldots}{arguments passed to the \code{plot()} function to set up the plot
    region and plot any symbols if required}
}
\details{
This plotting procedure builds on \code{linkMap.cross()} by adding the
QTL regions to the map and highlighting the appropriate markers obtained
from a fit to \code{wgaim}. If the linkage map is subsetted and QTL
regions fall outside the remaining map a warning will be given that
the QTL have been omitted from the display.

The \code{list.col} arguments \code{q.col}, \code{m.col} and
\code{t.col} have been added for personal colour highlighting of the QTL
regions, flanking markers and trait names. For greater flexibility the
procedure may also be given the usual \code{col} argument that will be
passed to the other markers.

The \code{list.cex} argument can be used to manipulate the character expansion of
the marker names using \code{m.cex} or the character expansion of the
\code{trait.labels} using \code{t.cex}. If a \code{"marker"} analysis has been
performed then \code{pch} is used to plot a symbol at the
location of the QTL. This character can be changed using the usual
arguments such as \code{pch} or \code{cex} that are passed through the
\ldots argument.
}
\value{For an \code{"interval"} analysis, the genetic linkage map is
  plotted with shaded QTL regions and highlighted flanking markers. For
  a \code{"marker"} analysis, a symbol is placed at the QTL locations and
  the markers are highlighted.
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
\author{Julian Taylor}
\seealso{\code{\link{linkMap.cross}}, \code{\link{wgaim}}}
\examples{

\dontrun{
# fit wgaim model

yield.qtl <- wgaim(yield.fm, intervalObj = genoRxK, merge.by = "Genotype",
                  trace = "trace.txt", na.action = na.method(x = "include"))

# plot QTL

linkMap(yield.qtl, genoRxK, list.col = list(m.col = "red"), col = "gray")

}
}
\keyword{hplot}
