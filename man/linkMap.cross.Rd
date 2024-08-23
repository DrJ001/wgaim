\name{linkMap.cross}
\alias{linkMap.cross}
\alias{linkMap}
\title{Plot a genetic linkage map}
\description{
Neatly plots the genetic linkage map with marker locations and marker names.
}
\usage{
\method{linkMap}{cross}(object, chr, chr.dist, marker.names = "markers",
     tick = FALSE, squash = TRUE, m.cex = 0.6, \ldots)
}
\arguments{
  \item{object}{object of class \code{"cross"}}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{chr.dist}{a list containing named elements \code{"start"} and
    \code{"end"} containing the start and end distances in cM the
    genetic map should be subsetted by. Each of these may also be
    a vector of distances equal to the length of the number of linkage
    groups to be plotted.}
  \item{marker.names}{a character string naming the type of marker
    information to plot. If \code{"dist"} then distances are
    plotted alongside each chromosome on the left. If \code{"markers"}
    then marker names are plotted instead. Defaults to \code{"markers"}}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names. Defaults to \code{FALSE}}
  \item{squash}{logical value. if \code{TRUE} then creates extra room on
    the left side of the chromosomes. This is useful for plotting trait
    names for QTL using \code{linkMap.wgaim()} and
    \code{linkMap.default()}}
  \item{m.cex}{the expansion factor to use for the marker names}
  \item{\ldots}{arguments passed to \code{plot()} function to set up the plot
    region. Arguments may also be passed to \code{text()} function for the
    manipulation of the marker names}
}
\details{
This plotting procedure provides a visual display of the
chromosomes without marker names overlapping vertically. The plotting
region will adjust itself to ensure that all marker names are in the region. For
this reason the value for \code{"m.cex"} is passed to the \code{text()}
function and should be manipulated until an aesthetic genetic map is reached.

For large maps with many chromosomes, marker names and adjacent
chromosomes will overlap horizontally. For the interest of readability
this has not been corrected. For this particular situation it is
suggested that the user horizontally maximise the plotting window
until no overlapping occurs or subset the genetic map to achieve the
desired result.
}
\value{This invisibly returns the following list for manipulation with
  \code{linkMap.wgaim()}
  \item{mt}{A list named by the chromosomes with each element containing the locations of
   the marker names after correcting for overlapping}
  \item{map}{A list named by the chromosomes with each element
    containing the locations of markers on the chromosomes}
  \item{chrpos}{The numerical position of the chromosomes on the
    plotting region}
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
\seealso{\code{\link{linkMap.wgaim}}}
\examples{

data(genoSxT, package = "wgaim")

## plot linkage map with marker names

linkMap(genoSxT, cex = 0.5)

## plot linkage map with distances

linkMap(genoSxT, cex = 0.5, marker.names = "dist")

}
\keyword{hplot}
