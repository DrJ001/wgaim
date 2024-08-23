\name{cross2int}
\alias{cross2int}
\title{Convert a cross genetic object to an interval object}
\description{
Converts an object of class "\code{cross}" to an object with class
"\code{interval}". The function also imputes missing markers.
}
\usage{
cross2int(object, impute = "MartinezCurnow", consensus.mark = TRUE,
     id = "id", subset = NULL)
}
\arguments{
  \item{object}{an object of class \code{"cross"} that inherits one of the class
    structures \code{"bc"}, \code{"dh"}, \code{"f2"}, \code{"riself"}.}
  \item{impute}{a character string determining how missing values in
    the linkage map should be imputed. If \code{"Broman"}, then missing values are
    imputed according to Bromans rules. If \code{"MartinezCurnow"} then
    missing values are imputed according to the rules of Martinez &
    Curnow (1994) (see reference list). The default is \code{"MartinezCurnow"} (see Details).}
  \item{consensus.mark}{logical value. If \code{TRUE} co-locating marker sets
    are condensed to form consensus markers (see Details). Defaults to
  \code{TRUE}.}
  \item{id}{a character string or name of the unique identifier for each row of genotype
    data (see Details). Defaults to \code{"id"}}
  \item{subset}{a possible character vector naming the subset of
    chromosomes to be returned. Defaults to \code{NULL} implying return
    all chromosomes.}
}
\details{
This function provides the conversion of genetic data objects that have
already been generated using \code{read.cross()} from Bromans \pkg{qtl}
package, to \code{"interval"} objects ready for use with
\code{wgaim}. Users should be aware that this function is restricted to
certain populations.  \code{object} must inherit one of the class
structures \code{"bc"}, \code{"dh"}, \code{"f2"}, \code{"riself"}.

During the conversion process three important linkage map attributes are
assessed.
\enumerate{
  \item The map may be subsetted using the \code{subset} argument
  \item If \code{consensus.mark = TRUE} then co-located marker sets are reduced
  to form single consensus markers before missing values are
  imputed. The marker similarity is determined by
  the genetic distances that are given in the map component for each linkage
  group. If a set of markers co-locate the name of the first marker is
  chosen and a single consensus marker is determined by coalescing the
  genetic information from all markers in the set. A "(C)" is placed
  after the marker name for easy identification. The markers removed
  from each set are returned with the object and placed under
  \code{"colocated.markers"} for inspection if required.
  \item Missing values are imputed according to
  the argument given by \code{impute}. This imputation results in a
  complete version of the marker data for each chromosome which is then
  used to create the interval data component "\code{interval.data}". The complete
  marker data for each chromosome can be obtained from the "\code{imputed.data}" element of the
  returned list. It is therefore also possible to perform whole genome \emph{marker}
  analysis using \code{wgaim}. See \code{wgaim.asreml} for more details.
}
}
\value{
a list of class \code{"cross"} that also inherits the class
\code{"interval"}. The list contains the following components
\item{geno}{
  A list with elements named by the corresponding names of the
  chromosomes. Each chromosome is itself a list with six
  elements: \code{"data"} is the actual estimated map matrix with rows
  as individuals named by \code{"id"} and markers as columns;
  \code{"map"} is a vector of marker positions on the
  corresponding chromosome; \code{"imputed.data"} is identical to \code{"data"}
  matrix but with all NAs replaced by imputed values according to the
  rules of \code{"impute"}; \code{"dist"} contains the genetic
  distance between adjacent markers or the genetic distances of the
  intervals; \code{"theta"} contains the recombination fractions for
  each interval; \code{"interval.data"} contains the recalculated intervals
  based on the recombination fractions and the missing marker information.}
\item{colocated.markers}{
If \code{consensus.mark = TRUE}, a four column data frame containing
stacked binned sets of co-located markers. In each binned set the first
marker is the unique consensus marker name used in the linkage map and the
others are the co-located marker names that were omitted. Additionally
for each binned set, the data frame also contains linkage group
identification and marker position information.}
\item{pheno}{
  A data.frame of phenotypic information with rows as individuals read
  in from \code{read.cross}. A copy of the column named by the \code{"id"} argument
  can be found here (see the help for the \code{read.cross() function}).}
}
\references{
  Martinez, O., Curnow. R. N. (1994) Missing markers when estimating
  quantitative trait loci using regression mapping. \emph{Heredity},
  \bold{73}, 198-206.

  Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
  Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
  Statistical Software}, \bold{40}(7), 1-18. URL \url{https://www.jstatsoft.org/v40/i07/}.

  Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
  by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 95-111.
}
\author{Julian Taylor and Ari Verblya}
\seealso{\code{\link[qtl]{read.cross}}}
\examples{

\dontrun{
# read in linkage map from a rotated .CSV file with "id" as the
# identifier for each unique row

wgpath <- system.file("extdata", package = "wgaim")
genoSxT <- read.cross("csvr", file="genoSxT.csv", genotypes=c("AA","BB"),
         na.strings = c("-", "NA"), dir = wgpath)
genoSxT <- cross2int(genoSxT, impute="MartinezCurnow", id = "id")

# plot linkage map

linkMap(genoSxT, cex = 0.5)

}
}
\keyword{regression}
