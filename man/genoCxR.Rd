\name{genoCxR}
\alias{genoCxR}
\docType{data}
\title{Genotypic marker data for Cascades x RAC875-2 doubled haploid
  population in R/qtl format}
\description{Linkage map marker data for the Cascades x RAC875-2 doubled
  haploid population in the form of an R/qtl cross object.
}
\usage{data(genoCxR)}
\format{This data relates to a linkage map of 663 markers genotyped on 93
  individuals. The linkage map consists of 42 linkage groups spanning
  the whole genome. Coincident markers have been removed reducing the
  linkage map to 458 markers. Map distances have been
  estimated using \code{read.cross()} with the haldane
  mapping function. The data object is therefore an R/qtl
  cross object. See \code{read.cross()} documentation for more details on the
  format of this object.
}
\examples{

data(genoCxR, package = "wgaim")
linkMap(genoCxR, cex = 0.5)

}
\keyword{datasets}
