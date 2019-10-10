\name{genoRxK}
\alias{genoRxK}
\docType{data}
\title{Genotypic marker data for RAC875 x Kukri doubled haploid
  population in R/qtl format}
\description{Linkage map marker data for the RAC875 x Kukri doubled haploid
  population in the form of an R/qtl cross object.}
\usage{data(phenoSxT)}
\format{
  This data relates to a linkage map of 500 genetic markers genotyped on 368
  individuals from the RAC875 x Kukri population. The linkage map
  consists of 21 linkage groups with varying numbers of markers. Map
  distances have been estimated using \code{read.cross()} with the kosambi
 mapping function. The data is therefore an R/qtl
 cross object. See \code{read.cross()} documentation for more details on the
 format of this object.
  }
\examples{

data(genoRxK, package = "wgaim")
linkMap(genoRxK, cex = 0.5)

}
\keyword{datasets}
