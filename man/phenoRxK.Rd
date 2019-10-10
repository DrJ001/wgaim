\name{phenoRxK}
\alias{phenoRxK}
\docType{data}
\title{Phenotypic RAC875 x Kukri trial data frame}
\description{Phenotype data arising from a field trial of a
  Doubled Haploid population involving a crossing of the wheat varities
  RAC875 and Kukri}
\usage{data(phenoRxK)}
\format{This data relates to a field trial conducted in 2007 at the
  Roseworthy Cmapus of the University of Adelaide. The trial consisted
  of 2 replicates of 254 Doubled Haploid lines from a cross between wheat
  varieties RAC875 and Kukri. The DH lines, the parents
  (RAC875, Kukri) and control varieties (ATIL, SOKOLL, WEEBIL) were
  randomly allocated to 520 plots using a randomized complete bloack
  design. The trial was laid out in a 20 by 26 rectangular array. The
  data frame consists of 520 rows with 9 columns described by:
  \describe{
    \item{Genotype:}{A 254 level factor containing a unique
      identification for the wheat varieties involved in the experiment.}
    \item{Type:}{A factor of four levels indexing the wheat varieties
      (Doubled Haploid, RAC875, Kukri, ATIL, SOKOLL, WEEBIL).}
    \item{Range:}{A factor of 20 numeric levels indexing the field Range.}
    \item{Row:}{A factor of 26 numeric levels indexing the field Rows.}
    \item{Block:}{A factor of 2 levels indexing the Blocks of the
      experiment.}
    \item{yield:}{A numeric vector of yield observations in kg/ha.}
    \item{tgw:}{A numeric vector of thousand grain weight observations.}
    \item{lrange:}{A centred numerical vector of the field Ranges.}
    \item{lrow:}{A centred numerical vector of the field Rows.}
     }
  }
\examples{
data(phenoRxK, package = "wgaim")

}
\keyword{datasets}
