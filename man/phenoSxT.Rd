\name{phenoSxT}
\alias{phenoSxT}
\docType{data}
\title{Phenotypic Sunco x Tasman trial data frame}
\description{Phenotype data arising from a two-phase experiment
  involving a Doubled Haploid population from a crossing of the wheat varities
  Sunco and Tasman}
\usage{data(phenoSxT)}
\format{This data relates to a two-phase epxeriment involving a set of 175
  Doubled Haploid lines. In the first phase DH lines were randomly
  allocated to plots using a complete block design with additional plots
  containing the parents (Sunco, Tasman) as well as commercial lines
  (Frame, Janz, Krichauff, Machete, RAC820, Trident). The trial was laid
  out in a rectangular array of 31 rows and 12 columns. In the second
  phase 23\% of the field samples were replicated in the milling
  process producing a total of 456 milling samples. These partially
  replicated field samples were then randomly allocated to 38
  mill days with 12 samples per mill day. The data frame
consists of 456 rows with 11 columns. These columns are
  \describe{
    \item{Expt:}{A one level of factor containing a unique
      identification for the experiment.}
    \item{Type:}{A factor of nine levels indexing the wheat variety
      (Doubled Haploid, Sunco, Tasman, (Frame, Janz, Krichauff, Machete, RAC820, Trident))}
    \item{id:}{A factor of 183 levels uniquely identifying the wheat
      varieties involved in the experiment.}
    \item{Range:}{A factor of 12 numeric levels indexing the field Range.}
    \item{Row:}{A factor of 31 numeric levels indexing the field Rows.}
    \item{Rep:}{A factor of 2 levels indexing the Block of the
      experiment}
    \item{Millday:}{A factor of 38 numeric levels indexing the milling
  day}
    \item{Millord:}{A factor of 12 levels indexing the milling order}
    \item{myield:}{A numeric vector of milling yield observations from
      the second phase of the experiment.}
    \item{lord:}{A centered numerical vector of milling orders, Millord}
    \item{lrow:}{A centered numerical vector of Rows}
    }
  }
\examples{
data(phenoSxT, package = "wgaim")

}
\keyword{datasets}
