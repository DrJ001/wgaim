\name{phenoCxR}
\alias{phenoCxR}
\docType{data}
\title{Phenotypic Cascades x RAC875-2 zinc experiment data frame}
\description{Zinc concentration data of a Doubled Haploid wheat population.
}
\usage{data(phenoCxR)}
\format{This data relates to a glasshouse experiment involving a set of 90
  Doubled Haploid (DH) lines from a crossing of Cascades x Rac875-2. The
  DH lines were allocated randomly to pots in the glasshouse using a
  randomised complete block design. There were also additional pots
  that contained 5 of each of the parents (Cascades and Rac-875-2). Two
  measurments were made, namely zinc concentration and shoot length.
  The data frame consists of 200 rows and 5 columns decribed below
  \describe{
    \item{id:}{A factor of 92 levels conataining the unique
      identification of the DH lines and parents.}
    \item{Block:}{A factor of two levels indexing the blocks in the experiment}
    \item{Type:}{A factor of 3 levels indexing the wheat variety (Doubled
      Haploid, Cascades, Rac875-2)}
    \item{shoot:}{A numerical variable of shoot lengths for each plant}
    \item{znconc:}{A numerical variable of zinc concentration levels for each plant}
    }
  }
\examples{
data(phenoCxR, package = "wgaim")
plot(phenoCxR$shoot, phenoCxR$znconc)

}
\keyword{datasets}
