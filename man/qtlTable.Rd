\name{qtlTable}
\alias{qtlTable}
\title{Stack QTL summary information into a super table}
\description{Stack QTL summary information into a super table ready for
  simple exporting.
}
\usage{
qtlTable(\ldots, intervalObj = NULL, labels = NULL, columns = "all")
}
\arguments{
  \item{\ldots}{list of objects of class \code{"wgaim"}. All
    models must have been analysed with the same \code{gen.type}
    (see help for \code{wgaim.asreml()}).}
  \item{intervalObj}{a genetic object of class \code{"interval"}
    reequire in a \code{wgaim} analysis  (see help for
    \code{wgaim.asreml()}). This is required to be non NULL}
  \item{labels}{a vector of character strings describing the trait
    names of each model QTL table.}
  \item{columns}{this can be either a numeric vector determining which columns of the QTL
    summaries should be outputted or \code{"all"} for all
    columns. The default is \code{"all"}.}
}
\details{
The super table is created by obtaining the QTL summaries for each model
in \code{\ldots} using \code{summary.wgaim()} and then row binding them
together. An extra column is created on the left hand side of the
super table for the trait names given in the \code{labels}
argument. If \code{labels = NULL} then trait names are extracted from
the left hand-side of the fixed component of the associated \code{wgaim}
model. The returned super table allows simple exporting to spreadsheet software packages
or with the R/LaTeX package \pkg{xtable}.
}
\value{A \code{data.frame} object with stacked QTL summaries
 }
 \references{
   Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
   Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
   Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{wgaim}}}
\examples{

\dontrun{

## fit wgaim models

rktgw.qtl <- wgaim(rktgw.asf, intervalObj = genoRxK, merge.by = "Genotype",
                   trace = "trace.txt", na.action = na.method(x = "include"))

rkyld.qtl <- wgaim(rkyld.asf, intervalObj = genoRxK, merge.by = "Genotype",
                   trace = "trace.txt", na.action = na.method(x = "include"))

## create super table and export

qtlt <- qtlTable(rktgw.qtl, rkyld.qtl, labels = c("Conc.", "Shoot"))
print(xtable(qtlt), file = "superQTL.tex", include.rownames = FALSE)
}
}

\keyword{regression}
