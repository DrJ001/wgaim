\name{outStat}
\alias{outStat}
\title{A faceted ggplot of the chromosome outlier statistics or the
  interval blups/outlier statistics obtained from specified iteratons of wgaim.}
\description{
  A faceted \code{ggplot()} of the chromosome outlier statistics or the
  interval blups/outlier statistics from specified iteratons of
  \code{wgaim}. The interval blups/outlier statistics appear as
a trace across the genome separated by chromosomes and appropriately
spaced by their cM distances.
}
\usage{
outStat(object, intervalObj, iter = NULL, chr = NULL, statistic =
        "outlier", plot.chr = FALSE, chr.lines = FALSE)
}
\arguments{
  \item{object}{object of class \code{"wgaim"}.
  }
  \item{intervalObj}{object of class \code{"interval"}.
  }
  \item{iter}{range of integers determining which iterations will be
    plotted.
  }
  \item{chr}{character vector naming the subset of chromosomes to
    plot.
  }
  \item{statistic}{character string naming the type of diagnostic statistic to be
    plotted. Default is \code{"outlier"} (outlier statistics). Other option is
    \code{"blups"} for the scaled empirical blups calculated during each
    iteration.
  }
  \item{plot.chr}{logical value, if \code{TRUE} then plot chromosome outlier
    statistics. If \code{FALSE} then plot interval outlier
    statistics (see Details). Defaults to \code{FALSE}.
  }
  \item{chr.lines}{logical value, if \code{TRUE} then plot vertical
     lines to show separation of linkage groups. This is only useful if
     \code{plot.chr = FALSE}. Defaults to \code{FALSE}.
   }
}
\details{
If \code{plot.chr = TRUE} then outlier statistics for each chromosome
are plotted in separate faceted panels for specified values of \code{chr} and
\code{iter}. This option requies \code{selection="chromosome"} to be set
in the \code{wgaim.asreml()} call. If \code{plot.chr = FALSE} then interval blups or
outlier statistics are plotted in separate faceted panels for specified
values of \code{chr} and \code{iter}.

Additionally, the set of significant QTL (chromosome and interval position) are
extracted from the model \code{object} and annotated on the plot in
their appropriate positions in each facet panel. Graphical aesthetics,
such as themes, text, font etc. can be further manipulated through
the inclusion of additional overlays to the returned \code{ggplot()} object.
}
\value{The blups or outlier statistics are plotted in a faceted \code{ggplot()}
  with information of significant QTL overlayed.
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
\seealso{\code{\link{tr.wgaim}}, \code{\link{wgaim}}}
\examples{

\dontrun{
# fit wgaim model

rkyld.qtl <- wgaim(rkyld.asf, intervalObj = genoRxK, merge.by = "Genotype",
                  trace = "trace.txt", na.action = na.method(x = "include"))

# plot QTL interval outlier statistics

outStat(rkyld.qtl, genoRxK, iter = 1:5)

}
}
\keyword{hplot}
