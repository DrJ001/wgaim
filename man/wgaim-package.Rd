\name{wgaim-package}
\alias{wgaim-package}
\docType{package}
\title{Whole Genome Average Interval Mapping (wgaim) for QTL detection
  and estimation
}
\description{This package provides an efficient computational implementation of the QTL
whole genome QTL analysis algorithm (wgaim) discussed in Verbyla et. al
(2007, 2012) using extensions of the functionality provided in the
linear mixed modelling R package \pkg{ASReml-R} V4.
}
\details{
\tabular{ll}{
Package: \tab wgaim\cr
Type: \tab Package\cr
Version: \tab 2.0-0\cr
Date: \tab 2019-08-12\cr
License: \tab GPL 2\cr
}
Welcome to version 2.0 of \pkg{wgaim}! The documentation given in this help file is
only brief and users should consult the vignette available with the
package by typing \code{vignette("wgaim")} at the prompt. Alternatively,
users can also consult the individual help files for the main functions of the
package.

The package provides a user friendly function \code{\link{cross2int}}
for the conversion of \code{"cross"} objects created using
\code{read.cross} in Bromans \pkg{qtl} package into an \code{"interval"}
object ready for use in \pkg{wgaim}. Specifically, \code{cross2int()} performs additional
calculations to impute missing marker values on each of the
chromosomes across the full linkage map and also provides users with
genetic distances and recombination fractions for the intervals. The
returned object retains the class structure
of an object created with \code{read.cross} and therefore allows further
use with the \pkg{qtl} package if desired.

The package also provides a function for the graphical display of the
chromosomes of a \code{"cross"} object. The method function
\code{\link{linkMap}} displays the full or subsetted linkage map
according to chromosome or distance as well as displays non-overlapping
marker names on the right hand side.

QTL analysis is conducted using the function \code{\link{wgaim}} which,
as its first argument, requires an \code{asreml} base model. High
dimensional genetic components are allowed (See \code{wgaim.asreml} for
more details). For convenience the default tracing of results from the
asreml models can be outputted to a file for further inspection. Outlier
statistics and marker effects from each iteration can be viewed using
\code{\link{outStat}}. Diagnostics of the likelihood ratio test
performed for each forward step can be displayed using
\code{\link{tr.wgaim}}. The function also displays an incremental
probability value matrix of the QTL ascertained at each forward step of the algorithm.

Summary and print methods are available for the returned \code{"wgaim"}
object and provide users with a detailed report on the QTL, their
size, their flanking markers, significance (including LOD score if
desired) and approximate contribution to the genetic variance. The
returned \code{"wgaim"} object may also be plotted using the method
function \code{\link{linkMap}}. This function plots the full linkage map
subsetted for chromosome and distance as well as provides shaded
QTL regions and highlighted flanking markers. Plotting of QTL for
multiple traits is also possible (see \code{\link{linkMap.default}}).

}
\author{Julian Taylor and Ari Verbyla
  Maintainer: Julian Taylor <julian.taylor@adelaide.edu.au>
}
\references{
 Verbyla, A. P & Taylor, J. D, Verbyla, K. L (2012). RWGAIM: An efficient high
  dimensional random whole genome average (QTL) interval mapping approach.
  \emph{Genetics Research}. \bold{94}, 291-306.

 Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
 Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
 Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.

 Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
 by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 95-111.
}
\keyword{ package }
\seealso{
  \code{\link[qtl]{qtl-package}}
}

