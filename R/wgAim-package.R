#' wgAim: Whole Genome Analyses via Integrated Modelling
#'
#' @description
#' The \pkg{wgAim} package provides a unified framework for three whole-genome
#' analysis methods in plant breeding populations, all implemented through the
#' ASReml-R linear mixed modelling engine:
#'
#' \describe{
#'   \item{\code{\link{QTLAim}}}{Forward-selection QTL detection in biparental
#'     mapping populations using whole-genome average interval mapping.}
#'   \item{\code{\link{GWASAim}}}{Forward-selection genome-wide association
#'     analysis in diversity panels using the same composite genome-wide model.}
#'   \item{\code{\link{GPAim}}}{Genomic best linear unbiased prediction
#'     (G-BLUP) delivering genomic estimated breeding values (GEBVs) via either
#'     a genomic relationship matrix (vm path) or direct marker effect
#'     estimation (mbf path).}
#' }
#'
#' @details
#' All three analyses share a common internal engine built around ASReml-R.
#' The unifying statistical idea is that the full complement of markers or
#' intervals is included simultaneously in the model as a composite
#' genome-wide random effect, represented either as a genomic relationship
#' matrix (\code{vm} path, when markers exceed lines) or as a marker-by-file
#' random effect (\code{mbf} path, when lines equal or exceed markers).
#'
#' \strong{Data preparation:}
#' \itemize{
#'   \item Biparental mapping populations: use \code{\link[wgaim]{cross2int}}
#'     from the legacy \code{wgaim} functions to prepare an \code{"interval"}
#'     object for \code{QTLAim} or \code{GPAim}.
#'   \item Diversity panels: use \code{\link{makePanel}} to prepare a
#'     \code{"panel"} object for \code{GWASAim} or \code{GPAim}.
#' }
#'
#' \strong{Legacy functions:}
#' The original \code{\link{wgaim}} QTL analysis function and its associated
#' utilities (\code{cross2int}, \code{fixMap}, \code{fineMap}, \code{linkMap},
#' \code{outStat}, etc.) remain fully available for backward compatibility.
#'
#' @references
#' Verbyla, A.P., Taylor, J.D. and Verbyla, K.L. (2012). RWGAIM: An efficient
#' high-dimensional random whole genome average (QTL) interval mapping approach.
#' \emph{Genetics Research}, \bold{94}, 291--306.
#'
#' Taylor, J. and Verbyla, A.P. (2011). R Package wgaim: QTL Analysis in
#' Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
#' Statistical Software}, \bold{40}(7), 1--18.
#'
#' @name wgAim-package
#' @aliases wgAim
#'
#' @import stats utils graphics grDevices qtl ggplot2
#'
#' @rawNamespace
#'   export(wgaim)
#'   export(wgaim.asreml)
#'   export(wgaim.default)
#'   export(cross2int)
#'   export(fineMap)
#'   export(fixMap)
#'   export(getQTL)
#'   export(linkMap)
#'   export(linkMap.cross)
#'   export(linkMap.default)
#'   export(linkMap.wgaim)
#'   export(outStat)
#'   export(print.wgaim)
#'   export(tr.wgaim)
#'   export(qtlSelect)
#'   export(qtlTable)
#'   export(summary.wgaim)
#'   S3method(summary, wgaim)
#'   S3method(print, wgaim)
#'   S3method(linkMap, wgaim)
#'   S3method(tr, wgaim)
#'   S3method(linkMap, cross)
#'   S3method(linkMap, default)
#'   S3method(wgaim, asreml)
#'   S3method(wgaim, default)
NULL
