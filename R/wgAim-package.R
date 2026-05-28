#' wgAim: Whole Genome Analyses via Integrated Modelling
#'
#' @description
#' The \pkg{wgAim} package provides a unified framework for whole-genome
#' analysis in plant breeding populations, built entirely around the
#' ASReml-R linear mixed modelling engine. It extends and supersedes the
#' legacy \pkg{wgaim} package, broadening the scope from QTL detection alone
#' into three integrated analytical pillars:
#'
#' \describe{
#'   \item{\code{\link{qtlAim}}}{Forward-selection QTL detection in biparental
#'     mapping populations using whole-genome average interval mapping. Supports
#'     fixed and random QTL effects, interval or marker selection, and
#'     multivariate models via the \code{Trait} argument (trials, environments,
#'     sites, treatments, or measured traits).}
#'   \item{\code{\link{gwasAim}}}{Forward-selection genome-wide association
#'     analysis in diversity panels. Uses the same composite genome-wide model
#'     as \code{qtlAim} with marker effects always treated as fixed. Supports
#'     multivariate models via the \code{Trait} grouping argument.}
#'   \item{\code{\link{gpAim}}}{Genomic best linear unbiased prediction
#'     (G-BLUP) delivering genomic estimated breeding values (GEBVs) and
#'     BLUP accuracies via either a genomic relationship matrix (\code{vm}
#'     path) or direct marker-by-file random effects (\code{mbf} path).
#'     Supports univariate and multivariate models via \code{Trait}.}
#' }
#'
#' The genomic prediction pillar is completed by \code{\link{selIndex}}, which
#' combines per-environment GEBVs from a multivariate \code{gpAim} fit into a
#' single selection index (weighted, Smith-Hazel, or Pesek-Baker
#' desired-gains), computes expected genetic gain, and provides diagnostic
#' plots for adaptation screening.
#'
#' @details
#' \strong{Shared engine:}
#' All three analyses share a common internal engine. The unifying statistical
#' idea is that the full complement of markers or intervals is included
#' simultaneously in the mixed model as a composite genome-wide random effect,
#' represented either as a genomic relationship matrix
#' (\code{vm} path, when markers exceed lines) or as a marker-by-file
#' random effect (\code{mbf} path, when lines equal or exceed markers).
#'
#' \strong{Multivariate extensions (\code{Trait} argument):}
#' All three analysis functions accept an optional \code{Trait} argument
#' naming a factor column in the phenotypic data that defines the grouping
#' structure for multivariate analysis — trials, environments, sites,
#' treatments, or measured traits. When supplied, a multivariate model is
#' fitted using mixture chi-squared LRT for significance testing and
#' Wald test pruning to classify QTL as main effects or group-by-genotype
#' interactions.
#' The variance structure of the additive genomic term is controlled by the
#' \code{str} argument (\code{"corh"}, \code{"fa1"}, etc.).
#'
#' \strong{Data preparation:}
#' \itemize{
#'   \item Biparental mapping populations: use \code{\link{primeCross}} to
#'     prepare a \code{"wgCross"} object from an R/qtl cross object.
#'     \code{\link{fixMap}} and \code{imputeGen()} are available as
#'     helpers for map and genotype quality control.
#'   \item Diversity panels: use \code{\link{primePanel}} to prepare a
#'     \code{"wgPanel"} object from a raw marker matrix.
#' }
#'
#' \strong{Supporting functions:}
#' \itemize{
#'   \item \code{\link{fineMap}} -- high-resolution scan around a detected QTL
#'     or association.
#'   \item \code{\link{aimTable}} -- stack summary tables from multiple
#'     \code{qtlAim} or \code{gwasAim} objects (multi-trait comparison).
#'   \item \code{\link{aimTrace}} -- trace the forward-selection algorithm:
#'     LRT statistics, QTL effect stability plots.
#'   \item \code{\link{linkMap}} -- genetic linkage map plots annotated with
#'     detected QTL or GWAS markers.
#'   \item \code{\link{getQTL}} -- extract QTL position and flanking marker
#'     information from a fitted object.
#' }
#'
#' \strong{Bundled datasets:}
#' Three biparental wheat populations with paired genotype/phenotype data are
#' included: RAC875 x Kukri (\code{genoRxK} / \code{phenoRxK}), Cascades x
#' RAC875-2 (\code{genoCxR} / \code{phenoCxR}), and Sunco x Tasman
#' (\code{genoSxT} / \code{phenoSxT}).
#'
#' @references
#' Verbyla, A.P., Cullis, B.R. and Thompson, R. (2007). The analysis of QTL
#' by simultaneous use of the full linkage map. \emph{Theoretical and Applied
#' Genetics}, \bold{116}, 95--111.
#'
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
#' @importFrom Rcpp sourceCpp
#' @useDynLib wgAim, .registration = TRUE
NULL
