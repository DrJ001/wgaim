#' Simulated DH genotype data for multivariate QTL mapping
#'
#' @description
#' A simulated Doubled Haploid (DH) mapping population for use in the
#' multivariate QTL mapping vignette. The population consists of 120 DH lines
#' genotyped at 100 markers on each of 6 chromosomes (600 markers total) at 1 cM
#' spacing, giving 594 marker intervals.
#'
#' The object is of class \code{"bc"} as produced by \code{qtl::sim.cross()} and
#' is intended to be passed directly to \code{\link{primeCross}()}.
#'
#' Four QTL were planted at the following positions:
#' \describe{
#'   \item{Chr1, interval 20}{Main effect — consistent positive effect across
#'     all eight trials.}
#'   \item{Chr3, interval 40}{Main effect — consistent negative effect across
#'     all eight trials.}
#'   \item{Chr2, interval 55}{G×E interaction — strong positive effect in
#'     early trials, fading to near zero by trial 8.}
#'   \item{Chr5, interval 70}{G×E interaction — symmetric crossover: negative
#'     in trials 1–4, positive in trials 5–8.}
#' }
#'
#' @format An object of class \code{c("bc", "cross")} with components:
#' \describe{
#'   \item{\code{$geno}}{Named list of 6 chromosome objects, each containing
#'     \code{$data} (marker genotype matrix) and \code{$map} (cM positions).}
#'   \item{\code{$pheno}}{Data frame with columns \code{id} (line identifier)
#'     and \code{phenotype} (placeholder; phenotypic data are in
#'     \code{\link{phenoMVsim}}).}
#' }
#'
#' @seealso \code{\link{phenoMVsim}}, \code{\link{primeCross}},
#'   \code{\link{qtlAim}}
#'
#' @examples
#' data(genoMVsim, package = "wgAim")
#' length(genoMVsim$geno)       # 6 chromosomes
#' nmar(genoMVsim)              # 100 markers per chromosome
"genoMVsim"

#' Simulated phenotypic data for multivariate QTL mapping
#'
#' @description
#' Phenotypic data for the simulated DH mapping population
#' \code{\link{genoMVsim}}, comprising 120 lines evaluated across eight field
#' trials in a two-replicate design (1920 observations total).
#'
#' Yield was simulated from a multivariate linear mixed model with:
#' \itemize{
#'   \item A polygenic background drawn from
#'     \eqn{\mathcal{N}(\mathbf{0},\, \mathbf{G}_a \otimes \mathbf{K})},
#'     where \eqn{\mathbf{G}_a} is an 8 \eqn{\times} 8 genetic (co)variance
#'     matrix and \eqn{\mathbf{K}} is the genomic relationship matrix.
#'   \item Four planted QTL (two main-effect, two G×E interaction) with effect
#'     sizes of approximately 3–4\% of the trait mean. See
#'     \code{\link{genoMVsim}} for QTL positions and types.
#'   \item Trial-specific residual variances of 0.80 and site means ranging
#'     from approximately 10 to 12.5 yield units.
#' }
#'
#' @format A data frame with 1920 rows and 4 columns:
#' \describe{
#'   \item{\code{id}}{Factor. Line identifier (\code{"DH001"} to
#'     \code{"DH120"}).}
#'   \item{\code{Trial}}{Factor with levels \code{Trial1} to \code{Trial8}.}
#'   \item{\code{Rep}}{Factor with levels \code{R1} and \code{R2}.}
#'   \item{\code{yield}}{Numeric. Simulated grain yield (arbitrary units).}
#' }
#'
#' @seealso \code{\link{genoMVsim}}, \code{\link{primeCross}},
#'   \code{\link{qtlAim}}
#'
#' @examples
#' data(phenoMVsim, package = "wgAim")
#' head(phenoMVsim)
#' tapply(phenoMVsim$yield, phenoMVsim$Trial, mean)
"phenoMVsim"
