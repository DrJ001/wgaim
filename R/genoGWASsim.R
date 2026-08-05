#' Simulated raw marker data for a diversity panel
#'
#' @description
#' A simulated diversity-panel genotype dataset used by the GWAS vignettes to
#' demonstrate the complete data-preparation workflow
#' \code{\link{checkPanel}()} -> \code{\link{filterPanel}()} ->
#' \code{\link{primePanel}()}.  The panel comprises 200 diversity lines
#' genotyped at approximately 1700 SNP markers spread across 10 chromosomes at
#' approximately 0.9 cM average spacing (chromosome length ~150 cM).  All
#' markers on the panel have minor allele frequency (MAF) greater than 0.15
#' -- rare variants are excluded at the sampling stage (AlphaSimR
#' \code{addSnpChip(minSnpFreq = 0.15)}) with a belt-and-braces post-random-
#' mating filter -- so the panel is dominated by informative common variants
#' with a healthy MAF spread (median MAF ~0.32).  Genotypes are stored in
#' \strong{raw 0/1/2 allele-count coding} exactly as they would be received
#' from a genotyping facility.
#'
#' The dataset deliberately contains a range of realistic data-quality
#' issues so that the panel-preparation functions can be demonstrated
#' end-to-end.  The embedded issues are:
#'
#' \describe{
#'   \item{Missingness}{Approximately 2\% baseline missingness across all
#'     cells, plus about 150 markers with 25--40\% missingness (failed
#'     genotyping assays) and 5 lines with 22--35\% missingness (poor-quality
#'     DNA).}
#'   \item{Monomorphic markers}{Fifteen markers forced to a single allele
#'     (MAF = 0).}
#'   \item{Duplicate lines}{Four pairs of lines with identical genotype
#'     profiles (accidentally submitted twice).}
#'   \item{Duplicate markers}{Approximately 40 pairs of markers with
#'     identical genotype vectors (very-close markers in perfect LD or
#'     re-assayed loci).  The observed count in \code{\link{checkPanel}()}
#'     is typically slightly higher because of natural duplicates arising
#'     from the LD structure.}
#'   \item{Unplaceable markers}{Twenty-five markers renamed
#'     \code{SNP_UNPLACED_001}..\code{SNP_UNPLACED_025} so they no longer
#'     match any entry in \code{\link{mapGWASsim}}; these are dropped by
#'     the map-consistency step of \code{\link{filterPanel}()}.}
#'   \item{Low-MAF markers}{None below the 0.15 sampling floor.  The
#'     default \code{\link{filterPanel}()} MAF filter (0.05) is therefore
#'     a no-op on this panel, but is still worth running for real datasets.}
#' }
#'
#' The panel is drawn as a single genetically diverse population with no
#' discrete subpopulation structure.  Within each chromosome, marker
#' genotypes are generated from an AR(1) latent liability with
#' \eqn{\lambda = 3} cM, giving realistic LD decay
#' (\eqn{r(\Delta = 1\,\text{cM}) \approx 0.72},
#' \eqn{r(\Delta = 3\,\text{cM}) \approx 0.37}) but no cross-chromosome
#' correlations beyond what is captured by the whole-genome kinship
#' matrix \eqn{K}.  This is deliberate: earlier versions injected
#' independent per-chromosome subpopulation liability shifts, which
#' produced within-chromosome subpop signal stronger than \eqn{K}'s
#' cross-chromosome average and confounded the multi-environment GWAS
#' with spurious hits on whichever chromosome drew the largest subpop
#' offsets.
#'
#' Five QTL were planted to drive detectable signal in the yield phenotype
#' \code{\link{phenoGWASsim}} (see that page for details).  All five QTL
#' markers survive the default \code{\link{filterPanel}()} step so the GWAS
#' pipeline can recover them from the prepared panel.
#'
#' @format A \code{data.frame} with 200 rows (one per diversity line) and
#' approximately 1700 columns (1 id column plus one column per marker; the
#' exact marker count varies slightly across regenerations depending on how
#' many markers drift below the 0.15 MAF floor after 1 generation of random
#' mating):
#' \describe{
#'   \item{\code{id}}{Character.  Line identifier (\code{"Line001"} to
#'     \code{"Line200"}).}
#'   \item{\code{SNP_C01_0001}, \code{SNP_C01_0002}, ...}{Integer.  Raw
#'     genotype values in 0/1/2 allele-count coding (0 = AA homozygote,
#'     1 = AB heterozygote, 2 = BB homozygote, \code{NA} = missing).
#'     Marker names follow the pattern \code{SNP_C<chr>_<index>}.
#'     Twenty-five columns are renamed
#'     \code{SNP_UNPLACED_001}..\code{SNP_UNPLACED_025} to simulate markers
#'     absent from the genetic map.}
#' }
#'
#' @seealso \code{\link{mapGWASsim}}, \code{\link{phenoGWASsim}},
#'   \code{\link{checkPanel}}, \code{\link{filterPanel}},
#'   \code{\link{primePanel}}, \code{\link{gwasAim}}
#'
#' @examples
#' data(genoGWASraw,  package = "wgAim")
#' data(mapGWASsim,   package = "wgAim")
#'
#' dim(genoGWASraw)         # 200 x ~1700
#' genoGWASraw[1:5, 1:6]
#'
#' # Diagnose the embedded issues
#' \dontrun{
#' chk <- checkPanel(genoGWASraw, mapGWASsim, id = "id")
#' print(chk)
#' }
"genoGWASraw"


#' Genetic map for the simulated diversity panel
#'
#' @description
#' The genetic map accompanying \code{\link{genoGWASraw}}.  It contains
#' one row per SNP marker (approximately 1700 total) across 10 chromosomes.
#' Marker positions are in centimorgans (cM) and each chromosome spans
#' approximately 0 to 150 cM, with average marker spacing of about 0.9 cM.
#' All markers have MAF greater than 0.15 in the shipped panel.
#'
#' Note that 25 markers present in \code{\link{genoGWASraw}} have been
#' renamed \code{SNP_UNPLACED_*} and therefore do not appear in this map.
#' The map-consistency step of \code{\link{filterPanel}()} drops those
#' unplaceable markers.
#'
#' @format A \code{data.frame} with approximately 1700 rows and 3 columns:
#' \describe{
#'   \item{\code{marker}}{Character.  Marker name in the format
#'     \code{SNP_C<chr>_<index>} (e.g. \code{"SNP_C01_0001"}).}
#'   \item{\code{chr}}{Character.  Chromosome identifier
#'     (\code{"Chr1"}..\code{"Chr10"}).}
#'   \item{\code{pos}}{Numeric.  Marker position in centimorgans.}
#' }
#'
#' @seealso \code{\link{genoGWASraw}}, \code{\link{phenoGWASsim}},
#'   \code{\link{checkPanel}}, \code{\link{filterPanel}},
#'   \code{\link{primePanel}}
#'
#' @examples
#' data(mapGWASsim, package = "wgAim")
#' head(mapGWASsim)
#' table(mapGWASsim$chr)     # ~170 markers per chromosome
"mapGWASsim"


#' Simulated multi-environment trial phenotypes for the diversity panel
#'
#' @description
#' Phenotypic data for the simulated diversity panel
#' \code{\link{genoGWASraw}}: a balanced multi-environment trial (MET) with
#' 200 lines evaluated across 6 sites in a 2-replicate randomised complete
#' block design, giving 2400 rows in total.  Yield values are generated
#' from a linear mixed model with:
#'
#' \itemize{
#'   \item \strong{Fixed effects.}  Per-site intercepts drawn from
#'     \eqn{U(4.5, 4.9)} t/ha plus \code{Site:Rep} block effects.
#'   \item \strong{Polygenic background.}  TWO Kronecker-structured
#'     genotype-by-environment components matching the wgAim multivariate
#'     analysis's \code{fa(Site,k):vm(id, K) + fa(Site,k):id} decomposition:
#'     \enumerate{
#'       \item \strong{K-weighted (trace):}
#'         \eqn{\mathrm{vec}(\mathbf{U}) \sim
#'         \mathrm{MVN}(\mathbf{0}, \mathbf{G}_{a,K} \otimes \mathbf{K})},
#'         where \eqn{\mathbf{K}} is the AlphaSimR-derived kinship
#'         (standardised mean diag = 1) and \eqn{\mathbf{G}_{a,K}} is a
#'         \eqn{6 \times 6} between-site covariance from a Wishart-style
#'         construction with mean diagonal \eqn{\approx 0.005}.  Kept
#'         deliberately small so that once the planted QTL are absorbed
#'         as fixed effects the composite genome-wide LRT drops below its
#'         significance threshold.  Absorbed by \code{fa(Site,k):vm(id, K)}.
#'       \item \strong{IID line (moderate):}
#'         \eqn{\mathrm{vec}(\mathbf{L}) \sim
#'         \mathrm{MVN}(\mathbf{0}, \mathbf{G}_{a,L} \otimes \mathbf{I})},
#'         equivalently each line draws its own MVN(\eqn{\mathbf{0}},
#'         \eqn{\mathbf{G}_{a,L}}) independently, with mean diagonal
#'         \eqn{\approx 0.02}.  This is a permanent-environment / line-
#'         unique deviation that has essentially zero projection onto
#'         \eqn{\mathbf{K}}, so it does not drive the \code{vm(id, K)}
#'         LRT but keeps the \code{fa(Site,k):id} analysis term well-
#'         conditioned.  Absorbed by \code{fa(Site,k):id}.
#'     }
#'     Total polygenic variance is \eqn{\approx 0.025} (\eqn{\sim 18\%}
#'     of residual, giving broad-sense heritability \eqn{\sim 15\%}).
#'     Generating both explicitly matches the classic BLUP decomposition
#'     (relatedness contribution + line-unique deviation) and gives each
#'     analysis term a clear, disjoint job.  A pure \eqn{\mathbf{G}_a
#'     \otimes \mathbf{K}} sim -- fine for MVsim's DH-cross interval K
#'     with strong block structure -- destabilises the \code{fa(Site,k)}
#'     upgrade with a diversity-panel K (off-diagonal SD \eqn{\approx
#'     0.085}), because \code{vm(id, K)} and \code{id} become near-
#'     substitutes and their FA loadings compete for the same signal.
#'   \item \strong{Residuals.}  IID Gaussian nugget error at plot level
#'     with SD ~0.377 t/ha per site (\eqn{\sim 8\%} of the trait mean); a
#'     simple random Row / Column model in the base LMM is sufficient
#'     (no residual spatial autocorrelation).
#'   \item \strong{Planted QTL.}  Five QTL added as fixed contributions
#'     to line x site yield values so both univariate and multivariate
#'     GWAS analyses have concrete ground truth to recover.
#' }
#'
#' Yields are on the tonnes-per-hectare (t/ha) scale, matching MVsim's
#' scale of ~1.  The kg/ha parameterisation used in earlier development
#' was numerically stable in principle but drove REML gradient checks in
#' ASReml close to their default tolerances; rescaling to t/ha avoids
#' those issues without changing any of the analysis results.
#'
#' The planted QTL (given here in t/ha per allele substitution) sit at
#' realistic breeding-data magnitudes (\eqn{2{-}4\%} of the ~4.7 t/ha
#' trait mean).  Their pooled genetic variance (\eqn{\approx 0.038} per
#' site across the five loci) dominates the polygenic background
#' (\eqn{\mathrm{diag}(\mathbf{G}_{a,K} + \mathbf{G}_{a,L}) \approx
#' 0.025}), so absorbing them as fixed effects removes the bulk of the
#' K-shaped signal driving the composite genome-wide LRT.  MAIN QTL are
#' sized so that \eqn{|\beta| / \mathrm{SE}_{\mathrm{per\text{-}site}}
#' \approx 4{-}5} (SE \eqn{\approx 0.03} t/ha with 200 lines and 2 reps
#' per site), which reliably keeps a truly-constant QTL classified as
#' MAIN by the internal Wald pruning step.  QTL 3 is deliberately left
#' small to demonstrate a borderline-detectable effect.
#'
#' \tabular{lllll}{
#'   \strong{QTL} \tab \strong{Chr} \tab \strong{cM} \tab \strong{Type} \tab \strong{Effect pattern} \cr
#'   1 \tab Chr2 \tab  45 \tab MAIN         \tab \eqn{+0.145}  t/ha across all 6 sites \cr
#'   2 \tab Chr5 \tab  90 \tab MAIN         \tab \eqn{-0.145}  t/ha across all 6 sites \cr
#'   3 \tab Chr7 \tab  30 \tab MAIN (weak)  \tab \eqn{+0.125}  t/ha across all 6 sites \cr
#'   4 \tab Chr3 \tab 110 \tab G x E fade   \tab \eqn{+0.230} in Env1 fading linearly to 0 in Env6 \cr
#'   5 \tab Chr9 \tab  60 \tab G x E cross  \tab \eqn{-0.180} in Env1--Env3, \eqn{+0.180} in Env4--Env6 \cr
#' }
#'
#' The corresponding QTL marker names (as they appear in
#' \code{\link{genoGWASraw}} and \code{\link{mapGWASsim}}) for the current
#' shipped simulation are \code{SNP_C02_0044}, \code{SNP_C05_0094},
#' \code{SNP_C07_0042}, \code{SNP_C03_0128}, and \code{SNP_C09_0068}
#' respectively.  These marker indices are not stable across re-runs of
#' \code{data-raw/simulate_GWASsim.R} because the underlying AlphaSimR
#' coalescent samples segregating sites at random -- the planted cM
#' positions (Chr2 @ 45, Chr5 @ 90, Chr7 @ 30, Chr3 @ 110, Chr9 @ 60) are
#' stable but the closest well-behaved marker chosen by the sim script can
#' shift by a few indices.
#'
#' The single-environment (univariate) GWAS vignette uses
#' \code{phenoGWASsim[phenoGWASsim$Site == "Env1", ]}.  The multi-environment
#' (multivariate) vignette uses all 6 sites with \code{Trait = "Site"}.
#'
#' @format A \code{data.frame} with 2400 rows and 6 columns:
#' \describe{
#'   \item{\code{Site}}{Factor with 6 levels
#'     (\code{"Env1"}..\code{"Env6"}).  Environment identifier.}
#'   \item{\code{id}}{Factor with 200 levels
#'     (\code{"Line001"}..\code{"Line200"}).  Line identifier.  This
#'     column name matches the \code{id} column of
#'     \code{\link{genoGWASraw}} and is the natural \code{merge.by}
#'     argument for \code{\link{gwasAim}()}.}
#'   \item{\code{Rep}}{Factor with levels \code{"1"} and \code{"2"}.
#'     Replicate within site.}
#'   \item{\code{Row}, \code{Column}}{Integer.  Plot coordinates within
#'     each site trial layout.  Can be modelled as random main effects
#'     to absorb row / column trend.}
#'   \item{\code{yield}}{Numeric.  Simulated grain yield in t/ha.
#'     Site means and per-site residual SDs vary within the ranges
#'     typical of a MET trial.}
#' }
#'
#' @seealso \code{\link{genoGWASraw}}, \code{\link{mapGWASsim}},
#'   \code{\link{gwasAim}}, \code{\link{primePanel}}
#'
#' @examples
#' data(phenoGWASsim, package = "wgAim")
#' head(phenoGWASsim)
#' table(phenoGWASsim$Site, phenoGWASsim$Rep)
#' tapply(phenoGWASsim$yield, phenoGWASsim$Site, mean)
"phenoGWASsim"
