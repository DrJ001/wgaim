# =============================================================================
# data-raw/simulate_GWASsim.R
#
# Generates the three package datasets used by the GWAS vignettes:
#   genoGWASraw    -- raw 0/1/2 marker data.frame WITH DATA-QUALITY ISSUES
#                     (missingness, monomorphic markers, duplicate lines and
#                     markers, unplaceable markers) so that the full pipeline
#                     checkPanel() -> filterPanel() -> primePanel() can be
#                     demonstrated end-to-end.
#   mapGWASsim     -- genetic map (marker, chr, pos in cM).  25 markers in
#                     genoGWASraw are renamed to SNP_UNPLACED_* so they no
#                     longer match the map (forces map-consistency filtering).
#   phenoGWASsim   -- balanced 6-site MET phenotypes with 5 planted QTL
#                     added as fixed contributions on top of an IID line-
#                     level polygenic term.
#
# Population generation
# ---------------------
# The panel is simulated with AlphaSimR's coalescent-based runMacs2() and a
# single generation of random mating.  This produces a diversity panel whose
# kinship matrix K has REALISTIC full-rank structure with a gradually
# decaying eigenvalue spectrum -- unlike a hand-rolled independent-draw
# construction (K nearly identity, vm(id, K) redundant with id) or a
# rank-few PC hack (K rank-deficient, all markers project onto the same
# 3-5 latent directions and forward selection fails to discriminate real
# QTL from noise).
#
# Concretely:
#   * 50 founder haplotypes via runMacs2 (Ne = 100, genLen = 1.5 Morgan ->
#     150 cM per chromosome), which produces realistic LD blocks and
#     segregating-site frequency spectra from a neutral coalescent.
#   * 1 generation of random mating (200 crosses among the 50 founders)
#     produces 200 diverse F1 lines with half-sib / full-sib relationships
#     -- the K matrix then has real off-diagonal structure that vm(id, K)
#     can exploit.
#   * ~150 SNPs per chromosome sampled from segregating sites with founder-
#     pop MAF > 0.15 (via minSnpFreq in AlphaSimR's addSnpChip); actively
#     drawing from the well-behaved end of the MAF spectrum stops the panel
#     being polluted by very-rare variants that produce spurious associations
#     in GWAS.  ~1500 total panel markers, further filtered by a post-random-
#     mating F1 MAF > 0.15 check (typically drops a small handful due to drift).
#
# Yield scale: TONNES per hectare (t/ha), not kg/ha.  Matching MVsim's
# scale (~1) avoids numerical drift in REML gradients that showed up with
# the previous 4700-kg/ha parameterisation: even though the analysis is
# invariant to global rescaling, ASReml's convergence heuristics and
# per-parameter tolerances behave more reliably when variance components
# sit near unit scale.
#
# Planted QTL (5 total; effects in t/ha, site mean ~ 4.7 t/ha):
#   QTL 1  Chr2   45 cM   MAIN         +0.145 t/ha across all 6 sites
#   QTL 2  Chr5   90 cM   MAIN         -0.145 t/ha across all 6 sites
#   QTL 3  Chr7   30 cM   MAIN weak    +0.125 t/ha across all 6 sites
#   QTL 4  Chr3  110 cM   GxE fade     +0.230 (Env1) fading linearly to 0 (Env6)
#   QTL 5  Chr9   60 cM   GxE cross    -0.180 (Env1-3), +0.180 (Env4-6)
#
# MAIN QTL are sized so that |beta| / per-site-SE ~ 4-5 (with n=200 lines and
# 2 reps per site the per-site SE for a single-marker fixed effect is ~0.03
# t/ha).  This is the minimum ratio that reliably keeps a truly-constant
# QTL classified as MAIN by .waldTest() -- with weaker effects, random per-
# site sampling variation across 6 sites can look like real G x E and the
# QTL gets promoted to INT.  The G x E QTL are sized so their peak per-site
# |beta| gives per-site LOD ~5-7 (comparable to the MAIN QTL) so the G x E
# pattern is unambiguous, not marginal.  QTL 3 (MAIN weak) is deliberately
# left at |beta|=0.100 so it demonstrates a borderline-detectable effect.
#
# Response model (independent of AlphaSimR's built-in trait mechanisms):
#   y_ijr = mu_j + rep_{jr} + sum_k beta_jk X_ik + U_ij + L_ij + eps_ijr
# where
#   * beta_jk X_ik are the planted-QTL fixed contributions
#   * vec(U) ~ MVN(0, Ga_K (x) K)  -- K-WEIGHTED polygenic (relatedness-
#     driven).  Ga_K is a 6 x 6 between-site covariance from a Wishart
#     construction; K is the AlphaSimR-derived kinship (mean diag = 1).
#     Absorbed by fa(Site,k):vm(id, K) in the analysis.
#   * vec(L) ~ MVN(0, Ga_L (x) I)  -- IID-LINE polygenic (line-unique
#     deviation from the K-predicted value).  Ga_L is a 6 x 6 covariance
#     also generated Wishart-style.  Equivalently: each line draws its
#     own MVN(0, Ga_L) independently.  Absorbed by fa(Site,k):id in the
#     analysis.
#   * eps_ijr ~ N(0, 0.142)  ->  residual SD ~0.377 t/ha per site (~8%
#     of trait mean).
#
# Why BOTH U and L?  With a diversity-panel K (off-diagonal SD ~0.085),
# fa(Site,k):vm(id, K) and fa(Site,k):id parameterise near-substitute
# 6x6 site covariances on the same line factor.  If the sim generates
# ONLY U (pure Ga (x) K) the analysis id term has no legitimate variance
# to fit and competes with vm for the same K-shaped signal, which
# destabilises the fa upgrade.  Generating BOTH components matches the
# classic BLUP decomposition (relatedness contribution + line-unique
# deviation) and gives each analysis term a clear disjoint job.  MVsim
# survives with only U because its DH-cross interval K has strong block
# structure that separates cleanly from I; our diversity panel K does
# not.
#
# Polygenic magnitudes (deliberately SMALL so the composite genome-wide
# LRT stops firing once all true QTL have been absorbed as fixed effects):
#
#   mean diag(Ga_K) ~ 0.005  -- trace K-weighted variance; provides just
#                               enough K-shape signal on top of the QTL
#                               marker columns to fire the initial LRT,
#                               but sub-threshold once the QTL are removed.
#   mean diag(Ga_L) ~ 0.020  -- moderate IID-line variance; keeps
#                               fa(Site,k):id from bottoming out at zero
#                               (which would destabilise the FA upgrade)
#                               without swamping the vm signal.
#
# Total polygenic variance is ~0.025 (vs residual 0.142 -> broad-sense h^2
# ~15%).  QTL variance sums to ~0.014 across the 5 planted loci, so QTL
# heritability is ~9% and TOTAL genetic heritability is ~22%.  These
# calibrations are empirical and may need tuning if the LRT still fires
# spuriously after all true QTL are absorbed -- shrink Ga_K further if so.
#
# Data-quality issues embedded in genoGWASraw (all raw, pre-filter):
#   ~150 markers with 25-40% missingness
#   ~5 lines   with 22-35% missingness
#   ~2% baseline missing rate across all cells
#   ~15 monomorphic markers (all AA)
#   ~4 duplicate line pairs
#   ~40 duplicate marker pairs
#   ~25 markers in geno but not in map (renamed to SNP_UNPLACED_*)
#
# Run interactively (needs AlphaSimR and wgAim available):
#   source("data-raw/simulate_GWASsim.R")
# =============================================================================

suppressPackageStartupMessages({
    library(AlphaSimR)
})
devtools::load_all()   # loads wgAim

set.seed(20260805)

# ---- Panel & genome parameters --------------------------------------------
n.lines    <- 200L
n.chr      <- 10L
n.mar.pc   <- 175L                  # markers per chromosome requested from
                                    # AlphaSimR (final count drops slightly
                                    # after the F1 MAF > 0.15 filter below)
maf.floor  <- 0.15                  # MAF floor applied both at panel-
                                    # selection time (via minSnpFreq in
                                    # addSnpChip) and post-hoc on the F1
                                    # panel to catch any drift from random
                                    # mating
n.founder  <- 50L
n.gen.mate <- 1L                    # random-mating generations after founders
chr.len.cM <- 150                   # cM per chromosome (genLen 1.5 Morgan)
n.sites    <- 6L
n.reps     <- 2L

# ---- Line IDs / env names -------------------------------------------------
w    <- floor(log10(n.lines)) + 1L
ids  <- sprintf(paste0("Line%0", w, "d"), seq_len(n.lines))

# biomAid-compatible site names padded to the minimum width needed
site.width <- nchar(as.character(n.sites))
env.names  <- sprintf(paste0("Env%0", site.width, "d"), seq_len(n.sites))

# ---- Simulate founder haplotypes via coalescent (runMacs2) ----------------
# runMacs2 gives full control over Ne, mutation rate, recombination rate.
# Ne = 100 (realistic effective population size for a diverse panel).
# genLen = 1.5 Morgan matches our chromosome length of 150 cM.
message(sprintf("Simulating %d founder haplotypes via runMacs2()...",
                n.founder))
founderPop <- runMacs2(
    nInd     = n.founder,
    nChr     = n.chr,
    segSites = 2000,                 # plenty for a 150-SNP chip per chr
    Ne       = 100,
    genLen   = chr.len.cM / 100      # Morgans (1.5 = 150 cM)
)

# ---- SimParam: register a n.mar.pc SNPs-per-chromosome panel --------------
# minSnpFreq = maf.floor tells AlphaSimR to ACTIVELY sample the SNP chip
# only from segregating sites whose founder-pop MAF exceeds the floor,
# rather than picking any segregating site at random (the default).  This
# stops the panel from being dominated by very-low-MAF markers that (a)
# have high per-marker score variance, (b) don't align well with K
# structure, and (c) tend to produce spurious associations in GWAS.
SP <- SimParam$new(founderPop)
SP$nThreads <- 1L                    # deterministic across runs
SP$addSnpChip(nSnpPerChr = n.mar.pc, minSnpFreq = maf.floor)

# ---- Expand to n.lines lines via random mating ----------------------------
# 300 random crosses among 50 founders -> 300 F1 lines, each with two
# founders as parents.  Half-sib / full-sib relationships give K real
# off-diagonal structure (as opposed to K ~= I for independent draws).
message(sprintf("Random-mating (%d generation%s) to expand to %d lines...",
                n.gen.mate, ifelse(n.gen.mate == 1L, "", "s"), n.lines))
pop <- newPop(founderPop, simParam = SP)
for (g in seq_len(n.gen.mate)) {
    pop <- randCross(pop, nCrosses = n.lines, simParam = SP)
}
stopifnot(nInd(pop) == n.lines)

# ---- Extract 0/1/2 genotypes and map --------------------------------------
geno.mat <- pullSnpGeno(pop, simParam = SP)   # n.lines x (n.chr * n.mar.pc)
storage.mode(geno.mat) <- "integer"
rownames(geno.mat) <- ids

# getSnpMap returns data.frame with columns id, chr (character), site (integer),
# pos (numeric, in Morgans on this chromosome).  Sort by chromosome then
# position and re-order geno.mat columns to match.
snp.map.raw         <- getSnpMap(simParam = SP)
snp.map.raw$chr.int <- as.integer(snp.map.raw$chr)
snp.map.raw         <- snp.map.raw[order(snp.map.raw$chr.int,
                                         snp.map.raw$pos), ]
geno.mat            <- geno.mat[, snp.map.raw$id, drop = FALSE]
n.pre               <- ncol(geno.mat)
message(sprintf("Raw panel from AlphaSimR: %d markers across %d chromosomes",
                n.pre, n.chr))

# ---- F1 MAF filter --------------------------------------------------------
# addSnpChip(minSnpFreq = maf.floor) selects the chip from founder-pop MAF >
# maf.floor, but 1 generation of random mating can shift some markers below
# that threshold by drift.  Recompute MAFs on the F1 panel and drop any
# below the floor to guarantee every shipped marker has F1 MAF > maf.floor.
f1_af  <- colMeans(geno.mat, na.rm = TRUE) / 2
f1_maf <- pmin(f1_af, 1 - f1_af)
keep   <- f1_maf > maf.floor
if (any(!keep)) {
    geno.mat    <- geno.mat[, keep, drop = FALSE]
    snp.map.raw <- snp.map.raw[keep, , drop = FALSE]
    f1_maf      <- f1_maf[keep]
    message(sprintf(
        "F1 MAF filter: dropped %d markers with F1 MAF <= %.2f",
        n.pre - ncol(geno.mat), maf.floor))
}
message(sprintf(
    "F1 MAF distribution:  min=%.3f  Q1=%.3f  median=%.3f  Q3=%.3f  max=%.3f",
    min(f1_maf), quantile(f1_maf, 0.25),
    median(f1_maf), quantile(f1_maf, 0.75), max(f1_maf)))
message(sprintf("Markers per chromosome after filter:"))
print(table(snp.map.raw$chr.int))

# ---- Re-index within chromosomes and build the shipped map ---------------
# After filtering some markers, re-number within each chromosome so marker
# names remain contiguous (SNP_CXX_0001, 0002, 0003, ... with no gaps).
within.chr.idx <- unlist(lapply(split(seq_len(nrow(snp.map.raw)),
                                       snp.map.raw$chr.int), seq_along),
                          use.names = FALSE)
mapGWASsim <- data.frame(
    marker = sprintf("SNP_C%02d_%04d", snp.map.raw$chr.int, within.chr.idx),
    chr    = paste0("Chr", snp.map.raw$chr.int),
    pos    = snp.map.raw$pos * 100,   # Morgans -> cM
    stringsAsFactors = FALSE
)
colnames(geno.mat) <- mapGWASsim$marker
message(sprintf("Shipped panel: %d lines x %d markers across %d chromosomes",
                nrow(geno.mat), ncol(geno.mat), n.chr))

# ---- Compute kinship K from the clean panel --------------------------------
# Standard VanRaden-style K: centre each marker column, then K = M M' / p,
# rescaled so mean(diag(K)) = 1.  Used both for the diagnostic block below
# AND for the polygenic-term Kronecker draw further down.
{
    M.raw <- geno.mat - 1L
    storage.mode(M.raw) <- "double"
    M.c   <- sweep(M.raw, 2L, colMeans(M.raw), "-")
    K     <- tcrossprod(M.c) / ncol(M.c)
    K     <- K / mean(diag(K))
    rm(M.raw, M.c)
}

# ---- Diagnostic: inspect K off-diagonal structure --------------------------
{
    off.tri <- K[lower.tri(K)]
    ev      <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    ev      <- sort(ev, decreasing = TRUE)
    message(sprintf(
        "K diag range [%.3f, %.3f], off-diag mean = %+.3f, SD = %.3f, range [%.3f, %.3f]",
        min(diag(K)), max(diag(K)),
        mean(off.tri), sd(off.tri), min(off.tri), max(off.tri)))
    message(sprintf(
        "K spectrum: top 5 = [%s], min = %.4f, top-5/total = %.1f%%",
        paste(round(head(ev, 5), 2), collapse = ", "),
        min(ev), 100 * sum(head(ev, 5)) / sum(ev)))
    rm(off.tri, ev)
}

# ---- Select QTL marker columns from the CLEAN matrix ----------------------
# Pick a well-behaved marker (MAF > 0.25, no missing yet) as close as
# possible to each planted QTL position.  A floor of 0.25 is used (rather
# than the more permissive 0.15) to keep the per-QTL variance contribution
# (2 * p * (1 - p) * beta^2) stable across sim regenerations -- with a
# 0.15 floor the coalescent occasionally lands the "closest" marker on
# MAF ~0.15 (per-site variance ~1.5x smaller than at MAF ~0.30), which
# can drop a small-effect QTL below the LRT detection threshold and
# force yet another round of effect-size retuning.  Trading a few cM of
# marker-to-target-position drift for MAF stability is worthwhile here.
qtl.plan <- list(
    list(chr = "Chr2", pos =  45,  type = "MAIN",     eff = rep( 0.145,      n.sites)),
    list(chr = "Chr5", pos =  90,  type = "MAIN",     eff = rep(-0.145,      n.sites)),
    list(chr = "Chr7", pos =  30,  type = "MAIN-wk",  eff = rep( 0.125,      n.sites)),
    list(chr = "Chr3", pos = 110,  type = "GxE-fade",
         eff = seq(0.230, 0, length.out = n.sites)),
    list(chr = "Chr9", pos =  60,  type = "GxE-cross",
         eff = c(rep(-0.180, n.sites %/% 2L), rep( 0.180, n.sites - n.sites %/% 2L)))
)
n.qtl <- length(qtl.plan)

find_qtl_marker <- function(chr, pos, geno.mat, map, maf.floor = 0.25) {
    mk.idx <- which(map$chr == chr)
    dist   <- abs(map$pos[mk.idx] - pos)
    ordr   <- order(dist)
    for (i in ordr) {
        j <- mk.idx[i]
        g <- geno.mat[, j]
        p <- mean(g, na.rm = TRUE) / 2
        if (min(p, 1 - p) > maf.floor) return(j)
    }
    # No marker on this chromosome hit maf.floor: pick the highest-MAF marker
    # available (rare fallback -- with 150 SNPs/chr and Ne=100 coalescent this
    # essentially never triggers at maf.floor = 0.25).
    mafs <- pmin(colMeans(geno.mat[, mk.idx, drop = FALSE], na.rm = TRUE) / 2,
                 1 - colMeans(geno.mat[, mk.idx, drop = FALSE], na.rm = TRUE) / 2)
    mk.idx[which.max(mafs)]
}
qtl.mk.idx <- vapply(qtl.plan,
                      function(q) find_qtl_marker(q$chr, q$pos, geno.mat, mapGWASsim),
                      integer(1L))
qtl.mk.names <- mapGWASsim$marker[qtl.mk.idx]

# QTL score matrix in +/-1 coding for effect application (heterozygotes = 0)
qtl.X <- geno.mat[, qtl.mk.idx, drop = FALSE] - 1L
storage.mode(qtl.X) <- "double"
qtl.eff <- do.call("rbind", lapply(qtl.plan, function(q) q$eff))
rownames(qtl.eff) <- qtl.mk.names
colnames(qtl.eff) <- env.names
qtl.contrib <- qtl.X %*% qtl.eff       # n.lines x n.sites (t/ha)

message("\nPlanted QTL summary:")
for (i in seq_len(n.qtl)) {
    q <- qtl.plan[[i]]; mk <- qtl.mk.names[i]
    p <- mean(geno.mat[, qtl.mk.idx[i]], na.rm = TRUE) / 2
    cat(sprintf("  [%s] %s @ %.1f cM  marker=%-16s  MAF=%.3f  eff=[%s]\n",
                q$type, q$chr, q$pos, mk, min(p, 1 - p),
                paste(sprintf("%+.3f", q$eff), collapse = ",")))
}

# ---- Between-site genetic covariance Ga (Wishart-style) --------------------
# Same construction pattern as MVsim: draw a small MVN sample from a
# compound-symmetric seed matrix S, take crossprod / df to get positive
# definite Wishart-distributed matrices, then rescale to a target
# mean-diagonal.  We build TWO 6 x 6 covariances so that the response
# has two disjoint polygenic components matching the wgAim MV analysis:
#
#   Ga_K  --  between-site covariance for the K-weighted term U
#             (analysis: fa(Site,k):vm(id, K))
#   Ga_L  --  between-site covariance for the IID-line term L
#             (analysis: fa(Site,k):id)
#
# Total polygenic variance ~ 0.025 (~18% of residual).  Deliberately
# SMALL so the composite genome-wide LRT stops firing once all true QTL
# have been absorbed as fixed effects (see header comment for the
# calibration argument).  Both use their own Wishart seed (independent
# set.seed) so their cross-site correlation structures differ, but both
# are moderate positive corgh.
target.GaK.diag <- 0.005                     # K-weighted per-site variance (trace)
target.GaL.diag <- 0.020                     # IID-line per-site variance (moderate)

# ---- Ga_K (K-weighted between-site covariance) ----------------------------
set.seed(42)
S.scale.K       <- matrix(0.55, n.sites, n.sites)
diag(S.scale.K) <- 1.00
Ga_K.raw <- MASS::mvrnorm(n.sites + 2L, mu = rep(0, n.sites), Sigma = S.scale.K)
Ga_K     <- crossprod(Ga_K.raw) / (n.sites + 2L)
Ga_K     <- Ga_K * (target.GaK.diag / mean(diag(Ga_K)))
dimnames(Ga_K) <- list(env.names, env.names)

message("\nTrue Ga_K diagonal (K-weighted per-site variance):")
print(setNames(round(diag(Ga_K), 5), env.names))
message("True Ga_K correlations:")
Dk.inv <- diag(1 / sqrt(diag(Ga_K)))
print(round(Dk.inv %*% Ga_K %*% Dk.inv, 3))

# ---- Ga_L (IID-line between-site covariance) ------------------------------
# Fresh seed and slightly different compound-symmetric prior so Ga_L has
# a distinct correlation pattern from Ga_K (both moderate positive).
# Seed changed from 1729 -> 1730: the 1729 Wishart draw put Env3's diagonal
# ~2.6x above the target (0.052 vs 0.020), which persistently added extra
# noise to Env3 observations and caused MAIN QTL on Chr5 to be misclassified
# as G x E because the Env3 effect estimate wandered while the other sites
# stayed near the true constant.
set.seed(1730)
S.scale.L       <- matrix(0.45, n.sites, n.sites)
diag(S.scale.L) <- 1.00
Ga_L.raw <- MASS::mvrnorm(n.sites + 2L, mu = rep(0, n.sites), Sigma = S.scale.L)
Ga_L     <- crossprod(Ga_L.raw) / (n.sites + 2L)
Ga_L     <- Ga_L * (target.GaL.diag / mean(diag(Ga_L)))
dimnames(Ga_L) <- list(env.names, env.names)

message("\nTrue Ga_L diagonal (IID-line per-site variance):")
print(setNames(round(diag(Ga_L), 5), env.names))
message("True Ga_L correlations:")
Dl.inv <- diag(1 / sqrt(diag(Ga_L)))
print(round(Dl.inv %*% Ga_L %*% Dl.inv, 3))

# ---- K-weighted polygenic  U ~ MVN(0, Ga_K (x) K) -------------------------
# Single draw of length n.lines * n.sites; reshape to n.lines x n.sites
# so U[i, j] is line i's K-driven genetic value at site j.  This puts
# real K-shape into the response so fa(Site,k):vm(id, K) has genome-wide
# signal to detect at iteration 1.
set.seed(101)
message(sprintf("\nDrawing U ~ MVN(0, Ga_K (x) K) of length %d ...",
                n.lines * n.sites))
Sigma_U <- kronecker(Ga_K, K)                # (n.lines*n.sites)^2  -- 1800^2
u.vec   <- MASS::mvrnorm(1L, mu = rep(0, n.lines * n.sites), Sigma = Sigma_U)
U       <- matrix(u.vec, nrow = n.lines, ncol = n.sites,
                  dimnames = list(ids, env.names))
rm(Sigma_U, u.vec)
message(sprintf("Empirical per-site variance of U: %s",
                paste(sprintf("%.5f", apply(U, 2L, var)), collapse = ", ")))
message(sprintf("Empirical mean U cross-site correlation: %.3f",
                mean(cor(U)[upper.tri(cor(U))])))

# ---- IID-line polygenic  L ~ MVN(0, Ga_L (x) I) ---------------------------
# kronecker(Ga_L, I_nlines) is block-diagonal with Ga_L in every block,
# so this is equivalent to drawing each line's per-site vector
# independently from MVN(0, Ga_L).  L is absorbed by fa(Site,k):id in the
# analysis (line-unique deviation, no K weighting).
set.seed(202)
message(sprintf("\nDrawing L ~ MVN(0, Ga_L (x) I) of length %d ...",
                n.lines * n.sites))
L <- MASS::mvrnorm(n.lines, mu = rep(0, n.sites), Sigma = Ga_L)
dimnames(L) <- list(ids, env.names)
message(sprintf("Empirical per-site variance of L: %s",
                paste(sprintf("%.5f", apply(L, 2L, var)), collapse = ", ")))
message(sprintf("Empirical mean L cross-site correlation: %.3f",
                mean(cor(L)[upper.tri(cor(L))])))

# ---- Site means and residual variances ------------------------------------
set.seed(4321)
mu.site  <- runif(n.sites, 4.5, 4.9)
names(mu.site) <- env.names
sigma2.e <- rep(0.142, n.sites)      # residual SD ~0.377 t/ha (~8% of mean)
message(sprintf("\nSite means: %s", paste(round(mu.site, 3), collapse = ", ")))

# ---- Assemble phenoGWASsim ------------------------------------------------
# Balanced 2-rep layout per site: 10 rows x 20 cols per rep, two reps
# stacked so Row 1..10 = Rep 1 and Row 11..20 = Rep 2 (20 unique Row
# levels, 20 unique Column levels per site, 400 plots per site).
n.rows.rep <- 10L
n.cols     <- n.lines %/% n.rows.rep
stopifnot(n.rows.rep * n.cols == n.lines)

set.seed(9182)
phenoGWASsim <- do.call("rbind", lapply(seq_len(n.sites), function(j) {
    q.j <- qtl.contrib[, j]   # named by id  (planted QTL contribution)
    U.j <- U[, j]             # named by id  (K-weighted polygenic, Ga_K x K)
    L.j <- L[, j]             # named by id  (IID-line polygenic, Ga_L x I)
    do.call("rbind", lapply(seq_len(n.reps), function(r) {
        perm     <- sample.int(n.lines)  # randomise plot -> line
        rows.vec <- (r - 1L) * n.rows.rep +
                    rep(seq_len(n.rows.rep), each = n.cols)
        cols.vec <- rep(seq_len(n.cols), n.rows.rep)
        line.ids <- ids[perm]
        data.frame(
            Site   = factor(env.names[j], levels = env.names),
            id     = factor(line.ids, levels = ids),
            Rep    = factor(as.character(r), levels = as.character(seq_len(n.reps))),
            Row    = as.integer(rows.vec),
            Column = as.integer(cols.vec),
            yield  = mu.site[j] + q.j[line.ids] + U.j[line.ids] + L.j[line.ids] +
                     rnorm(n.lines, 0, sqrt(sigma2.e[j])),
            stringsAsFactors = FALSE
        )
    }))
}))
phenoGWASsim <- phenoGWASsim[
    order(phenoGWASsim$Site, phenoGWASsim$Rep,
          phenoGWASsim$Row, phenoGWASsim$Column), ]
row.names(phenoGWASsim) <- NULL

message(sprintf("\nPheno assembled: %d rows  (%d sites x %d lines x %d reps)",
                nrow(phenoGWASsim), n.sites, n.lines, n.reps))
message(sprintf("Yield range: [%.3f, %.3f]  overall mean = %.3f",
                min(phenoGWASsim$yield), max(phenoGWASsim$yield),
                mean(phenoGWASsim$yield)))

# ---- Introduce data-quality issues ----------------------------------------
# QTL columns are protected from all subsequent issue injection so ground
# truth stays intact.
protected.cols <- qtl.mk.idx

# 1. Baseline ~2% missingness across the whole non-QTL matrix
non.qtl.cols <- setdiff(seq_len(ncol(geno.mat)), protected.cols)
n.cells <- n.lines * length(non.qtl.cols)
n.base  <- round(0.02 * n.cells)
flat.idx <- sample.int(n.cells, n.base)
row.idx  <- ((flat.idx - 1L) %% n.lines) + 1L
col.idx  <- non.qtl.cols[((flat.idx - 1L) %/% n.lines) + 1L]
geno.mat[cbind(row.idx, col.idx)] <- NA_integer_
message(sprintf("\n  Baseline missing:  %d cells (%.2f%%)",
                n.base, 100 * n.base / (n.lines * ncol(geno.mat))))

# 2. High-missingness markers (~150 markers, 25-40% miss each)
hi.miss.mk <- sample(non.qtl.cols, 150L)
for (j in hi.miss.mk) {
    mr    <- runif(1, 0.25, 0.40)
    nmiss <- rbinom(1L, n.lines, mr)
    geno.mat[sample.int(n.lines, nmiss), j] <- NA_integer_
}
message(sprintf("  High-miss markers: %d  (target 25-40%% each)", length(hi.miss.mk)))

# 3. High-missingness lines (5 lines, 22-35% miss each)
hi.miss.ln <- sample.int(n.lines, 5L)
for (i in hi.miss.ln) {
    mr    <- runif(1, 0.22, 0.35)
    ncols <- length(non.qtl.cols)
    nmiss <- rbinom(1L, ncols, mr)
    geno.mat[i, sample(non.qtl.cols, nmiss)] <- NA_integer_
}
message(sprintf("  High-miss lines:   %d  (target 22-35%% each)", length(hi.miss.ln)))

# 4. Monomorphic markers (~15 markers, forced all-AA)
mono.mk <- sample(setdiff(non.qtl.cols, hi.miss.mk), 15L)
geno.mat[, mono.mk] <- 0L
message(sprintf("  Monomorphic:       %d markers (forced AA)", length(mono.mk)))

# 5. Duplicate lines (4 pairs -- second of each pair copied from the first)
dup.line.donor <- sample.int(n.lines, 4L)
dup.line.copy  <- sample(setdiff(seq_len(n.lines), dup.line.donor), 4L)
for (k in seq_along(dup.line.donor))
    geno.mat[dup.line.copy[k], ] <- geno.mat[dup.line.donor[k], ]
message(sprintf("  Duplicate lines:   %d pairs (donor -> copy)", length(dup.line.donor)))

# 6. Duplicate markers (40 pairs -- avoid QTL, high-miss, mono; copy AFTER
#    line duplication so both copies remain identical)
avail <- setdiff(non.qtl.cols, c(hi.miss.mk, mono.mk))
dup.mk.donor <- sample(avail, 40L)
dup.mk.copy  <- sample(setdiff(avail, dup.mk.donor), 40L)
for (k in seq_along(dup.mk.donor))
    geno.mat[, dup.mk.copy[k]] <- geno.mat[, dup.mk.donor[k]]
message(sprintf("  Duplicate markers: %d pairs (donor -> copy)", length(dup.mk.donor)))

# 7. Unplaceable markers (~25 markers renamed so they no longer match the map)
unplaced.idx <- sample(setdiff(non.qtl.cols,
                                c(hi.miss.mk, mono.mk, dup.mk.donor, dup.mk.copy)),
                        25L)
new.names <- sprintf("SNP_UNPLACED_%03d", seq_along(unplaced.idx))
colnames(geno.mat)[unplaced.idx] <- new.names
message(sprintf("  Unplaceable:       %d markers (renamed SNP_UNPLACED_*)",
                length(unplaced.idx)))

# ---- Assemble raw genoGWASraw data.frame -----------------------------------
genoGWASraw <- data.frame(id = ids,
                           as.data.frame(geno.mat, check.names = FALSE),
                           check.names      = FALSE,
                           stringsAsFactors = FALSE)

# ---- Save -----------------------------------------------------------------
usethis::use_data(genoGWASraw, mapGWASsim, phenoGWASsim, overwrite = TRUE)

# ---- Verification summary --------------------------------------------------
cat("\n============================================================\n")
cat("  Simulation summary\n")
cat("============================================================\n")
cat(sprintf("  Lines            : %d\n", nrow(genoGWASraw)))
cat(sprintf("  Marker cols (raw): %d\n", ncol(genoGWASraw) - 1L))
cat(sprintf("  Markers in map   : %d\n", nrow(mapGWASsim)))
cat(sprintf("  NA in geno       : %d cells (%.2f%%)\n",
            sum(is.na(genoGWASraw[, -1L])),
            100 * mean(is.na(genoGWASraw[, -1L]))))
cat(sprintf("  Pheno rows       : %d  (%d sites x %d lines x %d reps)\n",
            nrow(phenoGWASsim), n.sites, n.lines, n.reps))
cat(sprintf("  Yield range      : [%.3f, %.3f]  mean=%.3f\n",
            min(phenoGWASsim$yield), max(phenoGWASsim$yield),
            mean(phenoGWASsim$yield)))
cat(sprintf("  Ga_K diag        : [%.5f, %.5f]  mean = %.5f\n",
            min(diag(Ga_K)), max(diag(Ga_K)), mean(diag(Ga_K))))
cat(sprintf("  Ga_K mean r      : %.3f\n",
            mean((diag(1 / sqrt(diag(Ga_K))) %*% Ga_K %*%
                  diag(1 / sqrt(diag(Ga_K))))[upper.tri(Ga_K)])))
cat(sprintf("  Ga_L diag        : [%.5f, %.5f]  mean = %.5f\n",
            min(diag(Ga_L)), max(diag(Ga_L)), mean(diag(Ga_L))))
cat(sprintf("  Ga_L mean r      : %.3f\n",
            mean((diag(1 / sqrt(diag(Ga_L))) %*% Ga_L %*%
                  diag(1 / sqrt(diag(Ga_L))))[upper.tri(Ga_L)])))
cat(sprintf("  Polygenic form   : U ~ MVN(0, Ga_K (x) K) + L ~ MVN(0, Ga_L (x) I)\n"))
cat(sprintf("  Residual SD      : %.3f t/ha per site (variance %.4f)\n",
            sqrt(sigma2.e[1L]), sigma2.e[1L]))

cat("\n  Planted QTL (marker actually chosen):\n")
for (i in seq_len(n.qtl)) {
    q <- qtl.plan[[i]]
    mk <- qtl.mk.names[i]
    p <- mean(geno.mat[, qtl.mk.idx[i]], na.rm = TRUE) / 2
    cat(sprintf("    %-12s %-4s @ %5.1f cM  marker=%-16s  MAF=%.3f\n",
                paste0("[", q$type, "]"), q$chr, q$pos, mk, min(p, 1 - p)))
}

cat("\n  Files written:\n")
cat("    data/genoGWASraw.rda\n")
cat("    data/mapGWASsim.rda\n")
cat("    data/phenoGWASsim.rda\n")
cat("============================================================\n")
