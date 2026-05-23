# =============================================================================
# data-raw/simulate_MVsim.R
#
# Generates genoMVsim and phenoMVsim — package datasets for the multivariate
# QTL mapping vignette.
#
# Population: 120 DH lines, 6 chromosomes x 100 markers @ 1 cM spacing.
#   594 intervals > 120 lines => vm() engine path throughout.
#
# Design: 8 field trials, 2 replicates per trial.
#   FA(2) for 8 trials: 8*2 + 8 = 24 parameters < unstructured 36.
#
# Planted QTL (4 total):
#   QTL 1  Chr1, interval 20  — main effect (consistent across all trials)
#   QTL 2  Chr3, interval 40  — main effect (consistent, opposite sign)
#   QTL 3  Chr2, interval 55  — G×E interaction (fades from T1 to T8)
#   QTL 4  Chr5, interval 70  — G×E interaction (crossover: negative T1-T4,
#                                positive T5-T8)
#
# Effect sizes ~3-4% of the trait mean — realistic and detectable but not
# dominant relative to the polygenic background.
#
# Run interactively (with devtools::load_all() or wgAim installed):
#   source("data-raw/simulate_MVsim.R")
# =============================================================================

library(qtl)
library(MASS)
devtools::load_all()   # loads wgAim for primeCross()

set.seed(2025)

# ---- Population parameters --------------------------------------------------
n.lines  <- 120L
n.chr    <- 6L
n.mar    <- 100L       # 99 intervals/chr x 6 = 594 > 120 lines → vm path
n.trials <- 8L         # FA(2): 24 params < unstructured 36
n.reps   <- 2L
ids      <- paste0("DH", sprintf("%03d", seq_len(n.lines)))

# ---- Simulate DH cross ------------------------------------------------------
map.sim <- sim.map(len       = rep((n.mar - 1L), n.chr),
                   n.mar     = rep(n.mar, n.chr),
                   eq.spac   = TRUE,
                   include.x = FALSE,
                   sex.sp    = FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))

genoMVsim <- sim.cross(map.sim, n.ind = n.lines, type = "bc",
                        error.prob = 0.005)
genoMVsim$pheno$id <- ids

# ---- Interval data for background and QTL scores ----------------------------
genObj <- primeCross(genoMVsim, type = "interval", id = "id")

M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids

K <- tcrossprod(M) / ncol(M)
K <- K / mean(diag(K))      # standardise so mean diagonal = 1

# ---- Planted QTL positions and genotype scores ------------------------------
qtl.info <- list(
    list(chr = "Chr1", int = 20L, type = "MAIN"),
    list(chr = "Chr3", int = 40L, type = "MAIN"),
    list(chr = "Chr2", int = 55L, type = "INT"),
    list(chr = "Chr5", int = 70L, type = "INT")
)
n.qtl <- length(qtl.info)

qtl.X <- matrix(0, nrow = n.lines, ncol = n.qtl,
                dimnames = list(ids,
                    sapply(qtl.info, function(q) paste(q$chr, q$int, sep = "."))))
for (i in seq_len(n.qtl))
    qtl.X[, i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]

# ---- QTL effect matrix (n.qtl x n.trials) -----------------------------------
# QTL1 & QTL2: constant across all 8 trials (main effects).
# QTL3: fades from strong (T1) to near zero (T8) — partial G×E.
# QTL4: symmetric crossover — negative T1-T4, positive T5-T8.
qtl.eff <- matrix(
    c( 0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  0.35,  # QTL1: main
      -0.30, -0.30, -0.30, -0.30, -0.30, -0.30, -0.30, -0.30,  # QTL2: main
       0.30,  0.26,  0.20,  0.14,  0.08,  0.05,  0.03,  0.02,  # QTL3: G×E fade
      -0.25, -0.18, -0.10, -0.03,  0.03,  0.10,  0.18,  0.25), # QTL4: crossover
    nrow    = n.qtl,
    ncol    = n.trials,
    byrow   = TRUE,
    dimnames = list(colnames(qtl.X), paste0("Trial", seq_len(n.trials)))
)

cat("Planted QTL effects per trial:\n")
print(round(qtl.eff, 3))
cat("\n")

# ---- Genetic (co)variance matrix Ga (Wishart-based PD) ----------------------
set.seed(42)
S.scale       <- matrix(0.035, n.trials, n.trials)
diag(S.scale) <- 0.055
Ga.raw <- MASS::mvrnorm(n.trials + 2L, mu = rep(0, n.trials), Sigma = S.scale)
Ga     <- crossprod(Ga.raw) / (n.trials + 2L)
Ga     <- Ga * (0.025 / mean(diag(Ga)))

cat("True Ga diagonal (genetic variances per trial):\n")
print(round(diag(Ga), 4))
cat("True genetic correlations:\n")
D.inv <- diag(1 / sqrt(diag(Ga)))
print(round(D.inv %*% Ga %*% D.inv, 3))
cat("\n")

# ---- Polygenic background from MVN(0, Ga ⊗ K) ------------------------------
set.seed(101)
u.vec <- MASS::mvrnorm(1, mu = rep(0, n.lines * n.trials),
                       Sigma = kronecker(Ga, K))
U <- matrix(u.vec, nrow = n.lines, ncol = n.trials,
            dimnames = list(ids, paste0("Trial", seq_len(n.trials))))

# ---- Site means and residual variances --------------------------------------
set.seed(42)
sigma2.e <- rep(0.80, n.trials)
mu.site  <- runif(n.trials, 9.5, 12.5)
cat("Site means:", round(mu.site, 2), "\n\n")

# ---- Assemble phenoMVsim ----------------------------------------------------
phenoMVsim <- do.call(rbind, lapply(seq_len(n.trials), function(j) {
    qtl.contrib <- rowSums(sweep(qtl.X, 2, qtl.eff[, j], "*"))
    do.call(rbind, lapply(seq_len(n.reps), function(r)
        data.frame(
            id    = factor(ids),
            Trial = factor(paste0("Trial", j),
                           levels = paste0("Trial", seq_len(n.trials))),
            Rep   = factor(paste0("R", r)),
            yield = mu.site[j] + qtl.contrib + U[, j] +
                    rnorm(n.lines, 0, sqrt(sigma2.e[j])),
            stringsAsFactors = FALSE
        )
    ))
}))
phenoMVsim <- phenoMVsim[order(phenoMVsim$Trial, phenoMVsim$Rep, phenoMVsim$id), ]
rownames(phenoMVsim) <- NULL

cat("pheno dimensions:", nrow(phenoMVsim), "rows x", ncol(phenoMVsim), "cols\n")
cat("yield range:", round(range(phenoMVsim$yield), 2), "\n")

# ---- Save -------------------------------------------------------------------
usethis::use_data(genoMVsim, phenoMVsim, overwrite = TRUE)
message("Saved data/genoMVsim.rda and data/phenoMVsim.rda")
