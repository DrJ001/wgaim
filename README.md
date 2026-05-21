# wgAim — Whole Genome Analyses via Integrated Modelling

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![Package version](https://img.shields.io/github/r-package/v/DrJ001/wgAim?color=orange&label=Package%20version&style=flat-square)](https://github.com/DrJ001/wgAim/commits/dev)
[![codecov](https://codecov.io/gh/DrJ001/wgAim/branch/dev/graph/badge.svg)](https://app.codecov.io/gh/DrJ001/wgAim)
<!-- badges: end -->

**Authors**: Julian Taylor & Ari Verbyla

> **This package is under active development.** The API may change
> between commits ahead of the first stable release. Three vignettes —
> one each for QTL mapping, GWAS, and genomic prediction & selection —
> are in preparation and will provide full worked examples of each
> pipeline.

---

## Overview

**wgAim** is a unified whole-genome analysis toolkit for plant breeding
populations, built entirely around the
[ASReml-R](https://vsni.co.uk/software/asreml-r/) linear mixed modelling
engine. It extends and supersedes the
[wgaim](https://github.com/DrJ001/wgaim-legacy) package, broadening
the scope from QTL detection alone into a three-pillar framework for
whole-genome inference:

| Pillar | Functions | Description |
|---|---|---|
| QTL mapping | `qtlAim()` | Forward-selection whole-genome marker or interval mapping in biparental populations |
| GWAS | `gwasAim()` | Forward-selection genome-wide association analysis in diversity panels |
| Genomic prediction & selection | `gpAim()` · `selIndex()` | G-BLUP delivering GEBVs and accuracies; `selIndex()` combines per-environment GEBVs into a weighted, Smith-Hazel, or desired-gains selection index |

All three analyses share a common internal engine and the same
unifying statistical idea: the full complement of markers or intervals
enters the mixed model simultaneously as a composite genome-wide random
effect, represented either as a genomic relationship matrix (`vm` path)
or as a direct marker-by-file random effect (`mbf` path). All methods
support multivariate (multi-environment) analysis through a `Trait`
argument that extends the model to capture genotype-by-environment
interaction.

Each pipeline follows the same three-stage workflow: **pre-analysis**
data preparation, **analysis**, and **post-analysis** interpretation.

---

## Part 1 — Pre-analysis: Data Preparation

Before any analysis can be run, raw genotype data must be prepared into
the structured objects that the analysis functions expect. **wgAim**
provides two preparation functions, one for each population type.

### Biparental mapping populations — `primeCross()`

`primeCross()` takes a `cross` object produced by the
[qtl](https://cran.r-project.org/package=qtl) package and returns a
`wgCross` object ready for `qtlAim()` or `gpAim()`. It handles
marker imputation and constructs the interval genotype scores used in
interval mapping.

```r
wgCross <- primeCross(cross, type = "interval", id = "Genotype")
```

Key arguments:

| Argument | Description |
|---|---|
| `type` | `"interval"` (default) for interval mapping using Haldane-weighted midpoint scores; `"marker"` for direct marker analysis |
| `id` | Name of the line identifier column in `cross$pheno` |
| `impute` | Missing genotype imputation method: `"Martinez"` (default) or `"Martini"` |

Two supporting functions assist with data quality:

- **`fixMap()`** — corrects marker ordering on a linkage map when
  markers have been placed out of sequence.
- **`imputeGen()`** — performs standalone genotype imputation on a
  `cross` object prior to calling `primeCross()`, useful when finer
  control over imputation is needed.

### Diversity panels — `primePanel()`

`primePanel()` prepares genotype and map data from a diversity panel
into a `wgPanel` object for `gwasAim()` or `gpAim()`.

```r
wgPanel <- primePanel(geno, map, id = "Genotype",
                      map.id = "marker", map.chr = "chr", map.pos = "pos",
                      encoding = "012", maf = 0.05)
```

Key arguments:

| Argument | Description |
|---|---|
| `encoding` | Input coding of the genotype matrix: `"012"` (default, AA/AB/BB counts) or `"-101"` (already centred) |
| `maf` | Minor allele frequency threshold; markers below this are removed |
| `impute` | Missing genotype imputation: `"none"` (default), `"mean"`, or `"mode"` |

Both `primeCross()` and `primePanel()` recode genotypes internally to
the ±1 scale used by the analysis engine and store the prepared
matrices in `$imputed.data` (and `$interval.data` for interval type).

---

## Part 2 — Analysis

Each analysis function takes a fitted **base ASReml-R model** (capturing
the experimental design) and a prepared genotype object, then
progressively builds a whole-genome model using forward selection.

### QTL mapping — `qtlAim()`

Designed for biparental mapping populations (`wgCross`). Performs
forward selection of QTL as fixed effects against a composite
genome-wide random background, iterating until no further significant
QTL are detected.

```r
qtl.fit <- qtlAim(baseModel, genObj = wgCross, merge.by = "Genotype",
                  TypeI = 0.05, selection = "interval")
```

Key arguments:

| Argument | Description |
|---|---|
| `selection` | `"interval"` (default) or `"marker"` or `"chromosome"` |
| `method` | `"fixed"` (default) or `"random"` — whether QTL enter as fixed or random effects |
| `gen.type` | `"interval"` or `"marker"` — genotype data source from `wgCross`; inferred from `primeCross()` type if `NULL` |
| `exclusion.window` | cM exclusion zone around detected QTL during search (default 20) |
| `TypeI` | Significance threshold for the likelihood ratio test (default 0.05) |
| `Trait` | Column name of a trial/environment factor for multivariate G×E analysis |
| `str` | Variance structure on the genomic term for multivariate models (e.g. `"corh"`, `"fa1"`) |

### GWAS — `gwasAim()`

Designed for diversity panels (`wgPanel`). Uses the same forward-
selection engine as `qtlAim()` but scans actual panel markers and
applies a panel-size-adjusted significance threshold.

```r
gwas.fit <- gwasAim(baseModel, genObj = wgPanel, merge.by = "Genotype",
                    TypeI = 0.05)
```

The interface is deliberately close to `qtlAim()`. The `method` and
`selection` arguments are fixed internally (`"fixed"` and `"interval"`
respectively); all other arguments are shared.

### Genomic prediction — `gpAim()`

Fits a whole-genome G-BLUP model and extracts genomic estimated
breeding values (GEBVs), prediction accuracies, and generalised
heritabilities for all genotyped lines. Works with both `wgCross`
and `wgPanel` objects.

```r
gp.fit <- gpAim(baseModel, genObj = wgCross, merge.by = "Genotype",
                gen.type = "marker")
```

Key outputs stored in `$GP`:

| Slot | Content |
|---|---|
| `$gebv` | Data frame: line, GEBV, SE, accuracy (Mrode formula), generalised H² (Cullis) |
| `$var.genetic` | Estimated additive genetic variance |
| `$heritability` | Narrow-sense h² |
| `$rel.matrix` | Genomic relationship matrix (vm path) |

### Selection index — `selIndex()`

Combines per-environment GEBVs from a multivariate `gpAim()` result
into a single selection index. Three index types are available:

```r
si <- selIndex(gp.fit, weights = c(E1 = 1, E2 = 2),
               type = "smith-hazel", prop.select = 0.10)
```

| `type` | Method |
|---|---|
| `"weighted"` | Weighted sum of GEBVs: **I** = **b′g** |
| `"smith-hazel"` | Classical Smith-Hazel index: **b** = **P**⁻¹**Ga w** |
| `"desired-gains"` | Pesek-Baker desired-gains index: **b** ∝ **Ga**⁻¹**d** |

The `$gain` slot reports the expected genetic gain per environment for
the selected proportion of lines.

---

## Part 3 — Post-analysis: Interpretation and Display

Once an analysis is complete, a consistent set of display and
interrogation functions is available for all three methods.

### Print and summary

All result objects have `print()` and `summary()` methods that provide
progressively more detail:

- **`print()`** — compact iteration-by-iteration digest: detected
  QTL/markers, effect sizes, and percentage variance explained.
- **`summary()`** — full tabular summary of the final model: chromosome
  positions, allele effects ± SE, percentage variance explained, and
  LOD scores. For multivariate models, results are reported per trial
  with a column indicating whether each QTL is a main effect or a
  G×E interaction.

### Plot methods

All three analysis objects support a `plot()` method with a `type`
argument selecting the display:

| `type` | Available for | Description |
|---|---|---|
| `"outlier"` | `qtlAim` | Genome-wide outlier statistics faceted by iteration |
| `"manhattan"` | `gwasAim` | Manhattan plot with chromosome background shading |
| `"blups"` | `qtlAim`, `gwasAim` | Scaled whole-genome BLUPs across the genome |
| `"effects"` | `qtlAim`, `gwasAim` | Lollipop plot: allele effect ± SE and % variance explained per QTL |
| `"contrast"` | `qtlAim`, `gwasAim` | Total genetic value scatter by allele class (requires phenotypic data) |
| `"heatmap"` | `qtlAim`, `gwasAim` | Genome × iteration heatmap of outlier statistics |
| `"chr"` | `qtlAim` | Chromosome-level bar chart (chromosome selection only) |
| `"caterpillar"` | `gpAim` | GEBV caterpillar plot with accuracy-based error bars |
| `"accuracy"` | `gpAim` | Distribution of prediction accuracies |
| `"index"` | `selIndex` | Index value scatter with selection threshold |
| `"heatmap"` | `selIndex` | GEBV × environment heatmap for selected lines |
| `"weights"` | `selIndex` | Index weight bar chart |

For multivariate models, `"effects"` and `"contrast"` produce grouped
displays with one panel or line per trial level.

### Iteration diagnostics — `aimTrace()`

`aimTrace()` traces the forward-selection iterations of a `qtlAim` or
`gwasAim` fit, showing how QTL were progressively added and how their
effects evolved.

```r
aimTrace(qtl.fit, plot = "lrt")       # LRT statistic at each iteration
aimTrace(qtl.fit, plot = "stability") # QTL effect ± 1 SE across iterations
aimTrace(qtl.fit, plot = "both")      # returns both plots as a list
```

For multivariate models, the stability plot shows one trajectory per
trial level.

### Fine mapping — `fineMap()`

`fineMap()` performs high-resolution single-QTL scanning in a window
around each detected QTL or GWAS marker, with all other detected
effects held in the model.

```r
fm <- fineMap(qtl.fit, genObj = wgCross, window = 50, step = 2)
print(fm)
plot(fm)
```

For `qtlAim` models a dense inferred grid is built across the window
using Haldane interpolation weights. For `gwasAim` models the actual
panel markers within the window are scanned. Both univariate and
multivariate models are supported: main-effect QTL use the standard
z-ratio test while interaction QTL use a joint Wald test across all
trial-level coefficients, with LOD scores on a common scale.

### Linkage map display — `linkMap()`

`linkMap()` produces publication-quality chromosome diagrams and has
methods for raw cross objects, fitted `qtlAim` results, and fitted
`gwasAim` results:

```r
linkMap(wgCross)                         # raw map
linkMap(qtl.fit, genObj = wgCross)       # map with QTL positions
linkMap(gwas.fit, genObj = wgPanel)      # map with GWAS marker positions
```

### Extracting and tabulating results — `getQTL()` and `aimTable()`

- **`getQTL()`** — extracts a tidy data frame of detected QTL/markers
  with positions, effect sizes, and percentage variance explained.
- **`aimTable()`** — formats results into a publication-ready table
  showing flanking markers, inferred positions, effect sizes, and
  LOD scores. Prints cleanly to the console and is straightforward
  to export.

---

## Vignettes

Three vignettes are in preparation, providing full worked examples from
data preparation through to interpretation:

| Vignette | Pipeline |
|---|---|
| *QTL mapping with wgAim* | `primeCross()` → `qtlAim()` → display and fine mapping |
| *GWAS with wgAim* | `primePanel()` → `gwasAim()` → display and fine mapping |
| *Genomic prediction and selection with wgAim* | `primeCross()` / `primePanel()` → `gpAim()` → `selIndex()` |

---

## Requirements

All three analysis functions require a valid **ASReml-R V4** licence,
available from [vsni.co.uk](https://vsni.co.uk/software/asreml-r/).

---

## Installation

Once the package reaches a stable release it will be submitted to CRAN.
In the meantime, the development version can be installed directly from
GitHub:

```r
# install.packages("devtools")
devtools::install_github("DrJ001/wgAim", ref = "dev")
```

---

## Legacy package

**wgAim** is the successor to
[wgaim](https://github.com/DrJ001/wgaim-legacy). The original `wgaim`
package — including its CRAN release and full QTL analysis
functionality — is preserved at the link above and remains fully
functional. Users relying on `wgaim` in existing workflows are
encouraged to stay on that version until `wgAim` reaches a stable
release.

---

## References

Verbyla, A.P., Cullis, B.R. & Thompson, R. (2007) The analysis of QTL
by simultaneous use of the full linkage map. *Theoretical and Applied
Genetics*, **116**, 95–111.

Verbyla, A.P., Taylor, J.D. & Verbyla, K.L. (2012) RWGAIM: An efficient
high-dimensional random whole genome average (QTL) interval mapping
approach. *Genetics Research*, **94**, 291–306.

Taylor, J. & Verbyla, A.P. (2011) R Package wgaim: QTL Analysis in
Bi-Parental Populations Using Linear Mixed Models. *Journal of
Statistical Software*, **40**(7), 1–18.
