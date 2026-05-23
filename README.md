# wgAim — Whole Genome Analyses via Integrated Modelling

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![Package version](https://img.shields.io/github/r-package/v/DrJ001/wgAim?color=orange&label=Package%20version&style=flat-square)](https://github.com/DrJ001/wgAim/commits/dev)
[![codecov](https://codecov.io/gh/DrJ001/wgAim/branch/dev/graph/badge.svg)](https://app.codecov.io/gh/DrJ001/wgAim)
<!-- badges: end -->

**Authors**: Julian Taylor & Ari Verbyla

> **This package is under active development.** The API may change
> between commits ahead of the first stable release. A vignette for the
> QTL mapping pipeline is now available; vignettes for GWAS and genomic
> prediction & selection are in preparation.

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
or as a contiguous block of markers (`mbf` path). All methods
support multivariate (multi-environment) analysis through a `Trait`
argument that extends the model to capture genotype-by-trait
interaction.

Each pipeline follows the same three-stage workflow: **pre-analysis**
data preparation, **analysis**, and **post-analysis** interpretation.

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

## Vignettes

Full worked univariate and multivariate examples of each pipeline, from data preparation through to
interpretation and display.

| Vignette | Description |
|---|---|
| [Univariate QTL Mapping with wgAim](https://drj001.github.io/wgAim/articles/qtlAim_pipeline.html) | A complete walkthrough of the biparental QTL mapping pipeline: preparing a `wgCross` object with `primeCross()`, fitting a whole-genome composite model with `qtlAim()`, and interrogating the results with `summary()`, `plot()`, `aimTrace()`, `linkMap()`, and `fineMap()`. |
| [Multivariate QTL Mapping with wgAim](https://drj001.github.io/wgAim/articles/qtlAim_mv_pipeline.html) | Multi-environment QTL mapping using a factor analytic variance structure: base model development from `diag(Trial)` to `fa(Trial, 2)`, forward selection with `qtlAim()`, and visualisation of main-effect and G×E interaction QTL with `summary()`, `plot()`, and `aimTrace()`. |
| *GWAS with wgAim* (in preparation) | Panel preparation with `checkPanel()`, `filterPanel()`, and `primePanel()`; genome-wide association analysis with `gwasAim()`; and post-analysis display and fine mapping. |
| *Genomic prediction and selection with wgAim* (in preparation) | Genomic prediction with `gpAim()` and construction of weighted, Smith-Hazel, and desired-gains selection indices with `selIndex()`. |

---

## Part 1 — Pre-analysis: Data Preparation functions

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
wgCross <- primeCross(object, type = "interval", impute = "MartinezCurnow",
                      consensus.mark = TRUE, id = "id", subset = NULL,
                      infer = "mid")
```

| Argument | Description |
|---|---|
| `object` | A `cross` object produced by the **qtl** package (`"bc"`, `"dh"`, `"f2"`, or `"riself"`) |
| `type` | `"interval"` (default) — Haldane-weighted midpoint scores; `"marker"` — imputed marker genotypes only |
| `impute` | Imputation method for missing marker genotypes: `"MartinezCurnow"` (default) or `"Broman"` |
| `consensus.mark` | If `TRUE` (default), co-located markers are collapsed to consensus genotypes before imputation |
| `id` | Name of the line identifier column in `object$pheno`. Default `"id"` |
| `subset` | Character vector of chromosome names to retain; `NULL` (default) uses all chromosomes |
| `infer` | `"mid"` (default) places interval markers at midpoints; a numeric cM step size generates a finer grid. Only used when `type = "interval"` |

### Diversity panels — `checkPanel()`, `filterPanel()`, and `primePanel()`

Preparing a diversity panel for analysis follows a three-step workflow:
diagnose data quality with `checkPanel()`, apply filters with
`filterPanel()`, then build the analysis-ready object with `primePanel()`.
Each function accepts the output of the previous step directly, so the
typical call sequence is:

```r
chk   <- checkPanel(geno, map)
print(chk)
clean <- filterPanel(chk, miss.marker = 0.20, miss.line = 0.20, maf = 0.05)
print(clean)
panel <- primePanel(clean, impute = "knn", knn.k = 5)
```

#### `checkPanel()` — diagnose data quality

`checkPanel()` is a **pure diagnostic** tool. It inspects the raw
genotype matrix and map, reports everything it finds, and returns a
`"checkPanel"` object without modifying any data. The results are
intended to guide the choice of thresholds for `filterPanel()`.

| Check | What is reported |
|---|---|
| Encoding validation | Whether genotype values lie within the expected range for the stated encoding |
| Map consistency | Markers present in `geno` but absent from `map` (unplaceable), and vice versa |
| Marker missingness | Per-marker proportion of missing genotypes; count above a flagging threshold |
| Line missingness | Per-line proportion of missing genotypes; count above a flagging threshold |
| MAF distribution | Min / median / mean / max MAF; count of monomorphic markers; count removed at standard thresholds (0.01, 0.02, 0.05, 0.10) |
| Heterozygosity | Per-line and per-marker heterozygosity rates; count above a flagging threshold |
| Duplicate lines | Lines with identical genotype profiles across all common markers |
| Duplicate markers | Markers with identical genotype profiles across all lines |
| Chromosome coverage | Number of markers and cM range per chromosome |

The `print()` method accepts `miss.line.thresh`, `miss.marker.thresh`,
and `het.thresh` to control which lines and markers are flagged in the
summary output.

#### `filterPanel()` — apply data quality filters

`filterPanel()` applies a sequential set of filters in a statistically
principled order. Each step operates on data already cleaned by the
preceding steps. Passing a `"checkPanel"` object as the first argument
automatically reuses the stored `geno`, `map`, `encoding`, and column-name
arguments.

| Step | Filter | Default |
|---|---|---|
| 1 | **Map consistency** — markers absent from `map` are dropped | always applied |
| 2 | **Marker missingness** — markers above `miss.marker` are removed first, before lines are assessed | `0.20` |
| 3 | **Line missingness** — lines above `miss.line` are removed, computed on the cleaned marker set | `0.20` |
| 4 | **Line heterozygosity** — lines above `het.line` are removed; excess het indicates a mislabelled or contaminated sample | `NULL` (skipped) |
| 5 | **Marker heterozygosity** — markers above `het.marker` are removed; high per-marker het suggests a paralogous locus | `NULL` (skipped) |
| 6 | **Duplicate lines** — second and subsequent copies removed after quality filters, so clean copies are retained | `TRUE` |
| 7 | **Duplicate markers** — second and subsequent copies removed | `FALSE` |
| 8 | **MAF** — markers below `maf` are removed last, on the fully cleaned dataset | `0.05` |

The `print()` method reports each step, the number of lines or markers
removed, and the final panel dimensions.

#### `primePanel()` — build the analysis-ready object

`primePanel()` converts the filtered (and optionally imputed) genotype
data into a `wgPanel` object for use with `gwasAim()` or `gpAim()`.
It accepts a `"filteredPanel"` object directly; all column-name and
encoding arguments are then extracted automatically.

```r
wgPanel <- primePanel(geno, map, id = "id", map.id = "marker",
                      map.chr = "chr", map.pos = "pos",
                      encoding = "012", impute = "none", knn.k = 5L)
```

| Argument | Description |
|---|---|
| `geno` | Matrix (lines × markers) with row names identifying lines, a data frame with a line identifier column named by `id`, or a `"filteredPanel"` object from `filterPanel()` |
| `map` | Data frame containing the genetic map with columns for marker name, chromosome, and cM position. Not required when `geno` is a `"filteredPanel"` object |
| `id` | Name of the line identifier column when `geno` is a data frame; ignored when `geno` is a matrix. Default `"id"` |
| `map.id` | Column name in `map` holding marker names. Default `"marker"` |
| `map.chr` | Column name in `map` holding chromosome labels. Default `"chr"` |
| `map.pos` | Column name in `map` holding cM positions. Default `"pos"` |
| `encoding` | Genotype coding in `geno`: `"012"` (default) — allele counts 0/1/2, also accepts dosage values in [0, 2]; `"pm1"` — already in ±1 coding |
| `impute` | Handling of missing values: `"none"` (default) — stops with an informative error if any NAs remain; `"knn"` — chromosome-wise k-nearest-neighbour imputation (recommended in-package option); `"mean"` — column-mean imputation (a warning is always issued; not recommended above 1–2% missingness) |
| `knn.k` | Number of nearest neighbours for `impute = "knn"`. Default `5`. Ignored otherwise |

Both `primeCross()` and `primePanel()` recode genotypes internally to
the ±1 scale used by the analysis engine and store the prepared
matrices in `$imputed.data`.

---

## Part 2 — Analysis functions

Each analysis function takes a fitted **base ASReml-R model** and a prepared
genotype object, then progressively builds a whole-genome model using forward selection.

### QTL mapping — `qtlAim()`

Designed for biparental mapping populations (`wgCross`). Performs
forward selection of QTL as fixed effects against a composite
genome-wide random background, iterating until no further significant
QTL are detected.

```r
qtl.fit <- qtlAim(baseModel, genObj, merge.by, fix.lines = TRUE,
                  gen.type = NULL, method = "fixed",
                  selection = "interval", force = FALSE,
                  exclusion.window = 20, breakout = -1,
                  TypeI = 0.05, trace = TRUE, verboseLev = 0,
                  Trait = NULL, str = NULL, ...)
```

| Argument | Description |
|---|---|
| `baseModel` | A converged ASReml-R model capturing the experimental design; the genetic line term must appear in the random formula |
| `genObj` | A `wgCross` object from `primeCross()` |
| `merge.by` | Name of the line identifier column shared by the model data and `genObj` |
| `fix.lines` | If `TRUE` (default), ungenotyped lines are absorbed into a fixed `Gomit` factor; variance is estimated from genotyped lines only |
| `gen.type` | `"interval"` or `"marker"` — which genotype scores to use; inferred from the `primeCross()` type if `NULL` |
| `method` | `"fixed"` (default) — QTL enter as fixed effects; `"random"` — QTL enter as random effects (recommended when many QTL are expected) |
| `selection` | `"interval"` (default) — globally best interval; `"chromosome"` — best interval within the best chromosome |
| `force` | If `TRUE`, forces the `mbf` engine regardless of the markers-to-lines ratio |
| `exclusion.window` | cM exclusion zone around each detected QTL in subsequent iterations. Default `20` |
| `breakout` | Positive integer to stop after that many iterations without adding the final QTL; `-1` (default) runs to completion |
| `TypeI` | Family-wise significance threshold for the likelihood ratio test. Default `0.05` |
| `trace` | `TRUE` (default) prints ASReml output to console; a file path redirects it to that file |
| `verboseLev` | `0` (default, silent) or `1` — prints per-interval outlier statistics at each iteration |
| `Trait` | Column name of a trial/environment factor for multivariate G×E QTL analysis |
| `str` | Variance structure on the multivariate genomic term: `NULL` (default, mirrors base model), `"corh"`, `"corgh"`, `"us"`, `"diag"`, `"fa1"`, `"fa2"`, etc. |
| `...` | Additional arguments passed to `update.asreml()` at each iteration |

### GWAS — `gwasAim()`

Designed for diversity panels (`wgPanel`). Uses the same forward-
selection engine as `qtlAim()` but scans actual panel markers and
applies a panel-size-adjusted significance threshold.

```r
gwas.fit <- gwasAim(baseModel, genObj, merge.by, fix.lines = TRUE,
                    force = FALSE, exclusion.window = 20, breakout = -1,
                    TypeI = 0.05, trace = TRUE, verboseLev = 0,
                    Trait = NULL, str = NULL, ...)
```

| Argument | Description |
|---|---|
| `baseModel` | A converged ASReml-R model capturing the experimental design, including any population structure covariates |
| `genObj` | A `wgPanel` object from `primePanel()` |
| `merge.by` | Name of the line identifier column shared by the model data and `genObj` |
| `fix.lines` | If `TRUE` (default), ungenotyped lines are absorbed into a fixed `Gomit` factor |
| `force` | If `TRUE`, forces the `mbf` engine regardless of panel dimensions |
| `exclusion.window` | cM exclusion zone around each detected marker in subsequent iterations. Default `20` |
| `breakout` | Positive integer to stop after that many iterations; `-1` (default) runs to completion |
| `TypeI` | Significance threshold for the LRT. Default `0.05` |
| `trace` | `TRUE` (default) prints ASReml output to console; a file path redirects it to that file |
| `verboseLev` | `0` (default) or `1` — prints per-marker outlier statistics at each iteration |
| `Trait` | Column name of a trial/environment factor for multivariate G×E GWAS |
| `str` | Variance structure on the multivariate genomic term: `NULL` (default, mirrors base model), `"corh"`, `"corgh"`, `"us"`, `"diag"`, `"fa1"`, `"fa2"`, etc. |
| `...` | Additional arguments passed to `update.asreml()` at each iteration |

Note: `method = "fixed"` and `selection = "interval"` are fixed internally.

### Genomic prediction — `gpAim()`

Fits a whole-genome G-BLUP model and extracts genomic estimated
breeding values (GEBVs), prediction accuracies, and generalised
heritabilities for all genotyped lines. Works with both `wgCross`
and `wgPanel` objects.

```r
gp.fit <- gpAim(baseModel, genObj, merge.by, fix.lines = TRUE,
                gen.type = "marker", force = FALSE,
                trace = TRUE, Trait = NULL, str = NULL, ...)
```

| Argument | Description |
|---|---|
| `baseModel` | A converged ASReml-R model capturing the experimental design; no genomic term should be present |
| `genObj` | A `wgCross` object from `primeCross()` or a `wgPanel` object from `primePanel()` |
| `merge.by` | Name of the line identifier column shared by the model data and `genObj` |
| `fix.lines` | If `TRUE` (default), ungenotyped lines are absorbed into a fixed `Gomit` factor |
| `gen.type` | `"marker"` (default) — uses `$imputed.data`; `"interval"` — uses `$interval.data` (requires `wgCross`) |
| `force` | If `TRUE`, forces the `mbf` engine regardless of the markers-to-lines ratio |
| `trace` | `TRUE` (default) prints ASReml output to console; a file path redirects it to that file |
| `Trait` | Column name of a trial/environment factor for multivariate genomic prediction |
| `str` | Variance structure on the multivariate genomic term: `NULL` (default, mirrors base model), `"corh"`, `"corgh"`, `"us"`, `"diag"`, `"fa1"`, `"fa2"`, etc. |
| `...` | Additional arguments passed to `update.asreml()` |

Key outputs stored in `$GP`:

| Slot | Content |
|---|---|
| `$gebv` | Data frame: line, GEBV, SE, accuracy (Mrode formula), generalised H² (Cullis) |
| `$var.genetic` | Estimated additive genetic variance |
| `$rel.matrix` | Genomic relationship matrix (vm path) |

### Selection index — `selIndex()`

Combines per-environment GEBVs from a multivariate `gpAim()` result
into a single selection index. Three index types are available:

```r
si <- selIndex(gp.fit, weights = c(E1 = 1, E2 = 2),
               type = "smith-hazel", prop.select = 0.10, selection = NULL)
```

| `type` | Method |
|---|---|
| `"weighted"` | Weighted sum of GEBVs |
| `"smith-hazel"` | Classical Smith-Hazel index |
| `"desired-gains"` | Pesek-Baker desired-gains index |

Users can use the various method to help select the top `prop.select` of
lines. A `selection` of lines can also be passed to the function. in both
instamces the `$gain` output reports the expected genetic gain per trait (or
environment) for the selection of lines.

---

## Part 3 — Post-analysis: Interpretation and Display functions

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

### Iteration diagnostics — `aimTrace()`

`aimTrace()` traces the forward-selection iterations of a `qtlAim` or
`gwasAim` fit, showing how QTL were progressively added and how their
effects evolved.

```r
aimTrace(qtl.fit, plot = "lrt")       # LRT statistic at each iteration
aimTrace(qtl.fit, plot = "stability") # QTL effect ± 1 SE across iterations
aimTrace(qtl.fit, plot = "both")      # returns both plots as a list
```

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
