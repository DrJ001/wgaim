# wgAim — Whole Genome Analyses via Integrated Modelling

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![Package version](https://img.shields.io/badge/Package%20version-0.1--0-orange.svg?style=flat-square)](https://github.com/DrJ001/wgAim/commits/dev)
[![codecov](https://codecov.io/gh/DrJ001/wgAim/branch/dev/graph/badge.svg)](https://app.codecov.io/gh/DrJ001/wgAim)
<!-- badges: end -->

**Authors**: Julian Taylor & Ari Verbyla

> ⚠️ **This package is under active development.** The API is not yet
> stable and breaking changes may occur between commits. Check back
> regularly for updates — a first release will be announced here when
> the package is ready for general use.

---

## Overview

**wgAim** is a unified whole-genome analysis toolkit for plant breeding
populations, built entirely around the
[ASReml-R](https://vsni.co.uk/software/asreml-r/) linear mixed modelling
engine. It extends and supersedes the
[wgaim](https://github.com/DrJ001/wgaim-legacy) package, broadening
the scope from QTL detection alone into a three-pillar framework for
whole-genome inference:

| Pillar | Function | Description |
|---|---|---|
| QTL detection | `qtlAim()` | Forward-selection whole-genome marker or interval mapping in biparental populations |
| GWAS | `gwasAim()` | Forward-selection genome-wide association analysis in diversity panels |
| Genomic prediction | `gpAim()` | Genomic best linear unbiased prediction (G-BLUP) delivering GEBVs |

All three analyses share a common internal engine and the same
unifying statistical idea: the full complement of markers or intervals
enters the mixed model simultaneously as a composite genome-wide random
effect, represented either as a genomic relationship matrix (`vm` path)
or as a direct marker-by-file random effect (`mbf` path).

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
