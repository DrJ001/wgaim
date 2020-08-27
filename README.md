---
output:
  html_document:
    keep_md: yes
---

<!-- README.md is generated from README.Rmd. Please edit that file -->





### R/wgaim: An R package for efficient whole genome QTL analysis

<!-- badges: start -->

[![Project Status: Active:  The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.0.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/wgaim)](https://cran.r-project.org/package=wgaim) 
[![packageversion](https://img.shields.io/badge/Package%20version-2.0--5-orange.svg?style=flat-square)](/commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--08--27-yellowgreen.svg)](/commits/master)
<!-- badges: end -->

**Authors**: Julian Taylor & Ari Verbyla

This is the public facing GitHub repository verison of the R package [wgaim](https://cran.r-project.org/package=wgaim) available on CRAN. This repository is updated more regularly than the CRAN version of the package and includes minor additions and small bug fixes. Please check the badges at the top of this README and the NEWS for the package. 

The package was built to implement the original wgaim algorithm in Verbyla et al. (2007, 2012). **R/wgaim** is a whole genome average interval mapping R package that uses ASReml-R V4 for it core linear mixed modelling functionality. To use full functionality of the package users will require a valid license for ASReml-R V4 and this can be obtained from [https://www.vsni.co.uk/software/asreml-r](https://www.vsni.co.uk/software/asreml-r). 

To install the package from GitHub you will need to do the following: 

1. Install the [devtools](https://cran.r-project.org/package=devtools) package. Do this by invoking R and then typing


```r
install.packages("devtools")
```

2. Install wgaim using 


```r
devtools::install_github("DrJ001/wgaim")
```

#### Getting Started

For a quick but complete introduction of the functionality of the package please visit the [vignette](https://cran.r-project.org/web/packages/wgaim/vignettes/wgaim_intro.html) on the CRAN package page.

#### References

Verbyla, A.P., Cullis, B...R. & Thompson, R. (2007) The analysis of QTL by simultaneous use of the of the full linkage map. *Theoretical and Applied Genetics*, **116**, 95-111.

Verbyla, A.P., Taylor, J.D. & Verbyla, K.L. (2012) RWGAIM: An efficient high dimensional random whole genome average (QTL) interval mapping approach. *Genetics Research*, **94**, 291-306.
