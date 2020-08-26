
### wgaim: An R package for efficient whole genome QTL analysis

Authors: Julian Taylor & Ari Verbyla

This is the public facing GitHub repository verison of the R package [wgaim](https://cran.r-project.org/package=wgaim) available on CRAN.

The package was built to implement the original wgaim algorithm in Verbyla et al. (2007, 2012). wgaim is a whole genome average interval mapping R package that uses ASReml-R V4 for it core linear mixed modelling functionality. To use full functionality of the package users will require a valid license for ASReml-R V4 and this can be obtained from \url{https://www.vsni.co.uk/software/asreml-r}. 

To install the package from GitHub you will need to do the following: 

1. You need to install the [devtools](https://cran.r-project.org/package=devtools) package. You can do this by invoking R and then typing

```
install.packages("devtools")
```

2. Load the devtools package.

```
library(devtools)
```

3. Install wgaim. In most cases, you just use install_github("author/package"). 

```
install_github("DrJ001/wgaim")
```

### Getting Started

For a quick but complete introduction of the functionality of the package please visit the [vignette](https://cran.r-project.org/web/packages/wgaim/vignettes/wgaim_intro.html) on the CRAN package page.



