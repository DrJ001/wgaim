# Release information for wgaim

## Changes in version 2.0-6

Version 2.0-6 of `wgaim` is a minor maintenance release to fix CRAN build warnings.

## Changes in version 2.0-x
### New Features

- Version 2.0-x of `wgaim` is here and it celebrates
a updated release of the package that utilizes the linear mixed
modelling functionality of the R package `ASReml-R` V4. It
should be noted that this version of `wgaim` is not compatible
with `ASReml-R` V3 and users should revert to version 1.4-11 if
a compatible version is required. Within this new version of
`wgaim` there have been many subtle changes to functions and
their arguments. Most of these changes have been documented below.
- The code in the main calling function `wgaim.asreml()` has been
significantly streamlined for better integration with new features
of `ASReml-R` V4. Many adjunct `wgaim` functions such as
`mergeData()` and `updateWgaim()` have been removed. Model
updating now occurs directly using `update.asreml()`.
- To ensure consistency between the phenotypic data used in
fitting the base model and the phenotypic data used in `wgaim.asreml()`,
the `phenoData` argument has been removed from the
`wgaim.asreml` call. The data is now recalled through
backwards evaluation of the base model call.
- The `wgaim` can now handle `"f2"` cross objects. This
includes the appropriate imputation of missing allele values through the
functionality of `cross2int()`.
- Some of the `cross2int()` functions arguments have been
changed to more appropriately reflect the nature of the task being
implemented. Specifically, argument `missgeno` has been changed to
`impute` and `rem.mark` has been changed to
`consensus.mark`.
- All `link.map.xxx` functions have changed to
`linkMap.xxx` for better naming consistency with S3 methods.
- The `out.stat` function has changed to `outStat` and has
been completely rewritten to use `ggplot2` functionality. See
`?outStat` for complete details.
- A vignette title "An Quick Introudction to wgaim QTL Analysis"
is available with the package and can be viewed using: `vignette("wgaim_intro")`

### Bug Fixes

- Fixed a chromosome labelling issue with `outStat()`. The
fix now ensures any graphic generated with the function produces the
correct chromosome labels in the appropriate position.
- Fixed a bug that caused the algorithm to crash when there is
only one marker on a chromosome and the argument `gen.type =
"interval"` is used in `wgaim.asreml()`.
- Fixed a bug that caused the algroithm to crash when
the number of markers in the model changes from being greater than
the number of lines to less than the number of lines (after
selection and exclusion).
- Fixed a bug that caused `summary.wgaim()` when only one
QTL was found.
- Removed `maxiter = 1` from internal
`predict.asreml()` to prevent spurious output of
non-convergence warnings.
- The use of "." in the chromosome names muddles regular
expression string matching in parts of the wgaim call. There is now
code in `cross2int()` to remove "."s from any chromosome names.
- Flanking marker highlighting in `linkMap.wgaim()` was
incorrect and has been amended.

## Changes in version 1.4-x
### New Features

- Restrictions on the layout of the plots in `out.stat()`
have been removed.
- The argument `flanking` has been added to the QTL
plotting functions ro ensure that only flanking markers or
linked markers are plotted and highlighted on the linkage map.
- The forward selection algorithm has been accelerated further
by smart matrix decomposition of the relationship matrix. Users
can expect around a 35\% reduction in computation time.
- Outlier statistics and BLUPs can now be returned for any
iteration of the algorithm regardless of whether a
significant QTL is detected. This now allows easy access to
outlier statistics and BLUPs for the first iteration when no QTL
are detected (see the `breakout` argument of `wgaim.asreml()`.
- The package now includes a PDF reference manual that is accessible by
navigating to the `"doc"` directory of the package. This can
be found on any operating system using the command

`> system.file("doc", package = "wgaim")`

The reference manual contains WGAIM theory and two thorough examples that show the
features of the package. It also contains a "casual walk through"
the package providing the user with a series of 5 steps to a successful wgaim
analysis.
- The package now includes three fully documented phenotypic
and genotypic data sets for users to explore. Two of these three
have been used in the manual and scripts that follow the
examples in the manual are available under the "doc" directory of the package.
- The package now provides very efficient whole genome QTL analysis of high
dimensional genetic marker data. All genetic marker data is passed
into `wgaim.asreml()` through the `"intervalObj"`
argument. Merging of genotypic and phenotypic data occurs within
`wgaim.asreml()`.
- `wgaim.asreml()` has several new arguments related
to selection of QTL. The `"gen.type"` argument allows the user to
choose a whole genome marker analysis or whole genome mid-point
interval analysis from Verbyla et. al (2007). The `"method"` argument gives you the choice of placing
QTL in the fixed part of the linear mixed model as in Verbyla et.al
(2007) or the random part of model as in Verbyla et. al
(2012). Finally, the `"selection"` argument allows you to choose whether QTL selection
is based on whole genome interval outlier statistics or a two stage process of
using chromosome outlier statistics and then interval outlier
statistics.
- A `"breakout"` argument is now also provided which allows
the user to breakout of the forward selection algorithm at any
stage. The current model along with any calculated QTL components are
all available for inspection.
- All linkage map plotting functions can be subsetted by
predefined distances. This includes a list of distances as long as the
number of linkage groups plotted.

### Bug Fixes

- Fixed a bug that created 0s for NAs in linkage groups with one
marker.
- Fixed a bug that caused `wgaim.asreml()` to bomb out if
the number of markers was less than the number of genotypes.
- Fixed a bug that outputted warning messages regarding a
`NaN` calculation from `sqrt(vatilde)` in `qtl.pick()`.
- Fixed a bug that caused wgaim to crash if `method =
"random"` was used with the new version of asreml.
- Fixed a bug that caused wgaim to crash with very recent
versions of asreml (04/05/2015).
- `cross2int()` now accepts R/qtl objects with cross type
`"bc","dh","riself"`.
- Fixed an issue with the internal function `fix.map()` that
allowed some co-located sets of markers to appear in the final
reduced linkage map.
- Fixed a long standing scoping issue with different versions of ASReml-R.
- Fixed an elusive problem that causes wgaim models to increase the size
of your .RData upon saving. This is actually an inherent problem with
using model formula within a function a returning the subsequent
model. There is now a function at the very tail of
`wgaim.asreml()` that quietly destroys the useless environments
that these formula contain.
- Fixed bug that caused `wgaim.asreml()` to crash when no QTL
were found.
- Fixed bug that caused `summary.wgaim()` to crash when one
QTL was found using `method = "random"`.
