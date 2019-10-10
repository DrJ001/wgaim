\name{wgaim-internal}
\alias{constructCM}
\alias{fixMap}
\alias{envFix}
\alias{getQTL}
\alias{addiag}
\alias{imputeGen}
\alias{wgaim.default}
\alias{qtlSelect}
\alias{vModify}
\alias{theme_scatter}
\title{Internal wgaim functions}
\description{
These are internal \code{wgaim} functions that are not intended to be
directly called by the user.
}
\usage{
constructCM(genoData, scale.method = "diag")
fixMap(full.data, rd = 3)
envFix(model, asremlEnv)
getQTL(object, intervalObj)
addiag(x = 1, di = 0, nrow.arg, ncol.arg = n)
imputeGen(theta, chr, dom = TRUE)
\method{wgaim}{default}(baseModel, \ldots)
qtlSelect(asm, phenoData, intervalObj, gen.type, selection,
          exclusion.window, state, verboseLev)
vModify(model, merge.by)
theme_scatter(base_size = 11, base_family = "")
}
\author{Julian Taylor and Ari Verbyla}
\keyword{internal}
