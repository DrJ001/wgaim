\name{news}
\title{Release information for wgaim}
\section{Changes in version 2.0-6}{
  \itemize{
    \item Version 2.0-6 of \pkg{wgaim} is a minor maintenance release
    to fix CRAN build errors.
    }
}
\section{Changes in version 2.0-x}{
\subsection{NEW FEATURES}{
  \itemize{
    \item Version 2.0-x of \pkg{wgaim} is here and it celebrates
    a updated release of the package that utilizes the linear mixed
    modelling functionality of the R package \code{ASReml-R} V4. It
    should be noted that this version of \code{wgaim} is not compatible
    with \code{ASReml-R} V3 and users should revert to version 1.4-11 if
    a compatible version is required. Within this new version of
    \code{wgaim} there have been many subtle changes to functions and
    their arguments. Most of these changes have been documented below.
    \item The code in the main calling function \code{wgaim.asreml()} has been
    significantly streamlined for better integration with new features
    of \code{ASReml-R} V4. Many adjunct \code{wgaim} functions such as
    \code{mergeData()} and \code{updateWgaim()} have been removed. Model
    updating now occurs directly using \code{update.asreml()}.
    \item To ensure consistency between the phenotypic data used in
    fitting the base model and the phenotypic data used in \code{wgaim.asreml()},
    the \code{phenoData} argument has been removed from the
    \code{wgaim.asreml} call. The data is now recalled through
    backwards evaluation of the base model call.
    \item The \pkg{wgaim} can now handle \code{"f2"} cross objects. This
    includes the appropriate imputation of missing allele values through the
    functionality of \code{cross2int()}.
    \item Some of the \code{cross2int()} functions arguments have been
    changed to more appropriately reflect the nature of the task being
    implemented. Specifically, argument \code{missgeno} has been changed to
    \code{impute} and \code{rem.mark} has been changed to
    \code{consensus.mark}.
    \item All \code{link.map.xxx} functions have changed to
    \code{linkMap.xxx} for better naming consistency with S3 methods.
    \item The \code{out.stat} function has changed to \code{outStat} and has
    been completely rewritten to use \pkg{ggplot2} functionality. See
    \code{?outStat} for complete details.
    \item A vignette title "An Quick Introudction to wgaim QTL Analysis"
    is available with the package and can be viewed using:

    \code{vignette("wgaim_intro")}
  }
}
\subsection{BUG FIXES}{
  \itemize{
    \item Fixed a chromosome labelling issue with \code{outStat()}. The
    fix now ensures any graphic generated with the function produces the
    correct chromosome labels in the appropriate position.
    \item Fixed a bug that caused the algorithm to crash when there is
    only one marker on a chromosome and the argument \code{gen.type =
    "interval"} is used in \code{wgaim.asreml()}.
    \item Fixed a bug that caused the algroithm to crash when
    the number of markers in the model changes from being greater than
    the number of lines to less than the number of lines (after
    selection and exclusion).
    \item Fixed a bug that caused \code{summary.wgaim()} when only one
    QTL was found.
    \item Removed \code{maxiter = 1} from internal
    \code{predict.asreml()} to prevent spurious output of
    non-convergence warnings.
    \item The use of "." in the chromosome names muddles regular
    expression string matching in parts of the wgaim call. There is now
    code in \code{cross2int()} to remove "."s from any chromosome names.
    \item Flanking marker highlighting in \code{linkMap.wgaim()} was
    incorrect and has been amended.
}
}
}
\section{Changes in version 1.4-x}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Restrictions on the layout of the plots in \code{out.stat()}
	have been removed.
      \item The argument \code{flanking} has been added to the QTL
	plotting functions ro ensure that only flanking markers or
	linked markers are plotted and highlighted on the linkage map.
      \item The forward selection algorithm has been accelerated further
      by smart matrix decomposition of the relationship matrix. Users
      can expect around a 35\% reduction in computation time.
      \item Outlier statistics and BLUPs can now be returned for any
      iteration of the algorithm regardless of whether a
      significant QTL is detected. This now allows easy access to
      outlier statistics and BLUPs for the first iteration when no QTL
      are detected (see the \code{breakout} argument of \code{wgaim.asreml()}.
      \item The package now includes a PDF reference manual that is accessible by
      navigating to the \code{"doc"} directory of the package. This can
      be found on any operating system using the command

      \code{> system.file("doc", package = "wgaim")}

      The reference manual contains WGAIM theory and two thorough examples that show the
      features of the package. It also contains a "casual walk through"
      the package providing the user with a series of 5 steps to a successful wgaim
      analysis.
      \item The package now includes three fully documented phenotypic
      and genotypic data sets for users to explore. Two of these three
      have been used in the manual and scripts that follow the
      examples in the manual are available under the "doc" directory of the package.
      \item The package now provides very efficient whole genome QTL analysis of high
      dimensional genetic marker data. All genetic marker data is passed
      into \code{wgaim.asreml()} through the \code{"intervalObj"}
      argument. Merging of genotypic and phenotypic data occurs within
      \code{wgaim.asreml()}.
      \item \code{wgaim.asreml()} has several new arguments related
      to selection of QTL. The \code{"gen.type"} argument allows the user to
  choose a whole genome marker analysis or whole genome mid-point
  interval analysis from Verbyla et. al (2007). The \code{"method"} argument gives you the choice of placing
  QTL in the fixed part of the linear mixed model as in Verbyla et.al
  (2007) or the random part of model as in Verbyla et. al
  (2012). Finally, the \code{"selection"} argument allows you to choose whether QTL selection
  is based on whole genome interval outlier statistics or a two stage process of
  using chromosome outlier statistics and then interval outlier
  statistics.
  \item A \code{"breakout"} argument is now also provided which allows
  the user to breakout of the forward selection algorithm at any
  stage. The current model along with any calculated QTL components are
  all available for inspection.
  \item All linkage map plotting functions can be subsetted by
  predefined distances. This includes a list of distances as long as the
  number of linkage groups plotted.
  }
}
\subsection{BUG FIXES}{
  \itemize{
    \item Fixed a bug that created 0s for NAs in linkage groups with one
      marker.
    \item Fixed a bug that caused \code{wgaim.asreml()} to bomb out if
      the number of markers was less than the number of genotypes.
    \item Fixed a bug that outputted warning messages regarding a
    \code{NaN} calculation from \code{sqrt(vatilde)} in \code{qtl.pick()}.
    \item Fixed a bug that caused wgaim to crash if \code{method =
      "random"} was used with the new version of asreml.
    \item Fixed a bug that caused wgaim to crash with very recent
      versions of asreml (04/05/2015).
    \item \code{cross2int()} now accepts R/qtl objects with cross type
    \code{"bc","dh","riself"}.
    \item Fixed an issue with the internal function \code{fix.map()} that
    allowed some co-located sets of markers to appear in the final
    reduced linkage map.
    \item Fixed a long standing scoping issue with different versions of ASReml-R.
    \item Fixed an elusive problem that causes wgaim models to increase the size
    of your .RData upon saving. This is actually an inherent problem with
    using model formula within a function a returning the subsequent
    model. There is now a function at the very tail of
    \code{wgaim.asreml()} that quietly destroys the useless environments
    that these formula contain.
    \item Fixed bug that caused \code{wgaim.asreml()} to crash when no QTL
    were found.
    \item Fixed bug that caused \code{summary.wgaim()} to crash when one
    QTL was found using \code{method = "random"}.
  }
}
}


