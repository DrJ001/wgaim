#' Extract QTL position information from a fitted model
#'
#' Extracts chromosome, flanking marker names, and cM positions for each
#' significant QTL (or marker) from a fitted \code{qtlAim} or \code{gwasAim}
#' object.
#'
#' @param object A fitted object of class \code{"qtlAim"} or \code{"gwasAim"}.
#' @param genObj The \code{wgCross} or \code{wgPanel} object used in the
#'   analysis, produced by \code{primeCross()} or \code{primePanel()}.
#'
#' @return A character matrix.  For marker-type QTL (4 columns):
#'   chromosome, interval index, marker name, marker position (cM).
#'   For interval-type QTL (8 columns): chromosome, interval index,
#'   left-flanking marker name and position, interval midpoint name and
#'   position, right-flanking marker name and position.
#'
#' @export
getQTL <- function(object, genObj)
{
    spe  <- strsplit(substring(unlist(names(object$QTL$effects)), 3), "\\.")
    wchr <- unlist(lapply(spe, function(el) el[1]))
    wint <- as.numeric(unlist(lapply(spe, function(el) el[2])))

    if (object$QTL$type == "marker") {
        qtlm <- matrix(ncol = 4, nrow = length(wchr))
        for (i in seq_along(wchr)) {
            lhmark <- genObj$geno[[wchr[i]]]$map[wint[i]]
            qtlm[i, 1:4] <- c(wchr[i], wint[i], names(lhmark), round(lhmark, 2))
        }
    } else {
        qtlm <- matrix(ncol = 8, nrow = length(wchr))
        for (i in seq_along(wchr)) {
            if (length(genObj$geno[[wchr[i]]]$inferred.map) > 1) {
                imark  <- genObj$geno[[wchr[i]]]$inferred.map[wint[i]]
                lhint  <- findInterval(imark, genObj$geno[[wchr[i]]]$map)
                lhmark <- genObj$geno[[wchr[i]]]$map[lhint]
                rhmark <- genObj$geno[[wchr[i]]]$map[lhint + 1]
            } else {
                lhmark <- rhmark <- imark <- genObj$geno[[wchr[i]]]$map[wint[i]]
            }
            qtlm[i, 1:8] <- c(wchr[i], wint[i],
                               names(lhmark), round(lhmark, 2),
                               names(imark),  round(imark,  2),
                               names(rhmark), round(rhmark, 2))
        }
    }
    qtlm
}
