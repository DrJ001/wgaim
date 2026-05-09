# =============================================================================
# linkMap.QTLAim.R
# S3 linkMap() method for QTLAim objects.
# Overlays detected QTL regions on the genetic linkage map.
# Note: the linkMap() generic and linkMap.cross() are defined in wgaim16.R.
# =============================================================================

#' @describeIn QTLAim Plot the genetic linkage map with detected QTL overlaid
#'   as coloured polygons (interval mode) or points (marker mode). QTL regions
#'   are annotated with flanking marker names and connected to trait labels by
#'   leader lines. Multiple traits can be overlaid using different colours by
#'   passing a list of \code{QTLAim} objects to \code{\link{linkMap.default}}.
#' @param object A \code{QTLAim} object.
#' @param intervalObj The \code{"cross"} or \code{"interval"} object used in
#'   the analysis.
#' @param chr Optional character vector of chromosome names to display. Default
#'   is all chromosomes carrying at least one detected QTL.
#' @param chr.dist Optional named list with elements \code{$start} and/or
#'   \code{$end} (numeric vectors of cM positions) to restrict the displayed
#'   map range per chromosome.
#' @param marker.names Character string controlling marker label display:
#'   \code{"markers"} (default, show marker names), \code{"dist"} (show cM
#'   positions), or \code{NULL} (no marker labels, compact display).
#' @param flanking Logical. If \code{TRUE} (default), only the flanking markers
#'   of detected QTL intervals are annotated, reducing clutter.
#' @param list.col Named list with colour specifications: \code{$q.col} (QTL
#'   fill colour, default \code{"light blue"}), \code{$m.col} (flanking marker
#'   label colour, default \code{"red"}), \code{$t.col} (trait label colour,
#'   default same as \code{q.col}).
#' @param list.cex Named list with text size specifications: \code{$t.cex}
#'   (trait label size, default \code{0.6}), \code{$m.cex} (marker label size,
#'   default \code{0.6}).
#' @param trait.labels Optional character vector of trait label(s) to display.
#'   Defaults to the response variable name from the model formula.
#' @param tick Logical. If \code{TRUE}, tick marks are drawn on the chromosome
#'   axis. Default \code{FALSE}.
#' @export
linkMap.QTLAim <- function (object, intervalObj, chr, chr.dist,
                            marker.names = "markers", flanking = TRUE,
                            list.col = list(q.col = "light blue",
                                            m.col = "red",
                                            t.col = "light blue"),
                            list.cex = list(t.cex = 0.6, m.cex = 0.6),
                            trait.labels = NULL, tick = FALSE, ...) {
    dots <- list(...)
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj must be of class \"cross\".")
    if (!length(object$QTL$effects)) {
        warning("There are no significant QTL's. Plotting map only...")
        linkMap(intervalObj, chr, chr.dist, marker.names = marker.names,
                tick = tick, squash = FALSE, ...)
        return(invisible())
    }
    qtlm <- getQTL(object, intervalObj)
    wchr <- qtlm[, 1]
    if (object$QTL$type == "interval")
        qtlm <- qtlm[, 3:6, drop = FALSE]
    else
        qtlm <- cbind(qtlm[, 3:4, drop = FALSE], qtlm[, 3:4, drop = FALSE])
    if (is.null(list.col$q.col)) list.col$q.col <- "light blue"
    if (is.null(list.col$t.col)) list.col$t.col <- list.col$q.col
    if (missing(chr))
        chr <- unique(wchr)[order(unique(wchr))]
    if (is.null(list.cex$m.cex)) list.cex$m.cex <- 0.6
    if (is.null(list.cex$t.cex)) list.cex$t.cex <- 0.6
    if (flanking)
        attr(intervalObj, "flanking") <- unique(c(qtlm[, 1], qtlm[, 3]))
    lmap <- linkMap(intervalObj, chr, chr.dist, marker.names = marker.names,
                    tick = tick, squash = TRUE, m.cex = list.cex$m.cex, ...)
    map <- lmap$map
    if (is.null(trait <- trait.labels))
        trait <- rep(as.character(object$call$fixed[[2]]), length(wchr))
    if (length(trait.labels) == 1)
        trait <- rep(trait.labels, length(wchr))
    qtrait <- unique(trait)
    qtlm   <- cbind.data.frame(qtlm, trait = factor(trait, levels = unique(trait)))
    qtlm[, 2] <- as.numeric(as.character(qtlm[, 2]))
    qtlm[, 4] <- as.numeric(as.character(qtlm[, 4]))
    qtlList <- lapply(split(qtlm, wchr), function(el) el[order(el[, 2]), ])
    qtlm    <- do.call("rbind", qtlList)
    wchr    <- wchr[order(wchr)]
    chr     <- names(map)
    if (!missing(chr)) {
        oc <- wchr %in% chr
        if (!all(oc)) {
            warning("Some QTL's exist outside chromosome(s) subset. Omitting QTL's...")
            qtlm <- qtlm[oc, ]
            wchr <- wchr[oc]
        }
    }
    if (!missing(chr.dist)) {
        mins <- unlist(lapply(map, min)); maxs <- unlist(lapply(map, max))
        mins <- mins[wchr]; maxs <- maxs[wchr]
        om <- qtlm[, 4] < maxs & qtlm[, 2] > mins
        if (!all(om)) {
            warning("Some QTL regions outside distances specified. Omitting QTL's...")
            qtlm <- qtlm[om, ]; wchr <- wchr[om]
        }
    }
    n.chr  <- length(map)
    maxlen <- max(unlist(lapply(map, max)))
    chrpos <- lmap$chrpos
    mt     <- lmap$mt
    if (!is.na(cind <- pmatch("col", names(dots)))) dots <- dots[-cind]
    if (is.null(dim(qtlm))) qtlm <- matrix(qtlm, nrow = 1, byrow = FALSE)
    tlis <- list()
    p.cex <- if (!is.na(pmatch("cex", names(dots)))) dots$cex else par("cex")
    dots$cex <- NULL
    for (i in 1:n.chr) {
        conv <- par("pin")[2] / maxlen
        ind  <- wchr %in% names(map)[i]
        if (any(ind)) {
            tlis[[i]] <- as.vector(qtlm[ind, 2] + qtlm[ind, 4]) / 2
            names(tlis[[i]]) <- as.character(qtlm[ind, 5])
            if (length(tlis[[i]]) > 1) {
                for (j in 1:(length(tlis[[i]]) - 1)) {
                    ch <- tlis[[i]][j + 1] * conv - (tlis[[i]][j] * conv +
                          10 * par("csi") * list.cex$t.cex / 9)
                    if (ch < 0) tlis[[i]][j + 1] <- (tlis[[i]][j + 1] * conv + abs(ch)) / conv
                }
            }
        }
    }
    qtld  <- qtlm[, 1:4]
    nodup <- !duplicated(do.call("paste", as.data.frame(qtld)))
    qtls  <- qtld[nodup, , drop = FALSE]
    whd   <- pmatch(do.call("paste", as.data.frame(qtld)),
                    do.call("paste", as.data.frame(qtls)), duplicates.ok = TRUE)
    dlis  <- split(as.character(qtlm[, 5]), whd)
    qtlm  <- qtls
    wchr  <- wchr[nodup]
    for (i in 1:n.chr) {
        ind <- wchr %in% names(map)[i]
        if (any(ind)) {
            ind <- (1:length(wchr))[ind]
            for (j in ind) {
                if (!is.null(marker.names)) {
                    wh <- mt[[i]][c(as.character(qtlm[j, 1]), as.character(qtlm[j, 3]))]
                    if (marker.names == "dist") {
                        dist  <- map[[i]][c(as.character(qtlm[j, 1]), as.character(qtlm[j, 3]))]
                        dlabs <- as.character(round(as.numeric(unlist(dist)), 2))
                        alis  <- list(x = chrpos[i] + 0.5, y = wh, labels = dlabs,
                                      adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                    } else {
                        alis <- list(x = chrpos[i] + 0.5, y = wh, labels = names(wh),
                                     adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                    }
                    do.call("text", c(alis, dots))
                }
                yv    <- c(qtlm[j, 2], qtlm[j, 4]); yv <- c(yv, rev(yv))
                dind  <- dlis[[j]]
                q.cols <- list.col$q.col[pmatch(dind, qtrait)]
                qind  <- 1:length(dind)
                if (length(dlis[[j]]) > 1) {
                    int <- seq(chrpos[i] - 0.2, chrpos[i] + 0.2, length = length(dind) + 1)
                    for (k in 1:length(dind)) {
                        xv <- c(rep(int[k], 2), rep(int[k + 1], 2))
                        if (object$QTL$type == "interval")
                            polygon(xv, y = yv, border = NA, col = q.cols[k])
                        else {
                            plis <- list(x = (xv[1] + xv[3]) / 2, y = yv[1],
                                         col = q.cols[k], cex = p.cex)
                            do.call("points", c(plis, dots))
                        }
                    }
                } else {
                    xv <- c(rep(chrpos[i] - 0.2, 2), rep(chrpos[i] + 0.2, 2))
                    if (object$QTL$type == "interval")
                        polygon(xv, y = yv, border = NA, col = q.cols)
                    else {
                        plis <- list(x = chrpos[i], y = yv[1], col = q.cols, cex = p.cex)
                        do.call("points", c(plis, dots))
                    }
                }
                midpt <- sum(yv[1:2]) / 2
                segments(chrpos[i] - 0.25, yv[1], chrpos[i] - 0.25, yv[2])
                segments(chrpos[i] - 0.25, midpt, chrpos[i] - 0.3, midpt)
                segments(chrpos[i] - 0.3, midpt, chrpos[i] - 0.4, tlis[[i]][qind])
                segments(chrpos[i] - 0.4, tlis[[i]][qind], chrpos[i] - 0.45, tlis[[i]][qind])
                t.cols <- if (length(list.col$t.col) > 1)
                              list.col$t.col[pmatch(dind, qtrait)]
                          else list.col$t.col
                text(chrpos[i] - 0.5, tlis[[i]][qind], names(tlis[[i]][qind]),
                     adj = c(1, 0.3), col = t.cols, cex = list.cex$t.cex)
                tlis[[i]] <- tlis[[i]][-qind]
            }
        }
        segments(chrpos[i] - 0.2, map[[i]], chrpos[i] + 0.2, map[[i]])
    }
    if (is.na(pmatch("main", names(dots))))
        title("Genetic Map with QTL")
}
