# =============================================================================
# theme_scatter.R
# Shared ggplot2 theme used by plot.qtlAim, plot.gwasAim, plot.gpAim,
# and manhattan.gwasAim.
# =============================================================================

#' @keywords internal
theme_scatter <- function(base_size = 11, base_family = "") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
        legend.position      = "none",
        panel.grid.minor     = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major     = element_line(colour = "grey90", linewidth = 0.8),
        panel.grid.minor.x   = element_blank(),
        panel.grid.major.x   = element_blank(),
        panel.background     = element_blank(),
        axis.text.x          = element_text(size = base_size),
        axis.ticks.x         = element_blank(),
        axis.text.y          = element_text(size = base_size),
        strip.text           = element_text(size = base_size)
    )
}
