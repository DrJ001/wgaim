# =============================================================================
# build_vignette.R
#
# Renders the wgAim vignette(s) from source and copies the pre-built outputs
# into inst/doc/ ready for CRAN submission.
#
# Usage (from the package root):
#   Rscript build_vignette.R
# or source it inside an R session:
#   source("build_vignette.R")
#
# Requirements:
#   - Quarto installed and on PATH (or set QUARTO_PATH below)
#   - ASReml-R installed with a valid licence
# =============================================================================

# --- Configuration -----------------------------------------------------------

pkg_root   <- here::here()   # or set explicitly: pkg_root <- "path/to/wgaim"
vig_dir    <- file.path(pkg_root, "vignettes")
doc_dir    <- file.path(pkg_root, "inst", "doc")
quarto_bin <- Sys.getenv("QUARTO_PATH",
                         "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe")

# Vignettes to build: named vector of <source qmd> -> <output stem>
vignettes <- c(
    "qtlAim_pipeline.qmd",
    "qtlAim_mv_pipeline.qmd"
)

# --- Helpers -----------------------------------------------------------------

stop_if <- function(cond, msg) if (cond) stop(msg, call. = FALSE)

render_vignette <- function(qmd, quarto) {
    src <- file.path(vig_dir, qmd)
    stop_if(!file.exists(src),      paste("Source not found:", src))
    stop_if(!file.exists(quarto),   paste("Quarto binary not found:", quarto,
                                          "\nSet QUARTO_PATH env var or edit quarto_bin above."))
    message("\n--- Rendering: ", qmd, " ---")
    ret <- system2(quarto, args = c("render", src), stdout = "", stderr = "")
    if (ret != 0L) stop("Quarto render failed for: ", qmd, call. = FALSE)
    message("    OK")
}

copy_to_doc <- function(qmd) {
    stem <- sub("\\.qmd$", "", qmd)
    html_src  <- file.path(vig_dir, paste0(stem, ".html"))
    html_dst  <- file.path(doc_dir,  paste0(stem, ".html"))
    stop_if(!file.exists(html_src),
            paste("Rendered HTML not found:", html_src))
    if (!dir.exists(doc_dir)) dir.create(doc_dir, recursive = TRUE)
    file.copy(html_src, html_dst, overwrite = TRUE)
    message("    Copied -> ", html_dst)
}

# --- Main --------------------------------------------------------------------

message("=== wgAim vignette build ===")
message("Package root : ", pkg_root)
message("Quarto binary: ", quarto_bin)

for (qmd in vignettes) {
    render_vignette(qmd, quarto_bin)
    copy_to_doc(qmd)
}

message("\n=== Done. inst/doc/ is up to date. ===\n")
