# ============================================================================
# Shared project setup for the MAST32001 final assignment (Tasks 1 and 2).
# Source this file at the top of either analysis script, then add the
# task-specific packages with ensure_packages():
#
#   source("R/project_setup.R")
#   ensure_packages(c("JADE"))      # Task 1
#   ensure_packages(c("posterior")) # Task 2
#
# It locates the project root, creates the output directories, loads the common
# packages, installs the global plotting theme and provides save_sci() for
# publication-ready PNG exports.
# All figures are exported as PNG only; the Quarto reports include these PNGs.
# ============================================================================

## ---- 0. project root ------------------------------------------------------
# The root is the first ancestor directory that contains R/task1.R or
# R/task2.R, so the same block works for both scripts, whether they are run
# with Rscript, source(), or from RStudio.
locate_root <- function() {
  p <- grep("^--file=", commandArgs(), value = TRUE)
  p <- if (length(p))
    sub("^--file=", "", p[1L])
  else
    ""
  # source() also works outside RStudio and from another working directory.
  for (frame in sys.frames()) {
    if (!is.null(frame$ofile))
      p <- frame$ofile
  }
  if (!nzchar(p) &&
      interactive() &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active <- rstudioapi::getSourceEditorContext()$path
    if (length(active) && nzchar(active))
      p <- active
  }
  # Try the resolved script directory first; if it is not under the project
  # (e.g. a launcher script in /tmp), fall back to the working directory.
  starts <- c(if (nzchar(p) &&
                  file.exists(p))
    dirname(normalizePath(p, winslash = "/")), getwd())
  for (d in unique(starts)) {
    repeat {
      if (any(file.exists(file.path(d, c(
        "R/task1.R", "R/task2.R"
      )))))
        return(d)
      parent <- dirname(d)
      if (identical(parent, d))
        break
      d <- parent
    }
  }
  stop("Could not find the project root (no R/task1.R or R/task2.R).")
}

root <- locate_root()
setwd(root)

## ---- 1. project directories -----------------------------------------------
data_dir   <- "data"            # raw inputs
table_dir  <- "results/tables"  # tables and csv
figure_dir <- "results/figures" # figures as PNG
result_dir <- "results"         # RData workspace
report_dir <- "report"          # rendered report, quarto output
for (folder in c(data_dir, table_dir, figure_dir, result_dir, report_dir))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

## ---- 2. packages ----------------------------------------------------------
# Common packages used by both tasks; missing ones are installed first.
common_pkgs <- c("tidyverse", "ggpubr")

ensure_packages <- function(extra = character()) {
  pkgs <- unique(c(common_pkgs, extra))
  installed <- rownames(installed.packages())
  to_install <- setdiff(pkgs, installed)
  if (length(to_install))
    install.packages(to_install)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
# Load the common packages right away: theme_set() and ggsave() below need
# ggplot2. Task-specific packages are added by ensure_packages(extra) in the
# analysis scripts.
ensure_packages()

## ---- 3. global plotting style ---------------------------------------------
# Okabe-Ito palette: blue, orange, green and purple first for contrast.
plot_palette <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7",
  "#E69F00",
  "#56B4E9",
  "#F0E442",
  "#000000"
)
# Backwards-compatible alias used by the Task 1 report code blocks.
PALETTE <- plot_palette

theme_set(
  theme_pubr(base_size = 11, base_family = "sans") +
    theme(
      palette.colour.discrete = plot_palette,
      palette.fill.discrete = plot_palette,
      axis.text = element_text(size = 10, colour = "#202020"),
      axis.title = element_text(size = 11),
      axis.line = element_line(linewidth = .3),
      axis.ticks = element_line(linewidth = .3),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.key.width = grid::unit(8, "mm"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(
        fill = "grey96",
        colour = "grey75",
        linewidth = .3
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
)

## ---- 4. PNG export helper ---------------------------------------------------
# Export at double (183 mm) or single (89 mm) column width, 600 dpi.
save_sci <- function(plot,
                     filename,
                     width = c("single", "double"),
                     height_mm = 68,
                     directory = figure_dir) {
  width <- match.arg(width)
  width_mm <- if (width == "single")
    89
  else
    183
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  
  out <- file.path(directory, paste0(filename, ".png"))
  ggsave(
    out,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    dpi = 600,
    device = grDevices::png
  )
  
  invisible(out)
}