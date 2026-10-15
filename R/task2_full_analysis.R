# ============================================================================
#
#
#
#
#
# ============================================================================

## ---- 0. output layout -----------------------------------------------------
data_dir <- "data"                  # raw inputs
table_dir <- "results/tables"       # tables ad csv
figure_dir <- "results/figures"     # figures as PDF + PNG
result_dir <- "results"             # RData workspace
report_dir <- "report"              # rendered report, quarto output
dir.create(table_dir, recursive = T, showWarnings = F)
dir.create(figure_dir, recursive = T, showWarnings = F)
dir.create(result_dir, recursive = T, showWarnings = F)
dir.create(report_dir, recursive = T, showWarnings = F)

## ---- 1. setup:libraries and visualization style ---------------------------
# list all packages needed
pkgs <- c("tidyverse", "scales", "ggpubr", "coda", "ggsci")

# get the names of all currently installed packages
installed <- rownames(installed.packages())
# find packages that are not installed yet
to_install <- setdiff(pkgs, installed)

# install only the missing packages
if (length(to_install) > 0) {
  install.packages(to_install)
}

# load the packages
# characer.only = TRUE lets library() accept a package name as a character string
invisible(lapply(pkgs, library, character.only = TRUE))

# --- publication ready style: global settings (fonts, sizes, colors) --------
