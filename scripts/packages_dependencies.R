# ====================================================================== #
# packages_dependencies.R
#
# Installs the R packages used by the analyses and figures of the
# manuscript. Packages that are already installed are skipped.
#
# Usage, from the project root:
#   Rscript scripts/packages_dependencies.R
# ====================================================================== #

packages <- c("here",
              "ggplot2",
              "ggpubr",
              "patchwork",
              "svglite",
              "dplyr",
              "data.table",
              "png")

# Install only the packages that are missing. The CRAN mirror is set here
# because Rscript cannot open the menu that asks for one.
missing_packages <- packages[!packages %in% rownames(installed.packages())]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos="https://cloud.r-project.org")
}
