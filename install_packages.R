packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "float",
  "stats",
  "mgcv",
  "MASS",
  "mvnfast",
  "tidyverse",
  "parallel"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg,repos = "https://cloud.r-project.org")
  }
}

invisible(sapply(packages, install_if_missing))
