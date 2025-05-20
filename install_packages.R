packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "float",
  "stats",
  "mgcv",
  "MASS",
  "mvnfast"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(sapply(packages, install_if_missing))
