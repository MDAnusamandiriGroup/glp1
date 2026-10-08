required_packages <- c(
  "deSolve", "haven", "dplyr", "survey", "ggplot2", "tidyr", "tidyselect",
  "testthat"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
