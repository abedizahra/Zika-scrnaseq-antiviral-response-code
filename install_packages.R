# Install R package dependencies for this repository
# Run from the repository root before executing scripts/run_all.R

cran_packages <- c('Cairo', 'Matrix', 'Seurat', 'data.table', 'dplyr', 'future', 'ggplot2', 'magick', 'patchwork', 'pheatmap', 'purrr', 'ragg', 'readr', 'scales', 'stringr', 'tibble', 'tidyr', 'tidyverse', 'writexl')

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

bioc_packages <- c('DESeq2')
if (length(bioc_packages) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

message("Dependency check complete.")
