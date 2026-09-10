## 00_Install_Packages.R

# CRAN
install.packages(c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "uwot",
  "harmony",
  "gridExtra",
  "scales",
  "janitor",
  "pheatmap",
  "circlize",
  "viridis",
  "hexDensity",
  "Matrix",
  "remotes"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "edgeR",
  "limma",
  "SingleCellExperiment",
  "SpatialExperiment",
  "SummarizedExperiment",
  "ComplexHeatmap"
))


remotes::install_github("ChenLaboratory/scider")
remotes::install_github("jinworks/CellChat")

