## 00_Install_Packages.R

# CRAN
install.packages(c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "ggrepel",
  "uwot",
  "harmony",
  "gridExtra",
  "scales",
  "pheatmap",
  "rjags",
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
  "SummarizedExperiment",
  "ComplexHeatmap",
  "scDblFinder",
  "infercnv",
  "org.Hs.eg.db",
  "AnnotationDbi"
))

