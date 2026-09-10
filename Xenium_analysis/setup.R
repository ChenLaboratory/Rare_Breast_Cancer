## setup.R
##

library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scider)
library(edgeR)
library(ggplot2)
library(uwot)
library(harmony)
library(gridExtra)
library(scales)
library(janitor)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)

options(width = 90, digits = 3)
options(future.globals.maxSize = 1.2e+10)

## ---- Colour palettes ----
col.pMedium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
    "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA")
col.pDark <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
col.pLight <- c("#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
    "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5")
col.p <- col.pMedium
col.p2 <- rep(c(col.pDark, col.pLight, col.pMedium), 2)

## Diverging palette used by all expression and enrichment heatmaps.
col.spec <- c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF",
    "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")

## Palette for the common niche clusters of the grid level integration.
col.nclst <- rep(c("#001F4D", "#DAA520", "red3", "#2F7C6D", "#56B4E9",
    "#316684", "#999999", "orchid1", "#E0E0DF"), 2)
names(col.nclst) <- seq_along(col.nclst) - 1

## ---- Cell type colour map ----
cols_Somi <- c(
    "#1F77B4", "#AEC7E8", "#67BF5C", "#FF7F0E", "#E377C2", "#8C564B", "#C49C94", "#C4AE94", "#8C564B",
    "#FFFF66", "#FFFFCC", "#A2A2A2", "#F7B6D1", "#9467BD", "#17BECF", "#9EDAE5", "#BCBD22", "#A2A2A2",
    "#DBDB8D", "#D62728", "#FF9896", "#33004D", "#731ACC", "#C5B0D5",
    "#ED97CA", "#CDCC5D", "#AD8BC9", "#C5B0D5", "#E6CCFF", "#7F7F7F",
    "#FF6666", "#FF6666", "#660000", "#CC0000", "#FFCCCC", "#FF9999",
    "#CCCC00", "#999900", "#999900", "#74BAAD",
    "#1F77B4", "#2F7C6D", "#56B4E9", "#49B39E", "#316684", "#99CCFF",
    "#AEC7E8", "#AEE8CC",
    "#1F77B4", "#2F7C6D", "#56B4E9", "#999999", "#99CCFF", "#AEE8CC",
    "#FF7F0E",
    "#D55E00", "#FF7F0E", "#FFB266", "#FFCC99", "#FFE5CC", "#E67300",
    "#001F4D", "#99E699", "#009900", "#67BF5C", "#C7C7C7", "#C7C7C7")
names(cols_Somi) <- c(
    "Tumor", "Cycling_Tumor", "Epithelial", "Fibroblast", "Adipocyte", "Endothelial", "Pericyte", "PVL", "VEC",
    "FDC", "Cycling_FDC", "Lymph.vessel", "Muscle", "Macrophage", "Plasma", "Cycling_Plasma", "B", "LEC",
    "Cycling_B", "T", "Cycling_T", "Dendritic", "pDC", "Cycling_Myeloid",
    "Mast", "Keratinocyte", "Myeloid", "TAM", "Neutrophil", "Erythroid",
    "Na\u00efve_T", "CD4_T", "CD8_T", "NK", "Treg", "Ex_T",
    "Na\u00efve_B", "GC_B", "Memory_B", "Monocyte",
    "Tumor_1", "Tumor_2", "Tumor_3", "Tumor_4", "Tumor_5", "Tumor_6",
    "Cycling_Tumor_1", "Cycling_Tumor_2",
    "Tumor_Mes", "Tumor_Epi", "Tumor_ML", "Tumor_Adi", "Tumor_Mix", "Tumor_Sq",
    "CAF",
    "iCAF", "myCAF", "apCAF", "intCAF", "cCAF", "FRC",
    "Basal", "LP", "ML", "Luminal", "Mixed", "Contamination")

## ---- Cell type groups ----
Stroma_Fibroblast_CT <- c("Fibroblast", "CAF", "iCAF", "myCAF", "apCAF",
    "intCAF", "cCAF", "FRC")

Tumor_CT <- c("Tumor", "Cycling_Tumor", "Tumor_Mes", "Tumor_Sq", "Tumor_ML",
    "Tumor_Adi", "Tumor_1", "Tumor_2", "Tumor_3", "Tumor_4", "Tumor_5", "Tumor_6")
Fibroblast_CT <- c("Fibroblast", "iCAF", "myCAF", "CAF")
T_CT <- c("Na\u00efve_T", "CD4_T", "CD8_T", "NK", "Treg", "T", "Naive_T")

## ---- CAF marker gene sets ----
iCAF_Genes <- c("IGF1", "CD34", "CXCL12", "PDGFRB", "PTGDS", "FBLN1",
    "CCDC80", "ADH1B", "EGR1")
myCAF_Genes <- c("LRRC15", "LOX", "POSTN", "GJB2", "ACTA2", "MYLK")
CAF_Genes <- c(iCAF_Genes, myCAF_Genes)
