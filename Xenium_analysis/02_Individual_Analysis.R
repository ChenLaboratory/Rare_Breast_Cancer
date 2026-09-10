## 02_Individual_Analysis.R
##
## Per sample Xenium analysis, in two parts:
##   (a) cell level clustering and manual cell type annotation;
##   (b) hexagonal bin (niche) analysis of the annotated tissue.
##
## Run once per sample. The sample and the stage can be given on the command
## line, or edited below when running interactively:
##   Rscript 02_Individual_Analysis.R APO_3
##   Rscript 02_Individual_Analysis.R APO_3 post_CAF
##
## STAGE = "main"
##   Cell type labels come from metadata/<Subtype>_celltype_Individual.txt.
## STAGE = "post_CAF"
##   The cell objects of the main pass are reused, the CAF clusters found by
##   05_Integration_Cell_CAF.R replace the pooled fibroblast and CAF labels,
##   and the niche analysis is repeated on the refined annotation.

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SAMPLE_ID <- if (length(args) >= 1) args[1] else "APO_3"
STAGE     <- if (length(args) >= 2) args[2] else "main"
stopifnot(STAGE %in% c("main", "post_CAF"))

tg <- targets_for_sample(SAMPLE_ID)
out_dir <- stage_dir(STAGE)
fig <- function(name) paste0(SAMPLE_ID, "_", name, if (STAGE == "post_CAF") "_post_CAF" else "")

n_dims <- 30


## ========================================================================
## Cell level analysis. The main stage derives the cell types from scratch,
## the post CAF stage refines the cell types of the main stage.
## ========================================================================

if (STAGE == "main") {

## ---- Read the ProSeg segmentation ----
 seg_out <- file.path(SEG_DIR, tg$Proseg, tg$Proseg, "outs")
  raw <- Read10X_h5(file.path(seg_out, "cell_feature_matrix.h5"))
  counts <- if (is.list(raw)) raw[[1]] else raw
  cell_info <- read.csv(file.path(seg_out, "cells.csv.gz"))
  rownames(cell_info) <- cell_info$cell_id
  cell_info$cell_id <- NULL
  so <- CreateSeuratObject(counts = counts, project = SAMPLE_ID, assay = "RNA",
                           meta.data = cell_info, min.cells = 0, min.features = 0)

## ---- Cell filtering ----
keep_cell <- so@meta.data$nCount_RNA >= tg$vcut & so@meta.data$nFeature_RNA >= tg$hcut
message(SAMPLE_ID, ": ", sum(keep_cell), " / ", ncol(so), " cells pass filtering")
so <- so[, keep_cell]

## ---- Normalisation, PCA and UMAP ----
so <- NormalizeData(so, scale.factor = 100, verbose = FALSE)
so <- ScaleData(so, features = rownames(so), verbose = FALSE)
so <- RunPCA(so, features = rownames(so), npcs = n_dims, verbose = FALSE)
so <- run_umap_from_knn(so, dims = 1:n_dims, min_dist = 0.7)

## ---- Cell clustering ----
so <- FindNeighbors(so, dims = 1:n_dims, verbose = FALSE)
so <- FindClusters(so, resolution = tg$resCell, verbose = FALSE)

## ---- Cluster markers ----
marker_all <- FindAllMarkers(so, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
marker_all <- marker_all[marker_all$p_val_adj < 0.05, ]
saveRDS(marker_all, file.path(GENELIST_DIR, paste0(SAMPLE_ID, "_cluster_marker.RDS")))
write.csv(marker_all, file.path(GENELIST_DIR, paste0(SAMPLE_ID, "_cluster_marker.csv")))

save_figure(DotPlot(so, features = top_markers(marker_all, n = 5), dot.scale = 8) +
              coord_flip(),
            fig("DotPlot_cluster_topMarkers"), width = 8, height = 13)

## ---- Manual cell type annotation ----
ct_map <- read_celltype_annotation(tg$Subtype)
ct_map <- ct_map[ct_map$sample == SAMPLE_ID, ]
cell_type <- ct_map$cell_type[match(as.character(so@meta.data$seurat_clusters),
                                    as.character(ct_map$cluster))]
ct_names <- names(cols_Somi)[names(cols_Somi) %in% cell_type]
so@meta.data$cell_type_manual <- factor(cell_type, levels = ct_names)
}


if (STAGE == "post_CAF") {

## ---- Reuse the cell objects of the main pass ----
# so <- readRDS(file.path(RDS_DIR, paste0("so_proseg_", SAMPLE_ID, ".RDS")))

## Replace the pooled fibroblast and CAF labels by the CAF clusters of
## 05_Integration_Cell_CAF.R, translated into CAF subtypes through
## metadata/<Subtype>_CAF_celltype.txt.
so_caf <- readRDS(file.path(RDS_DIR, paste0("so_cell_caf_", tg$Subtype, ".RDS")))
caf_map <- read_caf_annotation(tg$Subtype)

in_sample <- so_caf@meta.data$orig.ident == SAMPLE_ID
caf_cells <- sub(paste0("^", SAMPLE_ID, "-"), "", colnames(so_caf)[in_sample])
caf_label <- caf_map$cell_type[match(as.character(so_caf@meta.data$seurat_clusters[in_sample]),
                                     as.character(caf_map$cluster))]

## CAF cells dropped by the small cluster filter of script 05
## keep the cell type they were given in the main pass.
cell_type <- as.character(so@meta.data$cell_type_manual)
m <- match(caf_cells, colnames(so))
relabel <- !is.na(m) & !is.na(caf_label)
cell_type[m[relabel]] <- caf_label[relabel]

ct_names <- names(cols_Somi)[names(cols_Somi) %in% cell_type]
so@meta.data$cell_type_manual <- factor(cell_type, levels = ct_names)
print(t(t(table(so@meta.data$cell_type_manual))))

}


## ---- Cell type annotation on the UMAP and on the tissue ----
so@reductions$spatial <- so@reductions$umap
Key(so@reductions$spatial) <- "Sp_"
so@reductions$spatial@cell.embeddings[, 1] <- so@meta.data$x_centroid
so@reductions$spatial@cell.embeddings[, 2] <- so@meta.data$y_centroid

save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_manual",
                    raster = FALSE, cols = cols_Somi[ct_names], pt.size = 0.2) +
              ggtitle("Cell type annotation"),
            fig("UMAP_celltype"), width = 11, height = 7.5)

save_figure(DimPlot(so, reduction = "spatial", group.by = "cell_type_manual",
                    raster = FALSE, cols = cols_Somi[ct_names], pt.size = 0.03) +
              coord_fixed(),
            fig("Spatial_celltype"), width = tg$width, height = tg$height)

## ---- SpatialExperiment object ----
sce <- as.SingleCellExperiment(so)
colData(sce)$cell_type <- so@meta.data$cell_type_manual
colData(sce)$umap_1 <- Embeddings(so, "umap")[, 1]
colData(sce)$umap_2 <- Embeddings(so, "umap")[, 2]
spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = assay(sce, 1)),
  sample_id = SAMPLE_ID, colData = colData(sce),
  spatialCoordsNames = c("x_centroid", "y_centroid"))
spe$cell_type <- factor(spe$cell_type, levels = ct_names)

saveRDS(so,  file.path(out_dir, paste0("so_proseg_", SAMPLE_ID, ".RDS")))
saveRDS(spe, file.path(out_dir, paste0("spe_proseg_", SAMPLE_ID, ".RDS")))


## ---- Cell type density across the tissue ----
spe <- gridDensity(spe, grid.length.x = 50, bandwidth = 50)

CT <- list(present_or_default(Tumor_CT, ct_names, "Tumor"),
           present_or_default(Fibroblast_CT, ct_names, "Fibroblast"),
           "Myeloid",
           present_or_default(T_CT, ct_names, "T"))

## ---- Cellular niche profiles ----
## Cells are aggregated into hexagonal bins
## Bins with a small library size are dropped.
spe_gr <- gridSPE(spe, cell.count = TRUE)
keep_gr <- spe_gr@colData$LibSize >= tg$cutGrid
message(SAMPLE_ID, ": ", sum(keep_gr), " / ", ncol(spe_gr), " bins pass filtering")
spe_gr <- spe_gr[, keep_gr]

## ---- Niche clustering ----
so_gr <- CreateSeuratObject(spe_gr@assays@data$counts, min.cells = 0, min.features = 0)
so_gr <- NormalizeData(so_gr, scale.factor = 100, verbose = FALSE)
so_gr <- ScaleData(so_gr, features = rownames(so_gr), verbose = FALSE)
so_gr <- RunPCA(so_gr, features = rownames(so_gr), npcs = n_dims, verbose = FALSE)
so_gr <- run_umap_from_knn(so_gr, dims = 1:n_dims, min_dist = 0.7)
so_gr <- FindClusters(so_gr, resolution = tg$resGrid, verbose = FALSE)
print(table(Idents(so_gr)))

colData(spe_gr)$niche <- so_gr@meta.data$seurat_clusters
save_figure(plotGrid(spe_gr, group.by = "niche", cols = col.pMedium, pol.border = FALSE),
            fig("Spatial_niche"), width = tg$width * 1.4, height = tg$height * 1.4)

## ---- Cell type composition of the niches ----
ct_names_noMixed <- ct_names[tolower(ct_names) != "mixed"]
ct_clean <- janitor::make_clean_names(ct_names_noMixed)

cct <- rowsum(as.matrix(colData(spe_gr)$cell_count), colData(spe_gr)$niche)
print(t(cct))
perc <- t(cct[, ct_clean, drop = FALSE] / cct[, "overall"]) * 100

dat <- data.frame(
  Niche      = factor(rep(colnames(perc), each = length(ct_names_noMixed)),
                      levels = colnames(perc)),
  Celltype   = factor(rep(ct_names_noMixed, ncol(perc)), levels = ct_names_noMixed),
  Proportion = as.vector(perc))
p <- ggplot(dat, aes(fill = Celltype, y = Proportion, x = Niche)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cols_Somi[ct_names_noMixed]) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity",
           color = "black", linewidth = 0.2) +
  guides(fill = guide_legend(ncol = 2)) +
  coord_flip() + theme_classic()
save_figure(p, fig("Barplot_niche_cellProp"), width = 9, height = 4)

## ---- Cell type density across the niches ----
den_cols <- paste0("density_", ct_clean)
den_gr <- as.matrix(spe@metadata$grid_density[keep_gr, den_cols])
colnames(den_gr) <- ct_names_noMixed
den_gr <- rowsum(den_gr, colData(spe_gr)$niche) / as.vector(table(colData(spe_gr)$niche))
den_gr <- t(den_gr)

annot_col <- data.frame(Niche = levels(colData(spe_gr)$niche))
rownames(annot_col) <- colnames(den_gr)
ann_colors <- list(Niche = col.pMedium[seq_len(nlevels(colData(spe_gr)$niche))])
names(ann_colors$Niche) <- levels(colData(spe_gr)$niche)

save_figure(pheatmap::pheatmap(t(scale(t(den_gr))),
              color = colorRampPalette(col.spec[10:1])(100), border_color = NA,
              breaks = seq(-2, 2, length.out = 101),
              cluster_cols = TRUE, cluster_rows = TRUE, scale = "none",
              fontsize_row = 10, fontsize_col = 10, show_colnames = TRUE,
              treeheight_row = 70, treeheight_col = 70, cutree_cols = 1,
              annotation_col = annot_col, annotation_colors = ann_colors,
              main = "Niche", clustering_method = "ward.D2", silent = TRUE),
            fig("Heatmap_niche_density"), width = 7, height = 8)

## ---- Niche markers ----
niche_markers <- FindAllMarkers(so_gr, only.pos = TRUE)
niche_markers <- niche_markers[niche_markers$p_val_adj < 0.05, ]
saveRDS(niche_markers, file.path(GENELIST_DIR, paste0(SAMPLE_ID, "_niche_marker.RDS")))
write.csv(niche_markers, file.path(GENELIST_DIR, paste0(SAMPLE_ID, "_niche_marker.csv")))

save_figure(DotPlot(so_gr, features = top_markers(niche_markers, n = 5), dot.scale = 8) +
              coord_flip(),
            fig("DotPlot_niche_topMarkers"), width = 8, height = 13)

saveRDS(spe_gr, file.path(out_dir, paste0("spe_gr_", SAMPLE_ID, ".RDS")))
saveRDS(so_gr,  file.path(out_dir, paste0("so_gr_", SAMPLE_ID, ".RDS")))
