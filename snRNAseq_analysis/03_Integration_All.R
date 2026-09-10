## 03_Integration_All.R

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SUBTYPE <- if (length(args) >= 1) args[1] else "ER"

targets <- read_subtype_targets(SUBTYPE)
nSample <- nrow(targets)
resolution <- res_all[[SUBTYPE]]
n_dims <- 30

## ---- Pool the samples ----
meta_col <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
              "db_score", "db_type", "S.Score", "G2M.Score", "Phase",
              "seurat_clusters", "infercnv_score")
counts_all <- meta_all <- list()

for (i in seq_len(nSample)) {
  sample_id <- targets$Sample[i]
  so_i <- readRDS(file.path(RDS_DIR, paste0("so_", sample_id, ".RDS")))
  cnv <- readRDS(file.path(RDS_DIR, paste0(sample_id, "_instability_score.rds")))
  so_i@meta.data$infercnv_score <- cnv[match(colnames(so_i), names(cnv))]
  counts_all[[i]] <- so_i@assays$RNA$counts
  meta_all[[i]] <- so_i@meta.data[, meta_col]
}
names(counts_all) <- names(meta_all) <- targets$Sample

so <- CreateSeuratObject(counts = do.call(cbind, counts_all),
                         meta.data = do.call(rbind, meta_all),
                         assay = "RNA", min.cells = 0, min.features = 0)
so@meta.data$individual_clusters <- so@meta.data$seurat_clusters
so <- NormalizeData(so, verbose = FALSE)
so[["RNA"]] <- split(so[["RNA"]], f = so$orig.ident)
rm(counts_all)

## ---- Harmony integration ----
so <- integrate_cells(so, n_dims = n_dims, min_dist = 0.4)

save_figure(DimPlot(so, reduction = "umap", group.by = "orig.ident",
                    raster = FALSE, cols = col.p2, pt.size = 0.2, shuffle = TRUE),
            paste0(SUBTYPE, "_all_UMAP_sample"), width = 6.5, height = 5.5)

## ---- Cell clustering ----
so <- FindClusters(so, resolution = resolution, verbose = FALSE)
so@meta.data$clusters_L1 <- so@meta.data$seurat_clusters

save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 0.2),
            paste0(SUBTYPE, "_all_UMAP_cluster"), width = 6.5, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 0.2, split.by = "orig.ident", ncol = 3),
            paste0(SUBTYPE, "_all_UMAP_cluster_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)

save_figure(infercnv_boxplot(so, "seurat_clusters", title = SUBTYPE),
            paste0(SUBTYPE, "_all_inferCNV_byCluster"), width = 7, height = 4)

## ---- Cluster markers ----
so <- JoinLayers(so)
marker_all <- FindAllMarkers(so, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
marker_all <- marker_all[marker_all$p_val_adj < 0.05, ]
saveRDS(marker_all, file.path(GENELIST_DIR, paste0(SUBTYPE, "_all_marker.RDS")))
write.csv(marker_all, file.path(GENELIST_DIR, paste0(SUBTYPE, "_all_marker.csv")))

save_figure(DotPlot(so, features = top_markers(marker_all, n = 5), dot.scale = 8) +
              coord_flip(),
            paste0(SUBTYPE, "_all_DotPlot_topMarkers"), width = 8, height = 13)

## ---- Level 1 cell type annotation ----
ct_All <- read_celltype_annotation(SUBTYPE, "All")
cell_type <- ct_All$cell_type[match(as.character(so@meta.data$seurat_clusters),
                                    as.character(ct_All$cluster))]
ct_names <- names(cols_Somi)[names(cols_Somi) %in% cell_type]
so@meta.data$cell_type_L1 <- factor(cell_type, levels = ct_names)
print(t(t(table(so@meta.data$cell_type_L1))))

save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_L1",
                    raster = FALSE, cols = cols_Somi[ct_names], pt.size = 0.2),
            paste0(SUBTYPE, "_all_UMAP_celltype"), width = 7, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_L1",
                    raster = FALSE, cols = cols_Somi[ct_names], pt.size = 0.2,
                    split.by = "orig.ident", ncol = 3),
            paste0(SUBTYPE, "_all_UMAP_celltype_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)
save_figure(composition_barplot(so, "cell_type_L1", cols_Somi[ct_names]),
            paste0(SUBTYPE, "_all_Barplot_celltype"), width = 8, height = 5)

saveRDS(so, file.path(RDS_DIR, paste0("so_", SUBTYPE, "_all.RDS")))
