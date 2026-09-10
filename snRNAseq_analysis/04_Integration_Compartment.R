## 04_Integration_Compartment.R

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SUBTYPE     <- if (length(args) >= 1) args[1] else "ER"
COMPARTMENT <- if (length(args) >= 2) args[2] else "Immune"
stopifnot(COMPARTMENT %in% names(compartments))

cfg <- compartments[[COMPARTMENT]]
resolution <- res_compartment[[COMPARTMENT]][[SUBTYPE]]
nSample <- nrow(read_subtype_targets(SUBTYPE))
n_dims <- 30
tag <- paste0(SUBTYPE, "_", cfg$suffix)

## ---- Take the compartment out of the all cell integration ----
so <- readRDS(file.path(RDS_DIR, paste0("so_", SUBTYPE, "_all.RDS")))
so <- so[, so@meta.data$cell_type_L1 %in% cfg$cell_types]
so[["RNA"]] <- split(so[["RNA"]], f = so$orig.ident)
message(tag, ": ", ncol(so), " nuclei")

## ---- Harmony integration ----
so <- integrate_cells(so, n_dims = n_dims, min_dist = 0.4)

save_figure(DimPlot(so, reduction = "umap", group.by = "orig.ident",
                    raster = FALSE, cols = col.p2, pt.size = 0.2, shuffle = TRUE),
            paste0(tag, "_UMAP_sample"), width = 6.5, height = 5.5)

if (COMPARTMENT != "Immune") {
  save_figure(FeaturePlot(so, reduction = "umap", features = EMT_Genes,
                          order = TRUE, pt.size = 0.3, ncol = 5, raster = FALSE),
              paste0(tag, "_UMAP_EMT"), width = 20, height = 10.5)
}

## ---- Cell clustering ----
so <- FindClusters(so, resolution = resolution, verbose = FALSE)
so@meta.data$clusters_L2 <- so@meta.data$seurat_clusters

save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 0.2),
            paste0(tag, "_UMAP_cluster"), width = 6.5, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 0.2, split.by = "orig.ident", ncol = 3),
            paste0(tag, "_UMAP_cluster_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)
save_figure(composition_barplot(so, "seurat_clusters", col.p2, "Cluster"),
            paste0(tag, "_Barplot_cluster"), width = 8, height = 5)

save_figure(infercnv_boxplot(so, "seurat_clusters", title = SUBTYPE),
            paste0(tag, "_inferCNV_byCluster"), width = 7, height = 4)

## ---- Cluster markers ----
so <- JoinLayers(so)
marker_all <- FindAllMarkers(so, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
marker_all <- marker_all[marker_all$p_val_adj < 0.05, ]
saveRDS(marker_all, file.path(GENELIST_DIR, paste0(tag, "_marker.RDS")))
write.csv(marker_all, file.path(GENELIST_DIR, paste0(tag, "_marker.csv")))

save_figure(DotPlot(so, features = top_markers(marker_all, n = 5), dot.scale = 8) +
              coord_flip(),
            paste0(tag, "_DotPlot_topMarkers"), width = 8, height = 13)

## ---- Level 2 cell type annotation ----
ct <- read_celltype_annotation(SUBTYPE, COMPARTMENT)
cell_type <- ct$cell_type[match(as.character(so@meta.data$seurat_clusters),
                                as.character(ct$cluster))]
ct_names <- names(cols_Somi)[names(cols_Somi) %in% cell_type]
so@meta.data$cell_type_L2 <- factor(cell_type, levels = ct_names)
print(t(t(table(so@meta.data$cell_type_L2))))

save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_L2",
                    raster = FALSE, cols = cols_Somi[ct_names], pt.size = 0.2),
            paste0(tag, "_UMAP_celltype"), width = 7, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_L2",
                    raster = FALSE, cols = cols_Somi[ct_names], pt.size = 0.2,
                    split.by = "orig.ident", ncol = 3),
            paste0(tag, "_UMAP_celltype_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)
save_figure(composition_barplot(so, "cell_type_L2", cols_Somi[ct_names]),
            paste0(tag, "_Barplot_celltype"), width = 8, height = 5)

saveRDS(so, file.path(RDS_DIR, paste0("so_", tag, ".RDS")))
