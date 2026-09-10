## 05_Integration_Subpopulation.R

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SUBTYPE <- if (length(args) >= 1) args[1] else "ER"
SUBPOP  <- if (length(args) >= 2) args[2] else "Myeloid"
stopifnot(SUBPOP %in% names(subpopulations))

cfg <- subpopulations[[SUBPOP]]
resolution <- res_subpop[[SUBPOP]][[SUBTYPE]]
nSample <- nrow(read_subtype_targets(SUBTYPE))
n_dims <- 30
tag <- paste0(SUBTYPE, "_", cfg$suffix)

## ---- Take the sub-population out of its compartment ----
so <- readRDS(file.path(RDS_DIR, paste0("so_", SUBTYPE, "_", cfg$parent, ".RDS")))
so <- so[, so@meta.data$cell_type_L2 %in% cfg$cell_types]
so[["RNA"]] <- split(so[["RNA"]], f = so$orig.ident)
message(tag, ": ", ncol(so), " nuclei")

## ---- Harmony integration ----
so <- integrate_cells(so, n_dims = n_dims, min_dist = 0.4)

save_figure(DimPlot(so, reduction = "umap", group.by = "orig.ident",
                    raster = FALSE, cols = col.p2, pt.size = 0.2, shuffle = TRUE),
            paste0(tag, "_UMAP_sample"), width = 6.5, height = 5.5)

## The markers of the lineages this sub-population is expected to contain.
markers_lineage <- read_marker_panel(cfg$marker_types)
markers_lineage <- markers_lineage[markers_lineage %in% rownames(so)]
save_figure(FeaturePlot(so, reduction = "umap", features = markers_lineage,
                        order = TRUE, pt.size = 0.3, ncol = 5, raster = FALSE),
            paste0(tag, "_UMAP_markers"),
            width = 20, height = ceiling(length(markers_lineage) / 5) * 3.5)

## ---- Cell clustering ----
so <- FindClusters(so, resolution = resolution, verbose = FALSE)
so@meta.data$clusters_L3 <- so@meta.data$seurat_clusters

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

## ---- GO and KEGG analysis of the cluster markers ----
save_go_kegg(go_kegg_by_cluster(marker_all), tag)

saveRDS(so, file.path(RDS_DIR, paste0("so_", tag, ".RDS")))
