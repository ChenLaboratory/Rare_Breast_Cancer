## 05_Integration_Cell_CAF.R
##
## The CAF clusters found here are turned into CAF subtype labels in
## metadata/<Subtype>_CAF_celltype.txt and fed back into the post CAF pass of
## scripts 02 and 04.
##
## Run once per subtype:
##   Rscript 05_Integration_Cell_CAF.R CYS

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SUBTYPE <- if (length(args) >= 1) args[1] else "CYS"

targets <- read_subtype_targets(SUBTYPE)
nSample <- nrow(targets)
print(targets[, 1:3])

resolution <- caf_res[[SUBTYPE]]
n_dims <- 30
min_cluster_size <- 20

## ---- Collect the fibroblasts and CAFs of every sample ----
meta_col <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "seurat_clusters",
              "cell_type_manual")
counts_all <- meta_all <- spe_list <- list()

for (i in seq_len(nSample)) {
  sample_id <- targets$SampleID[i]
  so_i  <- readRDS(file.path(RDS_DIR, paste0("so_proseg_", sample_id, ".RDS")))
  spe_i <- readRDS(file.path(RDS_DIR, paste0("spe_proseg_", sample_id, ".RDS")))
  colnames(spe_i) <- paste(sample_id, colnames(spe_i), sep = "-")

  so_i@meta.data$orig.ident <- factor(sample_id)
  is_CAF <- so_i@meta.data$cell_type_manual %in% Stroma_Fibroblast_CT

  counts_all[[i]] <- so_i@assays$RNA$counts[, is_CAF]
  meta_all[[i]] <- so_i@meta.data[is_CAF, meta_col]
  rownames(meta_all[[i]]) <- colnames(counts_all[[i]]) <-
    paste(sample_id, rownames(meta_all[[i]]), sep = "-")
  spe_list[[i]] <- spe_i
}
names(counts_all) <- names(meta_all) <- names(spe_list) <- targets$SampleID
print(vapply(counts_all, ncol, integer(1)))

## ---- Pool the cells ----
so <- CreateSeuratObject(counts = do.call(cbind, counts_all),
                         meta.data = do.call(rbind, meta_all),
                         assay = "RNA", min.cells = 0, min.features = 0)
so@meta.data$individual_clusters <- so@meta.data$seurat_clusters
so <- NormalizeData(so, scale.factor = 100, verbose = FALSE)
so[["RNA"]] <- split(so[["RNA"]], f = so$orig.ident)
rm(counts_all)

## ---- Harmony integration ----
VariableFeatures(so) <- rownames(so)
so <- ScaleData(so, features = rownames(so), verbose = FALSE)
so <- RunPCA(so, features = rownames(so), npcs = n_dims, verbose = FALSE)
so <- IntegrateLayers(so, method = HarmonyIntegration,
                      orig.reduction = "pca", new.reduction = "integrated.harmony",
                      dims = 1:n_dims, k.anchor = 20, verbose = FALSE)

so <- run_umap_from_knn(so, reduction = "integrated.harmony",
                        dims = 1:n_dims, min_dist = 0.4)

save_figure(DimPlot(so, reduction = "umap", group.by = "orig.ident",
                    raster = FALSE, cols = col.p2, pt.size = 0.2, shuffle = TRUE),
            paste0(SUBTYPE, "_CAF_UMAP_sample"), width = 6.5, height = 5.5)

save_figure(FeaturePlot(so, reduction = "umap", features = CAF_Genes,
                        raster = FALSE, pt.size = 0.3, ncol = 5, order = TRUE),
            paste0(SUBTYPE, "_CAF_UMAP_CAFgenes"), width = 20, height = 10.5)

## ---- CAF clustering ----
so <- FindClusters(so, resolution = resolution, verbose = FALSE)
print(table(so@meta.data$seurat_clusters))

clst <- so@meta.data$seurat_clusters
keep_clst <- levels(clst)[table(clst) >= min_cluster_size]
so <- so[, clst %in% keep_clst]
so@meta.data$seurat_clusters <- Idents(so) <-
  factor(so@meta.data$seurat_clusters, levels = keep_clst)

save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 1),
            paste0(SUBTYPE, "_CAF_UMAP_cluster"), width = 9, height = 7.5)
save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 0.8, split.by = "orig.ident", ncol = 3),
            paste0(SUBTYPE, "_CAF_UMAP_cluster_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)

save_figure(VlnPlot(so, features = CAF_Genes, cols = col.p2, pt.size = 0, ncol = 3),
            paste0(SUBTYPE, "_CAF_Violin_CAFgenes"), width = 18, height = 20)

## ---- CAF clusters on the tissue ----
ncls <- nlevels(so@meta.data$seurat_clusters)
plist <- list()
for (i in seq_len(nSample)) {
  spe_i <- spe_list[[i]]
  metadata <- so@meta.data[so@meta.data$orig.ident == targets$SampleID[i], ]
  m <- match(rownames(metadata), colnames(spe_i))
  spe_i$CAF_cluster <- rep("Others", ncol(spe_i))
  spe_i$CAF_cluster[m] <- as.character(metadata$seurat_clusters)
  spe_i$CAF_cluster <- factor(spe_i$CAF_cluster,
                              levels = c("Others", seq_len(ncls) - 1))
  plist[[i]] <- plotSpatial(spe_i, group.by = "CAF_cluster",
                            cols = c("grey70", col.p2), pt.alpha = 0.8) +
                ggtitle(targets$SampleID[i])
}
save_figure(plist, paste0(SUBTYPE, "_CAF_Spatial_cluster"), width = 18,
            height = ceiling(nSample / 3) * 7, ncol = 3)

## ---- CAF cluster markers ----
so <- JoinLayers(so)
marker_all <- FindAllMarkers(so, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
marker_all <- marker_all[marker_all$p_val_adj < 0.05, ]
saveRDS(marker_all, file.path(GENELIST_DIR, paste0(SUBTYPE, "_CAF_marker.RDS")))
write.csv(marker_all, file.path(GENELIST_DIR, paste0(SUBTYPE, "_CAF_marker.csv")))

save_figure(DotPlot(so, features = top_markers(marker_all, n = 5), dot.scale = 8) +
              coord_flip(),
            paste0(SUBTYPE, "_CAF_DotPlot_topMarkers"), width = 8, height = 13)

## ---- Pseudo-bulk expression of the CAF clusters ----
y1 <- seurat2pb(so, sample = "orig.ident")
colnames(y1) <- gsub("cluster", "C", colnames(y1))
y2 <- sumTechReps(y1, ID = factor(paste0("C", y1$samples$cluster)))

write.csv(edgeR::cpm(y1, log = TRUE),
          file.path(GENELIST_DIR, paste0(SUBTYPE, "_CAF_logCPM_sample_cluster.csv")))
write.csv(edgeR::cpm(y2, log = TRUE),
          file.path(GENELIST_DIR, paste0(SUBTYPE, "_CAF_logCPM_cluster.csv")))

saveRDS(so, file.path(RDS_DIR, paste0("so_cell_caf_", SUBTYPE, ".RDS")))
