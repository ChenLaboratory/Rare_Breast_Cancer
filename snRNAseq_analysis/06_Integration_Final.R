## 06_Integration_Final.R

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SUBTYPE <- if (length(args) >= 1) args[1] else "ER"

nSample <- nrow(read_subtype_targets(SUBTYPE))

so <- readRDS(file.path(RDS_DIR, paste0("so_", SUBTYPE, "_all.RDS")))

## ---- Level 2 labels from the three compartments ----
meta_col <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
              "db_score", "db_type", "S.Score", "G2M.Score", "Phase",
              "seurat_clusters", "infercnv_score", "individual_clusters",
              "clusters_L1", "cell_type_L1", "clusters_L2", "cell_type_L2")

meta_L2 <- do.call(rbind, unname(lapply(compartments, function(cfg) {
  readRDS(file.path(RDS_DIR, paste0("so_", SUBTYPE, "_", cfg$suffix, ".RDS")))@meta.data[, meta_col]
})))

so <- so[, rownames(meta_L2)]
so@meta.data <- meta_L2

## ---- Level 3 clusters from the sub-populations ----
so@meta.data$clusters_L3 <- as.character(so@meta.data$clusters_L2)
for (cfg in subpopulations) {
  so_sub <- readRDS(file.path(RDS_DIR, paste0("so_", SUBTYPE, "_", cfg$suffix, ".RDS")))
  m <- match(colnames(so_sub), colnames(so))
  so@meta.data$clusters_L3[m] <- as.character(so_sub@meta.data$clusters_L3)
}

ct_names_L2 <- names(cols_Somi)[names(cols_Somi) %in% so@meta.data$cell_type_L2]
so@meta.data$cell_type_L2 <- factor(so@meta.data$cell_type_L2, levels = ct_names_L2)

save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_L2",
                    raster = FALSE, cols = cols_Somi[ct_names_L2], pt.size = 0.2,
                    shuffle = TRUE),
            paste0(SUBTYPE, "_final_UMAP_celltype"), width = 8, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_L2",
                    raster = FALSE, cols = cols_Somi[ct_names_L2], pt.size = 0.2,
                    split.by = "orig.ident", ncol = 3),
            paste0(SUBTYPE, "_final_UMAP_celltype_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)
save_figure(composition_barplot(so, "cell_type_L2", cols_Somi[ct_names_L2]),
            paste0(SUBTYPE, "_final_Barplot_celltype"), width = 8, height = 5)

saveRDS(so, file.path(RDS_DIR, paste0("so_", SUBTYPE, "_final.RDS")))

## ---- Relabel the tumour populations ----
cell_type_clean <- as.character(so@meta.data$cell_type_L2)
is_tumour <- so@meta.data$cell_type_L1 %in% Tumour_CT

for (label in names(tumour_relabel[[SUBTYPE]])) {
  sel <- is_tumour & so@meta.data$clusters_L1 %in% tumour_relabel[[SUBTYPE]][[label]]
  cell_type_clean[sel] <- label
}

if (SUBTYPE %in% tumour_rest_contamination) {
  cell_type_clean[cell_type_clean == "Tumor"] <- "Contamination"
}

## ---- Remove the contaminating clusters ----
contam <- read_contamination()
contam <- contam[contam$Subtype == SUBTYPE, ]

contam_level2 <- list(Immune = Immune_CT, Stroma = Stroma_CT, Tumour = Tumour_CT)
contam_level3 <- list(Immune_T = Immune_T_CT,
                      Immune_B_Plasma = Immune_B_Plasma_CT,
                      Immune_Myeloid = Immune_Myeloid_CT,
                      Stroma_Fibroblast = Stroma_Fibroblast_CT)

for (step in names(contam_level2)) {
  clusters <- contam$Cluster[contam$Step == step]
  if (length(clusters) == 0) next
  sel <- so@meta.data$cell_type_L1 %in% contam_level2[[step]] &
         so@meta.data$clusters_L2 %in% clusters
  cell_type_clean[sel] <- "Contamination"
}
for (step in names(contam_level3)) {
  clusters <- contam$Cluster[contam$Step == step]
  if (length(clusters) == 0) next
  sel <- so@meta.data$cell_type_L2 %in% contam_level3[[step]] &
         so@meta.data$clusters_L3 %in% clusters
  cell_type_clean[sel] <- "Contamination"
}

keep_clean <- cell_type_clean != "Contamination"

so@graphs$RNA_nn <- so@graphs$RNA_snn <- NULL
so <- so[, keep_clean]
cell_type_clean <- cell_type_clean[keep_clean]

ct_names_clean <- names(cols_Somi)[names(cols_Somi) %in% cell_type_clean]
so@meta.data$cell_type_clean <- factor(cell_type_clean, levels = ct_names_clean)
print(t(t(table(so@meta.data$cell_type_clean))))

save_figure(DimPlot(so, reduction = "umap", group.by = "clusters_L1",
                    raster = FALSE, cols = col.p2, pt.size = 0.2, shuffle = TRUE),
            paste0(SUBTYPE, "_clean_UMAP_cluster"), width = 6.5, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_clean",
                    raster = FALSE, cols = cols_Somi[ct_names_clean], pt.size = 0.2,
                    shuffle = TRUE),
            paste0(SUBTYPE, "_clean_UMAP_celltype"), width = 8, height = 5.5)
save_figure(DimPlot(so, reduction = "umap", group.by = "cell_type_clean",
                    raster = FALSE, cols = cols_Somi[ct_names_clean], pt.size = 0.2,
                    split.by = "orig.ident", ncol = 3),
            paste0(SUBTYPE, "_clean_UMAP_celltype_bySample"),
            width = 18, height = ceiling(nSample / 3) * 5)
save_figure(composition_barplot(so, "cell_type_clean", cols_Somi[ct_names_clean]),
            paste0(SUBTYPE, "_clean_Barplot_celltype"), width = 8, height = 5)

saveRDS(so, file.path(RDS_DIR, paste0("so_", SUBTYPE, "_clean.RDS")))

## ---- Markers and pathways of the tumour clusters ----
sub_clusters <- tumour_subset_clusters[[SUBTYPE]]
if (!is.null(sub_clusters)) {
  so_sub <- so[, so@meta.data$clusters_L1 %in% sub_clusters]
  markers <- FindAllMarkers(so_sub, logfc.threshold = 1, group.by = "clusters_L1",
                            min.pct = 0.1, only.pos = TRUE)
  markers <- markers[markers$p_val_adj < 0.05, ]
  saveRDS(markers, file.path(GENELIST_DIR, paste0(SUBTYPE, "_tumor_subset_marker.RDS")))
  write.csv(markers, file.path(GENELIST_DIR, paste0(SUBTYPE, "_tumor_subset_marker.csv")))

  save_go_kegg(go_kegg_by_cluster(markers, n_genes = 10),
               paste0(SUBTYPE, "_tumor_subset"))
}
