## 04_Integration_Grid.R
##
## Run once per subtype. The subtype and the stage can be given on the command line,
## or edited below when running interactively:
##   Rscript 04_Integration_Grid.R IMPC
##   Rscript 04_Integration_Grid.R IMPC post_CAF
##
## STAGE = "main"     reads the bins written by the main pass of script 02.
## STAGE = "post_CAF" reads the bins written by the post CAF pass of script 02

source("setup.R")
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
SUBTYPE <- if (length(args) >= 1) args[1] else "IMPC"
STAGE   <- if (length(args) >= 2) args[2] else "main"
stopifnot(STAGE %in% c("main", "post_CAF"))

targets <- read_subtype_targets(SUBTYPE)
nSample <- nrow(targets)
print(targets[, 1:3])

obj_dir <- stage_dir(STAGE)
resolution <- grid_res[[SUBTYPE]]
fig <- function(name) paste0(SUBTYPE, "_", name, if (STAGE == "post_CAF") "_post_CAF" else "")

n_dims_pca <- 30
n_dims_harmony <- 20
n_top <- 5

## ---- Read the per sample bins ----
meta_col <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "seurat_clusters")
counts_all <- meta_all <- spe_gr_list <- list()

for (i in seq_len(nSample)) {
  sample_id <- targets$SampleID[i]
  so_i  <- readRDS(file.path(obj_dir, paste0("so_gr_", sample_id, ".RDS")))
  spe_i <- readRDS(file.path(obj_dir, paste0("spe_gr_", sample_id, ".RDS")))

  if (sample_id %in% names(unwanted_niches)) {
    kp <- !so_i@meta.data$seurat_clusters %in% unwanted_niches[[sample_id]]
    so_i  <- so_i[, kp]
    spe_i <- spe_i[, kp]
  }

  so_i@meta.data$orig.ident <- factor(sample_id)
  counts_all[[i]] <- so_i@assays$RNA$counts
  meta_all[[i]] <- so_i@meta.data[, meta_col]
  rownames(meta_all[[i]]) <- colnames(counts_all[[i]]) <-
    paste(sample_id, rownames(meta_all[[i]]), sep = "-")
  spe_gr_list[[i]] <- spe_i
}
names(counts_all) <- names(meta_all) <- names(spe_gr_list) <- targets$SampleID

## ---- Pool the bins ----
so <- CreateSeuratObject(counts = do.call(cbind, counts_all),
                         meta.data = do.call(rbind, meta_all),
                         assay = "RNA", min.cells = 0, min.features = 0)
so@meta.data$individual_niche <- so@meta.data$seurat_clusters
so <- NormalizeData(so, scale.factor = 100, verbose = FALSE)
so[["RNA"]] <- split(so[["RNA"]], f = so$orig.ident)
rm(counts_all)

## ---- Harmony integration ----
VariableFeatures(so) <- rownames(so)
so <- ScaleData(so, features = rownames(so), verbose = FALSE)
so <- RunPCA(so, features = rownames(so), npcs = n_dims_pca, verbose = FALSE)
so <- IntegrateLayers(so, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "integrated.harmony",
    features = rownames(so),
    .options = harmony::harmony_options(
      orig       = "pca",
      npcs       = n_dims_harmony,
      max_iter   = 50,
      epsilon    = 1e-4,
      early_stop = TRUE),
    dims = 1:n_dims_harmony, k.anchor = 20, verbose = FALSE)

so <- run_umap_from_knn(so, reduction = "integrated.harmony",
                        dims = 1:n_dims_harmony, min_dist = 0.1)

save_figure(DimPlot(so, reduction = "umap", group.by = "orig.ident",
                    raster = FALSE, cols = col.p2, pt.size = 0.2, shuffle = TRUE),
            fig("UMAP_sample"), width = 6.5, height = 5.5)

## ---- Common niche clustering ----
so <- FindClusters(so, resolution = resolution, verbose = FALSE)
print(table(so@meta.data$seurat_clusters))

save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.nclst,
                    pt.size = 1),
            fig("UMAP_CmNch"), width = 9, height = 7.5)
save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.nclst,
                    pt.size = 0.8, split.by = "orig.ident", ncol = 3),
            fig("UMAP_CmNch_bySample"), width = 18, height = ceiling(nSample / 3) * 5)

## ---- Pseudo-bulk profiles of each sample by common niche ----
so <- JoinLayers(so)
y <- seurat2pb(so, sample = "orig.ident")
colnames(y) <- gsub("_cluster", "-CmNch_", colnames(y))

keep_pb <- y$samples$lib.size > pb_lib_cutoff_for(SUBTYPE)
print(table(keep_pb))
y <- y[, keep_pb]

lcpm <- edgeR::cpm(y, log = TRUE)
rownames(lcpm) <- y$genes$gene

grp <- factor(y$samples$sample, levels = targets$SampleID)
NicheClst <- factor(y$samples$cluster, levels = levels(so@meta.data$seurat_clusters))

annot_col <- data.frame(Sample = grp, NicheClst = NicheClst)
rownames(annot_col) <- colnames(y)
ann_colors <- list(Sample = col.p2[seq_len(nlevels(grp))],
                   NicheClst = col.nclst[seq_len(nlevels(NicheClst))])
names(ann_colors$Sample) <- levels(grp)
names(ann_colors$NicheClst) <- levels(NicheClst)

ha <- ComplexHeatmap::HeatmapAnnotation(df = annot_col, col = ann_colors,
    annotation_name_gp = grid::gpar(fontsize = 10))
col_fun <- circlize::colorRamp2(seq(-2, 2, length.out = 101),
    colorRampPalette(col.spec[10:1])(101))

heatmap_lcpm <- function(genes, name = "Expr") {
  ComplexHeatmap::Heatmap(t(scale(t(lcpm[genes, ]))),
    col = col_fun, name = name,
    column_split = NicheClst,
    cluster_columns = TRUE, cluster_rows = TRUE,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    show_column_names = TRUE, top_annotation = ha,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_dend_width = grid::unit(50, "pt"),
    column_dend_height = grid::unit(50, "pt"))
}

## ---- Canonical cell type signatures across the common niches ----
markers <- read_marker_panel()
markers <- unique(markers[markers %in% rownames(so)])
save_figure(heatmap_lcpm(markers), fig("Heatmap_signature"), width = 10, height = 19)

## ---- Common niches on the tissue ----
plist <- list()
for (i in seq_len(nSample)) {
  sample_id <- targets$SampleID[i]
  in_sample <- sub("-.*$", "", colnames(so)) == sample_id
  spe_gr_list[[i]]@colData$common_niche <- so@meta.data$seurat_clusters[in_sample]
  present <- table(spe_gr_list[[i]]$common_niche) != 0
  plist[[i]] <- plotGrid(spe_gr_list[[i]], group.by = "common_niche",
                         cols = col.nclst[present], pol.border = FALSE) +
                ggtitle(sample_id)
}
save_figure(plist, fig("Spatial_CmNch"), width = 18,
            height = ceiling(nSample / 3) * 7, ncol = 3)

saveRDS(spe_gr_list, file.path(obj_dir, paste0("spe_gr_list_", SUBTYPE, ".RDS")))

## ---- Cell type composition of the common niches ----
ct_clean_all <- unique(unlist(lapply(spe_gr_list,
    function(s) colnames(colData(s)$cell_count))))
ct_names <- names(cols_Somi)[janitor::make_clean_names(names(cols_Somi)) %in%
                             setdiff(ct_clean_all, "overall")]

ccnt <- list()
for (i in seq_along(spe_gr_list)) {
  m <- rowsum(as.matrix(colData(spe_gr_list[[i]])$cell_count),
              group = colData(spe_gr_list[[i]])$common_niche)
  m <- m[, colnames(m) != "overall", drop = FALSE]
  rownames(m) <- paste(names(spe_gr_list)[i], rownames(m), sep = "-")
  ccnt[[i]] <- m
}
names(ccnt) <- names(spe_gr_list)

RN <- unlist(lapply(ccnt, rownames))
ccnt_all <- matrix(0, length(RN), length(ct_names), dimnames = list(RN, ct_names))
for (i in seq_along(ccnt)) {
  m1 <- match(rownames(ccnt[[i]]), RN)
  m2 <- match(colnames(ccnt[[i]]), janitor::make_clean_names(ct_names))
  ccnt_all[m1, m2] <- ccnt[[i]]
}

ccnt_all <- ccnt_all[rownames(ccnt_all) %in%
                     paste(y$samples$sample, y$samples$cluster, sep = "-"), ]
ct_names_noMixed <- ct_names[ct_names != "Mixed"]
ccnt_all <- ccnt_all[, ct_names_noMixed, drop = FALSE]

perc <- t(ccnt_all / rowSums(ccnt_all)) * 100

dat <- data.frame(
  Niche      = factor(rep(rownames(ccnt_all), each = length(ct_names_noMixed)),
                      levels = rownames(ccnt_all)),
  Celltype   = factor(rep(ct_names_noMixed, ncol(perc)), levels = ct_names_noMixed),
  Proportion = as.vector(perc),
  Count      = as.vector(t(ccnt_all)))
dat$CommonNiche <- factor(paste0("CmNch_", gsub("^.*-", "", dat$Niche)))

write.csv(t(ccnt_all), file.path(GENELIST_DIR, paste0(SUBTYPE, "_cell_count_CmNch.csv")))
write.csv(perc, file.path(GENELIST_DIR, paste0(SUBTYPE, "_percentage_CmNch.csv")))

composition_barplot <- function(y_var) {
  p <- ggplot(dat, aes(fill = Celltype, y = .data[[y_var]], x = Niche)) +
    scale_fill_manual(values = cols_Somi[ct_names_noMixed]) +
    guides(fill = guide_legend(ncol = 2)) +
    coord_flip() +
    facet_grid(CommonNiche ~ ., scales = "free_y", space = "free_y") +
    theme_classic()
  if (y_var == "Proportion") {
    p <- p + scale_y_continuous(labels = scales::percent) +
      geom_bar(position = position_fill(reverse = TRUE), stat = "identity",
               color = "black", linewidth = 0.2)
  } else {
    p <- p + geom_bar(stat = "identity", color = "black", linewidth = 0.2)
  }
  p
}
save_figure(composition_barplot("Proportion"), fig("Barplot_CmNch_cellProp"),
            width = 9, height = 9)
save_figure(composition_barplot("Count"), fig("Barplot_CmNch_cellCount"),
            width = 9, height = 9)

## ---- Cell type enrichment across the common niches ----
grp_2 <- factor(gsub("-.*$", "", rownames(ccnt_all)), levels = targets$SampleID)
NicheClst_2 <- factor(gsub("^.*-", "", rownames(ccnt_all)),
                      levels = levels(so@meta.data$seurat_clusters))

annot_col_2 <- data.frame(Sample = grp_2, NicheClst = NicheClst_2)
rownames(annot_col_2) <- rownames(ccnt_all)
ann_colors_2 <- list(Sample = col.p2[seq_len(nlevels(grp_2))],
                     NicheClst = col.nclst[seq_len(nlevels(NicheClst_2))])
names(ann_colors_2$Sample) <- levels(grp_2)
names(ann_colors_2$NicheClst) <- levels(NicheClst_2)

sample_total <- vapply(ccnt, sum, numeric(1))
scale_factor <- exp(mean(log(sample_total))) / sample_total[as.character(grp_2)]
mat_2 <- t(scale(log1p(ccnt_all * scale_factor)))

ha_2 <- ComplexHeatmap::HeatmapAnnotation(df = annot_col_2, col = ann_colors_2,
    annotation_name_gp = grid::gpar(fontsize = 10))
save_figure(ComplexHeatmap::Heatmap(mat_2,
    col = col_fun, name = "Enrichment",
    column_split = NicheClst_2,
    cluster_columns = TRUE, cluster_rows = TRUE,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    show_column_names = TRUE, top_annotation = ha_2,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_dend_width = grid::unit(50, "pt"),
    column_dend_height = grid::unit(50, "pt")),
  fig("Heatmap_celltype_enrichment"), width = 10, height = 9)

## ---- Common niche markers ----
marker_all <- FindAllMarkers(so, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
marker_all <- marker_all[marker_all$p_val_adj < 0.05, ]
saveRDS(marker_all, file.path(GENELIST_DIR, paste0(SUBTYPE, "_all_marker.RDS")))
write.csv(marker_all, file.path(GENELIST_DIR, paste0(SUBTYPE, "_all_marker.csv")))

top <- top_markers(marker_all, n = n_top)
save_figure(DotPlot(so, features = top, dot.scale = 8) + coord_flip(),
            fig("DotPlot_topMarkers"), width = 8, height = 13)
save_figure(heatmap_lcpm(top), fig("Heatmap_topMarkers"), width = 10, height = 19)

saveRDS(so, file.path(obj_dir, paste0("so_gr_", SUBTYPE, "_all.RDS")))
