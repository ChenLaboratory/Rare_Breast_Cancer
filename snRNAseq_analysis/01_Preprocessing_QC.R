## 01_Preprocessing_QC.R

source("setup.R")
source("config.R")
library(scDblFinder)

args <- commandArgs(trailingOnly = TRUE)
SAMPLE_ID <- if (length(args) >= 1) args[1] else "SK01"

tg <- targets_for_sample(SAMPLE_ID)
n_dims <- 30

## ---- Read the CellRanger output ----
counts <- Read10X_h5(filename = cellranger_h5(SAMPLE_ID, tg$Run))
colnames(counts) <- paste(SAMPLE_ID, colnames(counts), sep = "-")
so <- CreateSeuratObject(counts = counts, project = SAMPLE_ID, assay = "RNA",
                         min.cells = 0, min.features = 0)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
print(dim(so))

## ---- Doublet detection ----
n_cells <- ncol(so)
dbr <- n_cells * dbr_per_1k(tg$Run) / 1e3

if (n_cells < dbl_split_size) {
  sce <- scDblFinder::scDblFinder(as.SingleCellExperiment(so), dbr = dbr)
  so$db_score <- sce$scDblFinder.score
  so$db_type <- sce$scDblFinder.class
} else {
  set1 <- sort(sample(seq_len(n_cells), size = floor(n_cells / 2)))
  set2 <- setdiff(seq_len(n_cells), set1)
  so$db_score <- NA_real_
  so$db_type <- NA_character_
  for (part in list(set1, set2)) {
    sce <- scDblFinder::scDblFinder(as.SingleCellExperiment(so[, part]), dbr = dbr)
    so$db_score[part] <- sce$scDblFinder.score
    so$db_type[part] <- as.character(sce$scDblFinder.class)
  }
}
so$db_type <- factor(so$db_type, levels = c("singlet", "doublet"))
print(table(so$db_type))

## ---- Cell filtering ----
vcut <- quantile(so@meta.data$nCount_RNA, 0.95)
hcut <- quantile(so@meta.data$nFeature_RNA, 0.95)

keep_cell <- so@meta.data$nFeature_RNA > 500 &
             so@meta.data$percent.mt < tg$MT &
             so@meta.data$nFeature_RNA < hcut &
             so@meta.data$nCount_RNA < vcut &
             so@meta.data$db_type == "singlet"
message(SAMPLE_ID, ": ", sum(keep_cell), " / ", ncol(so), " nuclei pass filtering")
so <- so[, keep_cell]

## ---- Normalisation, PCA and UMAP ----
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 5000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = n_dims, verbose = FALSE)
so <- run_umap_from_knn(so, dims = 1:n_dims, min_dist = 0.7)

## Cell cycle scores are carried through the integration steps.
so <- CellCycleScoring(so, s.features = cc.genes$s.genes,
                       g2m.features = cc.genes$g2m.genes)

## ---- Cell clustering ----
so <- FindClusters(so, resolution = tg$Resolution, verbose = FALSE)
print(table(so@meta.data$seurat_clusters))

save_figure(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2,
                    pt.size = 0.2),
            paste0(SAMPLE_ID, "_UMAP_cluster"), width = 9, height = 7.5)

saveRDS(so, file.path(RDS_DIR, paste0("so_", SAMPLE_ID, ".RDS")))
