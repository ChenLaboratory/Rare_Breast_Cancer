## config.R
##

if (!exists("cols_Somi")) source("setup.R")

set.seed(2025)

## ---- Directories ----
PROSEG_DIR       <- "Proseg"
SEG_DIR          <- "Proseg_mask"
METADATA_DIR     <- "metadata"
RDS_DIR          <- "RDS"
RDS_POST_CAF_DIR <- "RDS_post_CAF"
GENELIST_DIR     <- "GeneList"
FIGURE_DIR       <- "Figure"

for (d in c(RDS_DIR, RDS_POST_CAF_DIR, GENELIST_DIR, FIGURE_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

stage_dir <- function(stage = c("main", "post_CAF")) {
  stage <- match.arg(stage)
  if (stage == "main") RDS_DIR else RDS_POST_CAF_DIR
}

## ---- Sample metadata ----
read_targets <- function() {
  read.delim(file.path(METADATA_DIR, "targets_all.txt"), stringsAsFactors = FALSE)
}

targets_for_sample <- function(sample_id) {
  tg <- read_targets()
  i <- match(sample_id, tg$SampleID)
}

read_subtype_targets <- function(subtype) {
  read.delim(file.path(METADATA_DIR, paste0("targets_xenium_", subtype, ".txt")),
             stringsAsFactors = FALSE)
}

## Manual cluster to cell type assignment, one row per sample and cluster.
read_celltype_annotation <- function(subtype) {
  read.delim(file.path(METADATA_DIR, paste0(subtype, "_celltype_Individual.txt")),
             stringsAsFactors = FALSE)
}

## Manual CAF cluster to CAF subtype assignment, used by the post CAF pass of script 02.
## One row per CAF cluster of script 05, with columns cluster and cell_type.
read_caf_annotation <- function(subtype) {
  read.delim(file.path(METADATA_DIR, paste0(subtype, "_CAF_celltype.txt")),
             stringsAsFactors = FALSE)
}

## Marker gene panel used for the signature heatmap of the grid level integration
read_marker_panel <- function() {
  read.delim(file.path(METADATA_DIR, "HsMarkers_Xenium.txt"),
             stringsAsFactors = FALSE)$Genes
}

## ---- Per subtype parameters ----
## Clustering resolution for the common niche clusters (script 04) and for the CAF clusters (script 05).
grid_res <- c(IMPC = 0.1, PLC = 0.1, CYS = 0.1, MED = 0.1,
              APO = 0.1, ER = 0.1, TNBC = 0.1)
caf_res <- c(IMPC = 0.3, PLC = 0.3, CYS = 0.3, MED = 0.3,
             APO = 0.3, ER = 0.3, TNBC = 0.3)

## Per sample niche clusters dropped before the grid level integration. 
## These are the artefactual niches seen in the individual niche analysis of script 02 (tissue edges, folds and detached pieces).
unwanted_niches <- list(
  IMPC_1 = 7, IMPC_3 = 8, IMPC_4 = 8, IMPC_10 = 5,
  ER_6 = 4, TNBC_6 = 6, MED_2 = 5, MED_5 = 6
)

## Minimum pseudo-bulk library size for a sample by niche profile to enter the heatmaps of script 04.
pb_lib_cutoff <- c(IMPC = 1e5, TNBC = 1e5, ER = 2e5, MED = 4e4)
pb_lib_cutoff_for <- function(subtype) {
  if (subtype %in% names(pb_lib_cutoff)) pb_lib_cutoff[[subtype]] else 5e4
}

## ---- Helper: save a figure ----
save_figure <- function(x, name, width = 7, height = 7, ncol = 2, res = 300) {
  file <- file.path(FIGURE_DIR, paste0(name, ".png"))
  png(file, width = width, height = height, units = "in", res = res, pointsize = 15)
  on.exit(dev.off(), add = TRUE)
  if (inherits(x, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(x)
  } else if (is.list(x) && !inherits(x, c("ggplot", "pheatmap"))) {
    do.call(gridExtra::grid.arrange, c(x, ncol = ncol))
  } else {
    print(x)
  }
  invisible(file)
}

## ---- Helper: shared nearest neighbours and UMAP ----
run_umap_from_knn <- function(so, reduction = "pca", dims = 1:30, k = 30,
                              min_dist = 0.7, spread = 1.2, n_threads = 16,
                              reduction.name = "umap", seed = 42) {
  so <- FindNeighbors(so, reduction = reduction, dims = dims, k.param = k,
                      annoy.metric = "cosine", compute.SNN = FALSE,
                      return.neighbor = TRUE)
  so <- FindNeighbors(so, reduction = reduction, dims = dims, k.param = k,
                      annoy.metric = "cosine", compute.SNN = TRUE,
                      return.neighbor = FALSE)
  nbrs <- so[[paste0(DefaultAssay(so), ".nn")]]
  nn <- list(idx = Indices(nbrs), dist = Distances(nbrs))
  emb <- uwot::umap(X = NULL, nn_method = nn, n_neighbors = k,
                    min_dist = min_dist, spread = spread,
                    n_sgd_threads = n_threads, verbose = TRUE, seed = seed)
  key <- Key(object = reduction.name, quiet = TRUE)
  colnames(emb) <- paste0(key, seq_len(ncol(emb)))
  rownames(emb) <- rownames(nn$idx)
  so[[reduction.name]] <- CreateDimReducObject(embeddings = emb, key = key,
      assay = DefaultAssay(so), global = TRUE)
  so
}

## ---- Helper: renumber cluster labels ----
renumber_clusters <- function(labels) {
  labels <- as.character(labels)
  sp <- do.call(rbind, strsplit(labels, "_", fixed = TRUE))
  main <- as.integer(sp[, 1])
  sub <- if (ncol(sp) == 1) rep(-1L, length(labels)) else {
    ifelse(is.na(sp[, 2]), -1L, as.integer(sp[, 2]))
  }
  lev <- unique(labels[order(main, sub)])
  factor(labels, levels = lev, labels = seq_along(lev) - 1)
}

## ---- Helper: cell types present ----
present_or_default <- function(candidates, ct_names, default) {
  found <- intersect(ct_names, candidates)
  if (length(found) > 0) found else default
}

## ---- Helper: marker table to top genes ----
top_markers <- function(marker_table, n = 5) {
  by_cluster <- split(marker_table, marker_table$cluster)
  top <- do.call(rbind, lapply(by_cluster, head, n = n))
  unique(top$gene)
}

## ---- Helper: pseudo-bulk profiles ----
seurat2pb <- function(so, sample = "orig.ident", cluster = "seurat_clusters",
                      assay = "RNA") {
  if (inherits(so[[assay]], "Assay5") && length(Layers(so[[assay]])) > 1) {
    so <- JoinLayers(so)
  }
  edgeR::Seurat2PB(so, sample = sample, cluster = cluster, assay = assay)
}
