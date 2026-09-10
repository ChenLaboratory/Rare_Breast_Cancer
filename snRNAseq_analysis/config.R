## config.R

if (!exists("cols_Somi")) source("setup.R")

set.seed(42)

CELLRANGER_DIR <- "CellRanger"
METADATA_DIR   <- "metadata"
RDS_DIR        <- "RDS"
INFERCNV_DIR   <- "infercnv_output"
GENELIST_DIR   <- "GeneList"
FIGURE_DIR     <- "Figure"

## ---- Sample metadata ----
read_targets <- function() {
  read.delim(file.path(METADATA_DIR, "targets_all.txt"), stringsAsFactors = FALSE)
}

targets_for_sample <- function(sample_id) {
  tg <- read_targets()
  i <- match(sample_id, tg$Sample)
  if (is.na(i)) stop("Sample not listed in targets_all.txt: ", sample_id)
  tg[i, ]
}

read_subtype_targets <- function(subtype) {
  read.delim(file.path(METADATA_DIR, paste0("targets_", subtype, ".txt")),
             stringsAsFactors = FALSE)
}

## Manual cluster to cell type assignment. `step` is All, Immune, Stroma or Tumour
read_celltype_annotation <- function(subtype, step) {
  read.delim(file.path(METADATA_DIR, paste0(subtype, "_celltype_", step, ".txt")),
             stringsAsFactors = FALSE)
}

read_contamination <- function() {
  read.delim(file.path(METADATA_DIR, "Contamination.txt"), stringsAsFactors = FALSE)
}

read_marker_panel <- function(cell_types = NULL) {
  GOI <- read.delim(file.path(METADATA_DIR, "HsMarkers_snRNAseq.txt"),
                    stringsAsFactors = FALSE)
  if (!is.null(cell_types)) GOI <- GOI[GOI$Cell_Type %in% cell_types, ]
  unique(GOI$Genes)
}

## CellRanger per sample filtered matrix.
cellranger_h5 <- function(sample_id, run) {
  parent <- paste0("somiK_", run)
  file.path(CELLRANGER_DIR, parent, paste0(run, "_Combine"), "outs",
            "per_sample_outs", sample_id, "count",
            "sample_filtered_feature_bc_matrix.h5")
}

## ---- Subtypes ----
subtypes_all <- c("ER", "MED", "IMPC", "TNBC", "CYS", "APO", "PLC", "MpBC")
subtypes_tumour <- c("ER", "MED", "IMPC", "TNBC", "CYS", "APO", "PLC")
subtypes_caf    <- c("TNBC", "MED", "CYS", "APO", "ER", "IMPC", "PLC")

## ---- Clustering resolutions ----
res_all <- c(ER = 0.2, MED = 0.2, IMPC = 0.2, TNBC = 0.25,
             APO = 0.2, CYS = 0.2, PLC = 0.2, MpBC = 0.5)

res_compartment <- list(
  Immune = c(ER = 0.4, MED = 0.4, IMPC = 0.6, TNBC = 0.4,
             APO = 0.8, CYS = 0.4, PLC = 0.5, MpBC = 0.5),
  Stroma = c(ER = 0.2, MED = 0.3, IMPC = 0.2, TNBC = 0.2,
             APO = 0.2, CYS = 0.2, PLC = 0.2, MpBC = 0.5),
  Tumour = c(ER = 0.2, MED = 0.3, IMPC = 0.2, TNBC = 0.2,
             APO = 0.2, CYS = 0.2, PLC = 0.2, MpBC = 0.5)
)

res_subpop <- list(
  Myeloid    = c(ER = 0.4, MED = 0.4, IMPC = 0.4, TNBC = 0.4,
                 APO = 0.4, CYS = 0.4, PLC = 0.4, MpBC = 0.4),
  T          = c(ER = 0.4, MED = 0.6, IMPC = 0.4, TNBC = 0.4,
                 APO = 0.4, CYS = 0.2, PLC = 0.4, MpBC = 0.4),
  B_Plasma   = c(ER = 0.4, MED = 0.2, IMPC = 0.2, TNBC = 0.2,
                 APO = 0.2, CYS = 0.2, PLC = 0.2, MpBC = 0.2),
  Fibroblast = c(ER = 0.2, MED = 0.2, IMPC = 0.4, TNBC = 0.4,
                 APO = 0.2, CYS = 0.4, PLC = 0.4, MpBC = 0.05)
)

## ---- Analysis hierarchy ----
## Level 2: the three compartments taken out of the all cell integration.
compartments <- list(
  Immune = list(cell_types = Immune_CT, suffix = "immune"),
  Stroma = list(cell_types = Stroma_CT, suffix = "stroma"),
  Tumour = list(cell_types = Tumour_CT, suffix = "tumour")
)

## Level 3: the sub-populations taken out of a compartment.
subpopulations <- list(
  Myeloid = list(parent = "immune", cell_types = Immune_Myeloid_CT,
                 suffix = "immune_myeloid",
                 marker_types = c("Macrophage", "TAM", "Monocyte")),
  T = list(parent = "immune", cell_types = Immune_T_CT,
           suffix = "immune_t",
           marker_types = c("CD4_T", "CD8_T", "Treg", "Ex_T", "NK")),
  B_Plasma = list(parent = "immune", cell_types = Immune_B_Plasma_CT,
                  suffix = "immune_b_plasma",
                  marker_types = c("Naive_B", "Memory_B", "Plasma")),
  Fibroblast = list(parent = "stroma", cell_types = Stroma_Fibroblast_CT,
                    suffix = "stroma_fibroblast",
                    marker_types = c("Fibroblast", "iCAF", "myCAF", "apCAF"))
)

## Manual relabelling of the tumour compartment
tumour_relabel <- list(
  MED  = list(Tumor_1 = 0, Tumor_2 = 9, Epithelial = 11),
  TNBC = list(Tumor_1 = 0, Tumor_2 = 7, Tumor_3 = 10, Tumor_4 = 11,
              Cycling_Tumor = 5),
  ER   = list(Tumor_1 = 0, Tumor_2 = 4, Cycling_Tumor = 5),
  IMPC = list(Tumor_1 = 0, Tumor_2 = 4, Cycling_Tumor = 5),
  PLC  = list(Tumor_1 = 0, Tumor_2 = 4, Cycling_Tumor = 5, Tumor_3 = 8,
              Tumor_4 = c(11, 14), Tumor_5 = 15, Tumor_6 = 17),
  MpBC = list(Tumor_Mes = c(2, 7, 11, 12, 13, 15, 22, 26),
              Tumor_Epi = c(4, 5, 6), Cycling_Tumor = 3,
              Epithelial = c(19, 20))
)

## Subtypes where some tumour cells are treated as contamination.
tumour_rest_contamination <- c("TNBC", "ER", "IMPC", "PLC")
tumour_subset_clusters <- list(
  IMPC = c(0, 4, 5),
  MED  = c(0, 9),
  TNBC = c(0, 5, 7, 10, 11),
  ER   = c(0, 4, 5),
  PLC  = c(0, 4, 5, 8, 11, 14, 15, 17)
)

## ---- inferCNV ----
INFERCNV_REF <- "SK49"
GENE_ORDER_FILE <- file.path(METADATA_DIR, "gene_ordering_table.tsv")

## ---- Doublet detection ----
dbr_per_1k <- function(run) if (run %in% c("SK17", "SK18", "SK19")) 0.004 else 0.008
dbl_split_size <- 3e4

## ---- Pseudo-bulk ----
## Tumour pseudo-bulk profiles dropped after manual review.
pb_tumour_exclude <- list(
  list(subtype = "APO", cluster = "Tumor_3"),
  list(subtype = "ER",  cluster = "Tumor_2")
)

pb_tumour_lib_cutoff <- 2e5
caf_treat_lfc <- c(CYS = log2(1.6), MED = log2(3), TNBC = log2(2))
heatmap_top_n <- c(40, 30)
CAF_SIGNATURE_FILE <- file.path(METADATA_DIR, "CAF_Signatures.txt")

## ---- Helper: save a figure ----
save_figure <- function(x, name, width = 7, height = 7, ncol = 2, res = 300) {
  file <- file.path(FIGURE_DIR, paste0(name, ".png"))
  png(file, width = width, height = height, units = "in", res = res, pointsize = 15)
  on.exit(dev.off(), add = TRUE)
  if (is.list(x) && !inherits(x, c("ggplot", "pheatmap"))) {
    do.call(gridExtra::grid.arrange, c(x, ncol = ncol))
  } else if (inherits(x, "ggplot") || inherits(x, "pheatmap")) {
    print(x)
  }
  invisible(file)
}

## ---- Helper: shared nearest neighbours and UMAP ----
run_umap_from_knn <- function(so, reduction = "pca", dims = 1:30, k = 30,
                              min_dist = 0.4, spread = 1.2, n_threads = 16,
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

## ---- Helper: integrate one set of cells ----
integrate_cells <- function(so, n_dims = 30, n_features = 5000, min_dist = 0.4) {
  so <- FindVariableFeatures(so, nfeatures = n_features, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = n_dims, verbose = FALSE)
  so <- IntegrateLayers(so, method = HarmonyIntegration,
                        orig.reduction = "pca",
                        new.reduction = "integrated.harmony",
                        dims = 1:n_dims, k.anchor = 20, verbose = FALSE)
  run_umap_from_knn(so, reduction = "integrated.harmony", dims = 1:n_dims,
                    min_dist = min_dist)
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

## ---- Helper: marker table to top genes ----
top_markers <- function(marker_table, n = 5) {
  by_cluster <- split(marker_table, marker_table$cluster)
  top <- do.call(rbind, lapply(by_cluster, head, n = n))
  unique(top$gene)
}

## ---- Helper: composition bar plot ----
composition_barplot <- function(so, group.by, cols, legend_title = "Celltype") {
  tab <- table(so@meta.data$orig.ident, so@meta.data[[group.by]])
  print(colSums(tab))
  print(t(tab))
  perc <- t(round(100 * tab / rowSums(tab), 2))
  dat <- data.frame(
    Sample = factor(rep(colnames(perc), each = nrow(perc)), levels = colnames(perc)),
    Group  = factor(rep(rownames(perc), ncol(perc)), levels = rownames(perc)),
    Proportion = as.vector(perc))
  ggplot(dat, aes(fill = Group, y = Proportion, x = Sample)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = cols, name = legend_title) +
    geom_bar(position = position_fill(reverse = TRUE), stat = "identity",
             color = "black", linewidth = 0.2) +
    guides(fill = guide_legend(ncol = 2)) +
    coord_flip() + theme_classic()
}

## ---- Helper: inferCNV score bar plot ----
infercnv_boxplot <- function(so, group.by = "seurat_clusters", title = "") {
  dat <- so@meta.data[, c("infercnv_score", group.by)]
  ggplot(dat, aes(x = .data[[group.by]], y = infercnv_score,
                  color = .data[[group.by]])) +
    geom_boxplot(show.legend = FALSE) +
    scale_color_manual(values = col.p2) +
    labs(x = "", y = "inferCNV instability score", title = title) +
    coord_cartesian(ylim = c(0, min(0.008, max(dat$infercnv_score, na.rm = TRUE)))) +
    theme_classic() + Seurat::RotatedAxis()
}

## ---- Helper: GO and KEGG analysis of cluster markers ----
go_kegg_by_cluster <- function(marker_table, n_genes = Inf) {
  by_cluster <- split(marker_table, marker_table$cluster)
  res_GO <- res_KEGG <- list()
  for (i in names(by_cluster)) {
    genes <- head(by_cluster[[i]]$gene, n_genes)
    entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes,
                                    keytype = "SYMBOL", column = "ENTREZID")
    entrez <- unique(entrez[!is.na(entrez)])
    if (length(entrez) == 0) next
    tpg <- limma::topGO(limma::goana(entrez, species = "Hs"), n = Inf, p.value = 0.05)
    tpk <- limma::topKEGG(limma::kegga(entrez, species = "Hs"), n = Inf, p.value = 0.05)
    tpg$cluster <- tpk$cluster <- i
    res_GO[[i]] <- tpg
    res_KEGG[[i]] <- tpk
  }
  list(GO = res_GO, KEGG = res_KEGG)
}

## Writes the result of go_kegg_by_cluster to GeneList/.
save_go_kegg <- function(res, prefix) {
  saveRDS(res$GO,   file.path(GENELIST_DIR, paste0(prefix, "_GO.RDS")))
  saveRDS(res$KEGG, file.path(GENELIST_DIR, paste0(prefix, "_KEGG.RDS")))
  write.csv(do.call(rbind, res$GO),
            file.path(GENELIST_DIR, paste0(prefix, "_GO.csv")))
  write.csv(do.call(rbind, res$KEGG),
            file.path(GENELIST_DIR, paste0(prefix, "_KEGG.csv")))
}

## ---- Helper: pseudo-bulk profiles ----
seurat2pb <- function(so, sample = "orig.ident", cluster = "seurat_clusters",
                      assay = "RNA") {
  if (inherits(so[[assay]], "Assay5") && length(Layers(so[[assay]])) > 1) {
    so <- JoinLayers(so)
  }
  edgeR::Seurat2PB(so, sample = sample, cluster = cluster, assay = assay)
}

## ---- Helper: add Entrez ids ----
add_entrez <- function(y) {
  GeneID <- AnnotationDbi::mapIds(org.Hs.eg.db, rownames(y),
                                  keytype = "SYMBOL", column = "ENTREZID")
  y$genes <- data.frame(GeneID = GeneID, Symbol = rownames(y))
  keep <- !is.na(y$genes$GeneID)
  print(table(keep))
  y[keep, ]
}

## ---- Helper: volcano plot ----
plot_volcano <- function(qlf, title = "") {
  df <- topTags(qlf, n = Inf)$table[, c("logFC", "FDR", "Symbol")]
  df$DE <- "NS"
  df$DE[df$FDR < 0.05 & df$logFC > 0] <- "Up"
  df$DE[df$FDR < 0.05 & df$logFC < 0] <- "Down"
  label_df <- do.call(rbind, lapply(c("Up", "Down"), function(d) {
    sub <- df[df$DE == d, ]
    head(sub[order(sub$FDR), ], 20)
  }))
  ggplot(df, aes(x = logFC, y = -log10(FDR), color = DE)) +
    geom_point(alpha = 1, size = 1) +
    theme_void() +
    scale_color_manual(values = c(Down = "blue", NS = "grey", Up = "red")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14, angle = 90),
          axis.text.x = element_text(size = 11, margin = margin(t = 3)),
          axis.text.y = element_text(size = 11, margin = margin(r = 3)),
          axis.ticks = element_line(linewidth = 0.4),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = "right", legend.title = element_blank(),
          legend.text = element_text(size = 14)) +
    guides(colour = guide_legend(override.aes = list(size = 2, shape = 16, alpha = 1))) +
    labs(title = title, x = expression(log[2] ~ FC), y = expression(-log[10] ~ FDR)) +
    ggrepel::geom_text_repel(data = label_df, aes(label = Symbol), size = 3,
                             max.overlaps = Inf, box.padding = 0.3,
                             point.padding = 0.2)
}

## ---- Helper: heatmap of top DE genes ----
de_heatmap <- function(qlf_set, lcpm, group, keep_subtypes, subtype_cols,
                       top = 40, exclusive = TRUE, up_only = TRUE, title = "") {
  top_up <- lapply(qlf_set, function(q) {
    tbl <- topTags(q, n = Inf)$table
    if (up_only) tbl <- tbl[tbl$logFC > 0, ]
    head(tbl$Symbol, top)
  })
  if (exclusive) {
    n_hit <- table(unlist(top_up, use.names = FALSE))
    only_once <- names(n_hit)[n_hit == 1]
    sig <- unlist(lapply(top_up, function(g) g[g %in% only_once]), use.names = FALSE)
  } else {
    sig <- unique(unlist(top_up, use.names = FALSE))
  }

  keep_cols <- group %in% keep_subtypes
  col_order <- order(match(as.character(group[keep_cols]), keep_subtypes))
  mat <- t(scale(t(lcpm[sig, keep_cols][, col_order])))

  annot <- data.frame(Subtype = factor(group[keep_cols][col_order],
                                       levels = keep_subtypes))
  rownames(annot) <- colnames(mat)

  pheatmap::pheatmap(mat,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-2, 2, length.out = 101), border_color = NA,
    cluster_rows = !exclusive, cluster_cols = !exclusive, scale = "none",
    fontsize_row = 8, show_colnames = TRUE,
    treeheight_row = 70, treeheight_col = 70, clustering_method = "ward.D2",
    main = title,
    annotation_col = annot,
    annotation_colors = list(Subtype = subtype_cols[keep_subtypes]),
    silent = TRUE)
}
