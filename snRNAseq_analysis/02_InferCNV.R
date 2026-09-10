## 02_InferCNV.R

source("setup.R")
source("config.R")
library(infercnv)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
SAMPLE_ID <- if (length(args) >= 1) args[1] else "SK01"

so     <- readRDS(file.path(RDS_DIR, paste0("so_", SAMPLE_ID, ".RDS")))
so_ref <- readRDS(file.path(RDS_DIR, paste0("so_", INFERCNV_REF, ".RDS")))

## ---- Gene ordering ----
geneTable <- read.table(GENE_ORDER_FILE, sep = "\t")[, -1]
geneTable <- geneTable[!duplicated(geneTable[, 1]), ]
rownames(geneTable) <- geneTable[, 1]
colnames(geneTable) <- c("symbol", "chr", "start", "end")

genes_used <- Reduce(intersect, list(rownames(geneTable), rownames(so), rownames(so_ref)))
message(SAMPLE_ID, ": ", length(genes_used), " / ", nrow(geneTable),
        " ordered genes present in both objects")
geneTable <- geneTable[genes_used, ]

## ---- Counts of the sample and of the reference ----
as_dge <- function(object, group, new_group) {
  y <- DGEList(counts = as.matrix(object[["RNA"]]$counts)[genes_used, ])
  y$samples <- data.frame(new_group = new_group, group = group,
                          row.names = colnames(y))
  y
}

clusters <- so@meta.data$seurat_clusters
dge     <- as_dge(so, SAMPLE_ID, as.character(clusters))
dge_ref <- as_dge(so_ref, INFERCNV_REF, "Normal")

dge_combine <- cbind(dge, dge_ref)
dge_combine$samples$new_group <- factor(dge_combine$samples$new_group,
                                        levels = c(levels(clusters), "Normal"))
dge_combine$samples$group <- factor(dge_combine$samples$group,
                                    levels = c(SAMPLE_ID, INFERCNV_REF))
dge_combine$genes <- geneTable

## ---- Run inferCNV ----
out_dir <- file.path(INFERCNV_DIR, SAMPLE_ID)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = dge_combine$counts,
  annotations_file  = dge_combine$samples,
  gene_order_file   = dge_combine$genes[, 2:4],
  ref_group_names   = "Normal")

infercnv::run(infercnv_obj,
              resume_mode = FALSE,
              cutoff = 0.1,
              out_dir = out_dir,
              cluster_by_groups = TRUE,
              denoise = TRUE,
              HMM = FALSE)

## ---- Instability score per nucleus ----
cnv <- readRDS(file.path(out_dir, "run.final.infercnv_obj"))

grouped <- c(cnv@observation_grouped_cell_indices, cnv@reference_grouped_cell_indices)
ordered_groups <- data.frame(
  group    = rep(names(grouped), lengths(grouped)),
  position = unlist(grouped, use.names = FALSE))
ordered_groups$group <- factor(ordered_groups$group,
                               levels = levels(dge_combine$samples$new_group))

expr <- cnv@expr.data[, ordered_groups$position]
ordered_groups$sample <- dge_combine$samples[colnames(expr), "group"]

instability_score <- apply((expr - 1)^2, 2, mean)
saveRDS(instability_score,
        file.path(RDS_DIR, paste0(SAMPLE_ID, "_instability_score.rds")))

dat <- cbind.data.frame(instability_score, ordered_groups)
save_figure(
  ggplot(dat, aes(x = group, y = instability_score, color = group)) +
    geom_boxplot(show.legend = FALSE) +
    scale_color_manual(values = col.p2) +
    labs(x = "", y = "inferCNV instability score", title = SAMPLE_ID) +
    coord_cartesian(ylim = c(0, min(0.008, max(instability_score)))) +
    theme_classic() + Seurat::RotatedAxis(),
  paste0(SAMPLE_ID, "_inferCNV_instability"), width = 7, height = 4)

## ---- Copy number heatmap ----
ordered_groups$group <- droplevels(ordered_groups$group)
group_cols <- setNames(col.p2[seq_len(nlevels(ordered_groups$group))],
                       levels(ordered_groups$group))
sample_cols <- setNames(c("#4d4d4d", "#D3D3D3")[seq_len(nlevels(ordered_groups$sample))],
                        levels(ordered_groups$sample))

ht_list_raw <- lapply(levels(ordered_groups$group), function(g) {
  keep <- ordered_groups$group == g
  ha <- HeatmapAnnotation(group = ordered_groups$group[keep],
                          sample = ordered_groups$sample[keep],
                          col = list(group = group_cols, sample = sample_cols),
                          show_legend = TRUE, show_annotation_name = FALSE)
  Heatmap(matrix = expr[, keep],
          width = 7,
          top_annotation = ha,
          cluster_columns = TRUE, show_column_dend = FALSE,
          show_parent_dend_line = FALSE, show_heatmap_legend = FALSE,
          column_split = ordered_groups$sample[keep],
          column_title = NULL, column_gap = unit(0, "mm"),
          cluster_rows = FALSE, row_split = cnv@gene_order$chr, row_title_rot = 0,
          use_raster = TRUE,
          show_column_names = FALSE, show_row_names = FALSE)
})
names(ht_list_raw) <- levels(ordered_groups$group)
ht_list <- Reduce(`+`, ht_list_raw)

save_figure(ComplexHeatmap::draw(ht_list, ht_gap = unit(1, "mm")),
            paste0(SAMPLE_ID, "_inferCNV_heatmap"), width = 10, height = 10)

saveRDS(ht_list,     file.path(RDS_DIR, paste0(SAMPLE_ID, "_heatmap_list_combined.rds")))
saveRDS(ht_list_raw, file.path(RDS_DIR, paste0(SAMPLE_ID, "_heatmap_list_raw.rds")))
