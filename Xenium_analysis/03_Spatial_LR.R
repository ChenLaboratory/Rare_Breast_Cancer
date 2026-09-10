## 03_Spatial_LR.R
##
## Run once per sample:
##   Rscript 03_Spatial_LR.R APO_3

source("setup.R")
source("config.R")
library(CellChat)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
SAMPLE_ID <- if (length(args) >= 1) args[1] else "APO_3"

tg <- targets_for_sample(SAMPLE_ID)
spe_gr <- readRDS(file.path(RDS_DIR, paste0("spe_gr_", SAMPLE_ID, ".RDS")))

n_perm <- 1e5 - 1
n_top_lr <- 24

## ---- Ligand receptor pairs covered by the Xenium panel ----
LR_db <- CellChat::CellChatDB.human$interaction
on_panel <- LR_db$ligand.symbol %in% rownames(spe_gr) &
            LR_db$receptor.symbol %in% rownames(spe_gr)
LR <- LR_db[on_panel, c("ligand.symbol", "receptor.symbol")]
message(SAMPLE_ID, ": ", nrow(LR), " ligand receptor pairs on the panel")

## ---- Local spatial co-expression of each pair ----
logcounts <- log2(as.matrix(spe_gr@assays@data$counts) + 1)
spe_gr@assays@data$logcounts <- logcounts
spe_gr <- findNbrsGrid(spe_gr)

empty <- matrix(0, ncol(spe_gr), nrow(LR), dimnames = list(NULL, rownames(LR)))
res <- list(cluster = empty, pvalue = empty, lisa = empty)

for (i in seq_len(nrow(LR))) {
  cc <- localMoran(spe_gr,
                   data1 = logcounts[LR[i, 2], ], data2 = logcounts[LR[i, 1], ],
                   permutations = n_perm, hhonly = FALSE)
  res$cluster[, i] <- as.character(cc$cluster)
  res$pvalue[, i]  <- cc$p.value
  res$lisa[, i]    <- cc$lisa
}

## ---- Rank the pairs by the number of significant bins ----
ranking <- order(colSums(res$cluster == "High-High"), decreasing = TRUE)
res <- lapply(res, function(m) m[, ranking, drop = FALSE])
LR <- LR[ranking, ]

saveRDS(res, file.path(RDS_DIR, paste0("res_LR_", SAMPLE_ID, ".RDS")))

df <- data.frame(LR = factor(colnames(res$cluster), levels = colnames(res$cluster)),
                 nSig = colSums(res$cluster == "High-High"))
save_figure(ggplot(df, aes(x = LR, y = nSig)) + geom_point() + theme_classic() +
              labs(x = "LR pair", y = "# Significant bins") +
              theme(axis.text.x = element_text(angle = 60, hjust = 1)),
            paste0(SAMPLE_ID, "_LR_ranking"), width = 12, height = 6)

## ---- Significant regions of the top pairs ----
gradient_fn <- colorRampPalette(c("yellow", "red3"))

plot_lr_grid <- function(i) {
  sel <- res$cluster[, i] == "High-High"
  cols_gr <- rep("grey80", ncol(spe_gr))
  pval <- -log10(res$pvalue[sel, i])
  scaled <- (pval - min(pval)) / (max(pval) - min(pval))
  cols_gr[sel] <- gradient_fn(100)[as.integer(scaled * 99) + 1]
  ## Invisible points, drawn only to obtain a continuous colour bar.
  legend_dummy <- data.frame(x = 1:100, y = 1:100,
                             val = seq(min(pval), max(pval), length.out = 100))
  plotGrid(spe_gr, cols = cols_gr) + guides(fill = "none") +
    geom_point(data = legend_dummy, aes(x = x, y = y, color = val), alpha = 0) +
    scale_color_gradient(name = "-log10(pval)", low = "yellow", high = "red3") +
    guides(color = guide_colorbar(title.position = "top")) +
    theme(legend.position = "right") +
    ggtitle(colnames(res$cluster)[i])
}

top_lr <- seq_len(min(n_top_lr, nrow(LR)))
for (page in split(top_lr, ceiling(top_lr / 4))) {
  save_figure(lapply(page, plot_lr_grid),
              paste0(SAMPLE_ID, "_LR_spatial_", page[1]),
              width = tg$width * 1.4, height = tg$height * 1.4, ncol = 2)
}

## ---- Cell types carrying each interaction ----
ct_assays <- setdiff(grep("^counts_", assayNames(spe_gr), value = TRUE), "counts_mixed")
ct_labels <- names(cols_Somi)[match(sub("^counts_", "", ct_assays),
                                    janitor::make_clean_names(names(cols_Somi)))]

lr_signal <- function(i) {
  sig <- matrix(1 * (res$cluster[, i] == "High-High"), ncol = 1)
  cnt <- vapply(ct_assays, function(a) {
    as.vector(as.matrix(assay(spe_gr, a)[c(LR[i, 1], LR[i, 2]), , drop = FALSE]) %*% sig)
  }, numeric(2))
  colnames(cnt) <- ct_labels
  0.5 * log2(t(cnt[1, , drop = FALSE]) %*% cnt[2, , drop = FALSE] + 1)
}

heatmap_lr <- function(i) {
  mat <- lr_signal(i)
  ComplexHeatmap::Heatmap(mat, name = "Signal",
    col = circlize::colorRamp2(seq(1, max(mat, na.rm = TRUE), length.out = 100),
                               viridis::viridis(100)),
    cluster_columns = TRUE, cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete", clustering_method_columns = "complete",
    row_title = paste0(LR[i, 1], " (Ligand)"),
    column_title = paste0(LR[i, 2], " (Receptor)"),
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    row_dend_width = grid::unit(50, "pt"),
    column_dend_height = grid::unit(50, "pt"))
}

grab_heatmap <- function(ht) {
  grid::grid.grabExpr(ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                                           merge_legends = FALSE))
}

for (page in split(top_lr, ceiling(top_lr / 6))) {
  save_figure(lapply(page, function(i) grab_heatmap(heatmap_lr(i))),
              paste0(SAMPLE_ID, "_LR_celltype_", page[1]),
              width = 13, height = 18, ncol = 2)
}
