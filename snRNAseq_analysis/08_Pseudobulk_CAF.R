## 08_Pseudobulk_CAF.R

source("setup.R")
source("config.R")

subtypes <- subtypes_caf

## ---- Pseudo-bulk profiles of the CAF clusters ----
so_list <- lapply(setNames(subtypes, subtypes), function(s) {
  readRDS(file.path(RDS_DIR, paste0("so_", s, "_stroma_fibroblast.RDS")))
})

save_figure(lapply(subtypes, function(s) {
    DimPlot(so_list[[s]], pt.size = 0.2, raster = FALSE, cols = col.p2) + ggtitle(s)
  }), "CAF_UMAP_bySubtype", width = 14, height = 13, ncol = 3)

y_list <- lapply(so_list, function(so) {
  y <- seurat2pb(so, sample = "orig.ident", cluster = "seurat_clusters")
  colnames(y) <- gsub("cluster", "C", colnames(y))
  y$genes <- y$genes[, 1, drop = FALSE]
  y
})
for (s in subtypes) y_list[[s]]$samples$subtype <- s

y_all <- do.call(cbind, y_list)
y <- sumTechReps(y_all, ID = y_all$samples$sample)
y <- add_entrez(y)

lcpm_all <- edgeR::cpm(normLibSizes(y), log = TRUE)
rownames(lcpm_all) <- y$genes$Symbol
write.csv(lcpm_all, file.path(GENELIST_DIR, "logCPM_CAF.csv"))

## ---- Filtering and normalisation ----
group <- y$samples$group <- factor(y$samples$subtype, levels = subtypes)
keep <- filterByExpr(y)
print(table(keep))
y <- y[keep, , keep = FALSE]
y <- normLibSizes(y)
print(y$samples)

lcpm <- edgeR::cpm(y, log = TRUE)
rownames(lcpm) <- y$genes$Symbol
subtype_cols <- setNames(col.p2[seq_along(levels(group))], levels(group))

mds_sets <- list(all = subtypes,
                 TNBC = c("TNBC", "MED", "CYS", "APO"),
                 ER = c("ER", "IMPC", "PLC"))
for (nm in names(mds_sets)) {
  set <- mds_sets[[nm]]
  save_figure({
      keep_cols <- group %in% set
      mds <- plotMDS(y[, keep_cols], plot = FALSE)
      plotMDS(mds, col = subtype_cols[group[keep_cols]],
              main = paste("MDS by subtype:", paste(set, collapse = " / ")))
      legend("bottomleft", legend = set, col = subtype_cols[set], pch = 16)
    }, paste0("CAF_MDS_", nm), width = 9, height = 9)
}

## ---- Quasi-likelihood GLM fit ----
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
fit <- glmQLFit(y, design, robust = TRUE)
save_figure(plotQLDisp(fit), "CAF_QLDisp", width = 6, height = 6)

## ---- Pairwise comparisons ----
contr <- makeContrasts(
  TNBC_vs_MED = TNBC - MED,
  TNBC_vs_CYS = TNBC - CYS,
  TNBC_vs_APO = TNBC - APO,
  ER_vs_IMPC  = ER - IMPC,
  ER_vs_PLC   = ER - PLC,
  IMPC_vs_PLC = IMPC - PLC, levels = design)

qlf_list <- lapply(seq_len(ncol(contr)), function(i) glmQLFTest(fit, contrast = contr[, i]))
names(qlf_list) <- colnames(contr)

for (nm in names(qlf_list)) {
  qlf <- qlf_list[[nm]]
  keg <- kegga(qlf, geneid = y$genes$GeneID, species = "Hs", FDR = 0.2)

  write.csv(topTags(qlf, n = Inf),
            file.path(GENELIST_DIR, paste0("DEG_CAF_", nm, ".csv")))
  write.csv(topKEGG(keg, n = Inf),
            file.path(GENELIST_DIR, paste0("KEGG_CAF_", nm, ".csv")))
  save_figure(plotMD(qlf), paste0("CAF_MD_", nm), width = 7, height = 7)
}

heatmap_sets <- list(
  all  = list(subtypes = subtypes,
              contrasts = names(qlf_list)),
  TNBC = list(subtypes = c("TNBC", "MED", "CYS", "APO"),
              contrasts = c("TNBC_vs_MED", "TNBC_vs_CYS", "TNBC_vs_APO")),
  ER   = list(subtypes = c("ER", "IMPC", "PLC"),
              contrasts = c("ER_vs_IMPC", "ER_vs_PLC", "IMPC_vs_PLC"))
)
for (nm in names(heatmap_sets)) {
  spec <- heatmap_sets[[nm]]
  save_figure(
    de_heatmap(qlf_list[spec$contrasts], lcpm, group, spec$subtypes, subtype_cols,
               top = 20, exclusive = FALSE, up_only = FALSE,
               title = paste("Top DE genes:", paste(spec$subtypes, collapse = " / "))),
    paste0("CAF_Heatmap_pairwise_", nm), width = 9, height = 12)
}

## ---- Each subtype against all others ----
vs_others <- lapply(setNames(subtypes, subtypes), function(s) {
  others <- subtypes[subtypes != s]
  formula <- paste0(s, " - (", paste(others, collapse = " + "), ")/", length(others))
  contr_vs <- makeContrasts(contrasts = formula, levels = design)
  colnames(contr_vs) <- paste0(s, "_vs_Others")

  lfc <- if (s %in% names(caf_treat_lfc)) caf_treat_lfc[[s]] else log2(1.3)
  qlf <- glmQLFTest(fit, contrast = contr_vs)
  trt <- glmTreat(fit, contrast = contr_vs, lfc = lfc)
  keg <- kegga(qlf, geneid = y$genes$GeneID, species = "Hs", FDR = 0.2)

  write.csv(topTags(qlf, n = Inf),
            file.path(GENELIST_DIR, paste0("DEG_CAF_", s, "_vs_Others.csv")))
  write.csv(topKEGG(keg, n = Inf),
            file.path(GENELIST_DIR, paste0("KEGG_CAF_", s, "_vs_Others.csv")))
  write.csv(topTags(trt, n = Inf),
            file.path(GENELIST_DIR, paste0("DEG_CAF_", s, "_vs_Others_TREAT.csv")))

  list(qlf = qlf, trt = trt, keg = keg)
})

for (s in subtypes) {
  res <- vs_others[[s]]
  save_figure(plotMD(res$qlf), paste0("CAF_MD_", s, "_vs_Others_QLF"),
              width = 7, height = 7)
  save_figure(plotMD(res$trt), paste0("CAF_MD_", s, "_vs_Others_TREAT"),
              width = 7, height = 7)
}

for (n_top in heatmap_top_n) {
  for (test in c("qlf", "trt")) {
    for (excl in c(FALSE, TRUE)) {
      if (test == "trt" && !excl) next  ## the TREAT results are only shown exclusive
      qlf_set <- lapply(vs_others, function(r) r[[test]])
      label <- paste0(if (test == "trt") "TREAT " else "",
                      "Top ", n_top, if (excl) " exclusive " else " ",
                      "up-regulated DE genes: each subtype vs others")
      save_figure(
        de_heatmap(qlf_set, lcpm, group, subtypes, subtype_cols,
                   top = n_top, exclusive = excl, title = label),
        paste0("CAF_Heatmap_vsOthers_", test, if (excl) "_exclusive" else "",
               "_top", n_top),
        width = 9, height = if (n_top == 40) 22 else 18)
    }
  }
}

## ---- CAF signature scores ----
sig_df <- read.delim(CAF_SIGNATURE_FILE, stringsAsFactors = FALSE)
sig_list <- split(sig_df$Gene, factor(sig_df$CellType, levels = unique(sig_df$CellType)))

score_mat <- sapply(sig_list, function(genes) {
  colMeans(lcpm_all[intersect(genes, rownames(lcpm_all)), , drop = FALSE])
})

long <- do.call(rbind, lapply(colnames(score_mat), function(sig) {
  data.frame(sample = rownames(score_mat),
             subtype = factor(group, levels = subtypes),
             Signature = sig,
             Score = score_mat[, sig])
}))
long$Signature <- factor(long$Signature, levels = colnames(score_mat))

save_figure(
  ggplot(long, aes(x = subtype, y = Score, fill = subtype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    facet_wrap(~ Signature, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = subtype_cols) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.text = element_text(face = "bold")) +
    labs(x = "Subtype", y = "Signature score (mean logCPM)",
         title = "CAF signature scores across cancer subtypes"),
  "CAF_Signature_scores", width = 12, height = 8)
