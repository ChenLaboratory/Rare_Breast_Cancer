## 07_Pseudobulk_Tumour.R

source("setup.R")
source("config.R")

subtypes <- subtypes_tumour

## ---- Pseudo-bulk profiles of the tumour populations ----
dge_list <- lapply(setNames(subtypes, subtypes), function(s) {
  so <- readRDS(file.path(RDS_DIR, paste0("so_", s, "_clean.RDS")))
  y <- seurat2pb(so, sample = "orig.ident", cluster = "cell_type_clean")
  colnames(y) <- gsub("_cluster", "_", colnames(y))
  y$samples$subtype <- s
  y
})
saveRDS(dge_list, file.path(GENELIST_DIR, "dge_list_Tumour.RDS"))

y <- do.call(cbind, dge_list)
y <- y[, grepl("^Tumor$|^Tumor_", y$samples$cluster)]

## Two populations are dropped after manual review (see config.R).
drop <- Reduce(`|`, lapply(pb_tumour_exclude, function(ex) {
  y$samples$cluster == ex$cluster & y$samples$subtype == ex$subtype
}))
print(table(keep = !drop))
y <- y[, !drop]

y <- sumTechReps(y, ID = factor(y$samples$sample))
y <- add_entrez(y)

keep_sample <- y$samples$lib.size > pb_tumour_lib_cutoff
print(table(keep_sample))
y <- y[, keep_sample]

lcpm_all <- edgeR::cpm(normLibSizes(y), log = TRUE)
rownames(lcpm_all) <- y$genes$Symbol
write.csv(lcpm_all, file.path(GENELIST_DIR, "logCPM_Tumour.csv"))

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

save_figure({
    mds <- plotMDS(y, plot = FALSE)
    plotMDS(mds, col = subtype_cols[group], main = "MDS by subtype")
    legend("bottom", legend = levels(group), col = subtype_cols, pch = 16)
  }, "Tumour_MDS", width = 8, height = 8)

## ---- Quasi-likelihood GLM fit ----
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
fit <- glmQLFit(y, design, robust = TRUE)
save_figure(plotQLDisp(fit), "Tumour_QLDisp", width = 6, height = 6)

## ---- Helper: run one set of contrasts ----
run_contrast_set <- function(contr, prefix = "Tumour") {
  qlf_list <- lapply(seq_len(ncol(contr)), function(i) {
    glmQLFTest(fit, contrast = contr[, i])
  })
  names(qlf_list) <- colnames(contr)

  for (nm in names(qlf_list)) {
    qlf <- qlf_list[[nm]]
    keg <- kegga(qlf, geneid = y$genes$GeneID, species = "Hs")

    cat("\n\n#### ", nm, "\n")
    print(topTags(qlf))
    print(summary(decideTests(qlf)))
    print(topKEGG(keg, sort = "up", truncate = 45))
    print(topKEGG(keg, sort = "down", truncate = 45))

    write.csv(topTags(qlf, n = Inf),
              file.path(GENELIST_DIR, paste0("DEG_", prefix, "_", nm, ".csv")))
    write.csv(topKEGG(keg, n = Inf),
              file.path(GENELIST_DIR, paste0("KEGG_", prefix, "_", nm, ".csv")))

    save_figure(plotMD(qlf), paste0(prefix, "_MD_", nm), width = 6, height = 6)
    save_figure(plot_volcano(qlf, title = nm),
                paste0(prefix, "_Volcano_", nm), width = 7, height = 7)
  }
  qlf_list
}

## ---- Pairwise comparisons ----
contr <- makeContrasts(
  TNBC_vs_MED = TNBC - MED,
  TNBC_vs_APO = TNBC - APO,
  TNBC_vs_CYS = TNBC - CYS,
  ER_vs_IMPC  = ER - IMPC,
  ER_vs_PLC   = ER - PLC, levels = design)
qlf_pairwise <- run_contrast_set(contr)

## ---- Each subtype against the others of its group ----
group_contrasts <- list(
  ER = list(
    subtypes = c("ER", "IMPC", "PLC"),
    contrasts = list(ER_vs_IMPC_PLC = "ER   - (IMPC + PLC)/2",
                     IMPC_vs_ER_PLC = "IMPC - (ER   + PLC)/2",
                     PLC_vs_ER_IMPC = "PLC  - (ER   + IMPC)/2")),
  TNBC = list(
    subtypes = c("TNBC", "CYS", "MED", "APO"),
    contrasts = list(TNBC_vs_CYS_MED_APO = "TNBC - (CYS + MED + APO)/3",
                     CYS_vs_TNBC_MED_APO = "CYS  - (TNBC + MED + APO)/3",
                     MED_vs_TNBC_CYS_APO = "MED  - (TNBC + CYS + APO)/3",
                     APO_vs_TNBC_CYS_MED = "APO  - (TNBC + CYS + MED)/3"))
)

for (grp in names(group_contrasts)) {
  spec <- group_contrasts[[grp]]
  contr_grp <- do.call(makeContrasts, c(spec$contrasts, list(levels = design)))
  colnames(contr_grp) <- names(spec$contrasts)
  qlf_grp <- run_contrast_set(contr_grp)

  for (n_top in heatmap_top_n) {
    save_figure(
      de_heatmap(qlf_grp, lcpm, group, spec$subtypes, subtype_cols,
                 top = n_top, exclusive = TRUE,
                 title = paste0("Top ", n_top, " exclusive up-regulated DE genes: ",
                                paste(spec$subtypes, collapse = " / "))),
      paste0("Tumour_Heatmap_", grp, "_top", n_top), width = 7, height = 14)
  }
}
