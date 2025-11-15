# ===============================
# 0) Setup & project parameters
# ===============================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
})

# ---- paths (edit base_dir only) ----
base_dir <- normalizePath("C:/Users/mmsid/Documents/My Project/scRNA Analysis/GSE111429_RAW", mustWork = TRUE)
raw_dir  <- base_dir                     # your barcodes.tsv, genes.tsv, matrix.mtx live here
out_dir  <- file.path(base_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE)
dirs <- list(
  qc   = file.path(out_dir, "qc"),
  pca  = file.path(out_dir, "pca"),
  emb  = file.path(out_dir, "embedding"),
  mrk  = file.path(out_dir, "markers"),
  obj  = file.path(out_dir, "objects")
)
lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# ---- analysis knobs ----
min_genes        <- 1000       # filter: ≥ this many detected genes
regress_percent_mt <- TRUE     # set FALSE to skip slow regression
pcs_for_embed    <- 1:12       # match paper-style (can change based on elbow)
umap_or_tsne     <- "tsne"     # "umap" or "tsne"
res_clustering   <- 0.25       # adjust to merge/split clusters
k_param_neighbors <- 30

# ===============================
# 1) Load matrix & basic QC
# ===============================
mtx <- Read10X(data.dir = raw_dir)
obj <- CreateSeuratObject(mtx, project = "GSE111429")

mito_pat <- if (any(grepl("^MT-", rownames(obj)))) "^MT-" else "^mt-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pat)

# Save raw QC violins (before filtering)
p_before <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)
ggsave(file.path(dirs$qc, "QC_violin_before.png"), p_before, width = 9, height = 3.5, dpi = 300)

# Filter
obj_filt <- subset(obj, subset = nFeature_RNA >= min_genes)

# Save after QC violins
p_after <- VlnPlot(obj_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)
ggsave(file.path(dirs$qc, "QC_violin_after.png"), p_after, width = 9, height = 3.5, dpi = 300)

# Combined before/after panel
obj$QC_stage <- "Before"; obj_filt$QC_stage <- "After"
mm <- merge(obj, y = obj_filt)
mm$nGene <- mm$nFeature_RNA; mm$nUMI <- mm$nCount_RNA; mm$percent.mito <- mm$percent.mt/100
df <- FetchData(mm, vars = c("QC_stage","nGene","nUMI","percent.mito")) %>%
  mutate(QC_stage = factor(QC_stage, levels = c("Before","After"))) %>%
  pivot_longer(cols = c(nGene, nUMI, percent.mito), names_to = "metric", values_to = "value")
p_combo <- ggplot(df, aes(QC_stage, value)) +
  geom_violin(fill = NA, color = "black", linewidth = 0.6, scale = "width") +
  geom_jitter(width = 0.25, size = 0.2, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", color = "red", size = 1.2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(x = "Identity", y = NULL) + theme_bw() +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(dirs$qc, "QC_violin_before_after.png"), p_combo, width = 9, height = 3.2, dpi = 300)

saveRDS(obj_filt, file.path(dirs$obj, "step1_qc_filtered.rds"))

# ===============================
# 2) Normalize, HVGs, Scale, PCA
# ===============================
obj_filt <- NormalizeData(obj_filt, normalization.method = "LogNormalize", scale.factor = 1e4)
obj_filt <- FindVariableFeatures(obj_filt, selection.method = "vst", nfeatures = 2000)
hvg <- VariableFeatures(obj_filt)

# HVG scatter
p_hvg <- LabelPoints(plot = VariableFeaturePlot(obj_filt), points = head(hvg, 15), repel = TRUE)
ggsave(file.path(dirs$pca, "HVG_scatter.png"), p_hvg + theme_bw(), width = 7, height = 5, dpi = 300,  bg = "white")

# Scale (optionally regress percent.mt)
obj_filt <- ScaleData(
  obj_filt,
  features = rownames(obj_filt),
  vars.to.regress = if (regress_percent_mt) "percent.mt" else NULL
)

# PCA
obj_filt <- RunPCA(obj_filt, features = hvg, npcs = 50, verbose = FALSE)

# PCA diagnostics
png(file.path(dirs$pca, "PCA_loadings_PC1-3.png"), width = 2000, height = 700, res = 200)
print(VizDimLoadings(obj_filt, dims = 1:3, reduction = "pca")); dev.off()

p_elbow <- ElbowPlot(obj_filt, ndims = 50)
ggsave(file.path(dirs$pca, "ElbowPlot.png"), p_elbow, width = 6, height = 4, dpi = 300)

# PC heatmaps 1–12 (combined + individual)
png(file.path(dirs$pca, "PC_heatmaps_1-12.png"), width = 2400, height = 1800, res = 200)
print(DimHeatmap(obj_filt, dims = 1:12, cells = 500, balanced = TRUE, fast = TRUE)); dev.off()
for (k in 1:12) {
  png(file.path(dirs$pca, sprintf("PC_%02d_heatmap.png", k)), width = 1200, height = 900, res = 200)
  print(DimHeatmap(obj_filt, dims = k, cells = 500, balanced = TRUE, fast = TRUE)); dev.off()
}

saveRDS(obj_filt, file.path(dirs$obj, "step2_afterPCA.rds"))

# ===============================
# 3) Neighbors, clustering, embedding
# ===============================
obj_filt <- FindNeighbors(obj_filt, dims = pcs_for_embed, k.param = k_param_neighbors)
obj_filt <- FindClusters(obj_filt, resolution = res_clustering)

if (umap_or_tsne == "umap") {
  obj_filt <- RunUMAP(obj_filt, dims = pcs_for_embed)
  p_emb <- DimPlot(obj_filt, reduction = "umap", label = TRUE) + ggtitle(sprintf("UMAP (PCs %s)", paste(range(pcs_for_embed), collapse = "–")))
  ggsave(file.path(dirs$emb, "UMAP.png"), p_emb, width = 7, height = 6, dpi = 300)
} else {
  obj_filt <- RunTSNE(obj_filt, dims = pcs_for_embed, perplexity = 30)
  p_emb <- DimPlot(obj_filt, reduction = "tsne", label = TRUE) + ggtitle(sprintf("tSNE (PCs %s)", paste(range(pcs_for_embed), collapse = "–")))
  ggsave(file.path(dirs$emb, "tSNE.png"), p_emb, width = 7, height = 6, dpi = 300)
}

# marker genes overlay (edit list as needed)
marker_panel <- c("Krt5","Krt14","Trp63","Krt8","Krt18","Epcam")
p_feat <- FeaturePlot(obj_filt, features = marker_panel, ncol = 3, reduction = ifelse(umap_or_tsne=="umap","umap","tsne"))
ggsave(file.path(dirs$emb, "FeaturePlot_marker_panel.png"), p_feat, width = 10, height = 7, dpi = 300)

# QC metrics by cluster
p_qc_by_clust <- VlnPlot(obj_filt, c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
ggsave(file.path(dirs$emb, "QC_violin_by_cluster.png"), p_qc_by_clust, width = 12, height = 4, dpi = 300)

saveRDS(obj_filt, file.path(dirs$obj, "step3_clustered.rds"))

# ===============================
# 4) Differential markers
# ===============================
DefaultAssay(obj_filt) <- "RNA"
markers <- FindAllMarkers(
  obj_filt,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)
write.csv(markers, file.path(dirs$mrk, "markers_per_cluster.csv"), row.names = FALSE)

# top10 heatmap
top10 <- markers %>% group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
p_heat <- DoHeatmap(obj_filt, features = unique(top10$gene), group.by = "seurat_clusters") +
  ggtitle("Top10 markers per cluster")
ggsave(file.path(dirs$mrk, "markers_top10_heatmap.png"), p_heat, width = 10, height = 7, dpi = 300)

saveRDS(obj_filt, file.path(dirs$obj, "step4_markers_done.rds"))

# ===============================
# 5) Session info (reproducibility)
# ===============================
writeLines(c(capture.output(sessionInfo()), ""), con = file.path(out_dir, "sessionInfo.txt"))


## =========================================================
## FeaturePlot with cluster labels overlaid (tSNE embedding)
## =========================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# 1) Choose the genes you want to show
marker_genes <- c("Krt5","Krt14","Trp63","Krt8","Krt18","Epcam")

# 2) Get tSNE coordinates and cluster IDs, then compute cluster centers
coords <- Embeddings(obj_filt, "tsne") %>%
  as.data.frame() %>%
  mutate(cluster = Idents(obj_filt))

centers <- coords %>%
  group_by(cluster) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2),
    .groups = "drop"
  )

# 3) Make the basic FeaturePlot panel (one subplot per gene)
p_panel <- FeaturePlot(
  obj_filt,
  features  = marker_genes,
  reduction = "tsne",
  ncol      = 3
)

# 4) Add cluster labels (0,1,2,...) to every subplot
for (i in seq_along(p_panel)) {
  p_panel[[i]] <- p_panel[[i]] +
    geom_text(
      data = centers,
      aes(x = tSNE_1, y = tSNE_2, label = cluster),
      color    = "black",
      size     = 4,
      fontface = "bold"
    )
}

# 5) Show in RStudio
p_panel

# 6) Save to file (adjust outdir if needed)
outdir <- "C:/Users/mmsid/Documents/My Project/scRNA Analysis/GSE111429_RAW"
ggsave(
  filename = file.path(outdir, "FeaturePlot_marker_panel_withClusters.png"),
  plot     = p_panel,
  width    = 10,
  height   = 7,
  dpi      = 300
)



markers <- FindAllMarkers(
  obj_filt,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(markers)

markers %>%
  filter(grepl("^Krt", gene)) %>%
  arrange(cluster, desc(avg_log2FC))

