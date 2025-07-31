# 10_heatmaps_and_featureplots.R

# Description:
# Generates heatmaps and feature plots for selected marker genes across annotated clusters.

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load Seurat object
mono_object <- readRDS("data/Mono_object.rds")

# Directory for plots
dir.create("plots/heatmaps", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/features", recursive = TRUE, showWarnings = FALSE)

# Top markers per cluster (if not precomputed)
if (!"markers" %in% names(mono_object@misc)) {
  mono_object <- FindAllMarkers(mono_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top_genes <- mono_object %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) %>%
    pull(gene)
} else {
  top_genes <- mono_object@misc$markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) %>%
    pull(gene)
}

# Heatmap of top markers
heatmap <- DoHeatmap(mono_object, features = unique(top_genes)) + NoLegend()
ggsave("plots/heatmaps/top_markers_heatmap.png", plot = heatmap, width = 12, height = 10, dpi = 300)

# Feature plots for selected genes
genes_to_plot <- c("CD14", "FCGR3A", "CLEC4C", "IFI44", "IFI6")

for (gene in genes_to_plot) {
  p <- FeaturePlot(mono_object, reduction = "umap_cca", features = gene, cols = c("#D3D3D3", "red")) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()

  ggsave(paste0("plots/features/featureplot_", gene, ".png"), plot = p, width = 8, height = 8, dpi = 300)
}