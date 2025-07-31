# 04_dimensionality_reduction_and_umaps.R

# Description:
# This script runs UMAPs on each of the integrated reductions (Harmony, CCA, scVI)
# and exports the corresponding plots.

# Load libraries
library(Seurat)
library(ggplot2)

# Load integrated object
mono_object <- readRDS("data/Mono_object.rds")

# Ensure plots directory exists
dir.create("plots/batch_effects", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# UMAP - Harmony
# -----------------------------
mono_object <- RunUMAP(
  mono_object,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  a = 1, b = 1
)

DimPlot(mono_object, reduction = "umap_harmony", group.by = "tissue")
ggsave("plots/batch_effects/umap_harmony.jpg", units = "cm", height = 18, width = 18)

DimPlot(mono_object, reduction = "umap_harmony", split.by = "batch", ncol = 5) & NoLegend()
ggsave("plots/batch_effects/umap_harmony_split_by_batch.jpg", units = "cm", height = 16, width = 28)

# -----------------------------
# UMAP - CCA
# -----------------------------
mono_object <- RunUMAP(
  mono_object,
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap_cca",
  a = 1, b = 0.8
)

DimPlot(mono_object, reduction = "umap_cca", group.by = "tissue")
ggsave("plots/batch_effects/umap_cca_tissue.jpg", units = "cm", height = 18, width = 18, dpi = 600)

# Custom color scheme for condition
custom_colors <- c("MI" = "#FF91A4", "CCD" = "#4169E1")
DimPlot(mono_object, reduction = "umap_cca", group.by = "condition") + 
  scale_color_manual(values = custom_colors)
ggsave("plots/batch_effects/umap_cca_condition.jpg", units = "cm", height = 18, width = 18, dpi = 600)

DimPlot(mono_object, reduction = "umap_cca", split.by = "batch", ncol = 5) & NoLegend()
ggsave("plots/batch_effects/umap_cca_split_by_batch.jpg", units = "cm", height = 16, width = 28)

# -----------------------------
# UMAP - scVI
# -----------------------------
mono_object <- RunUMAP(
  mono_object,
  reduction = "integrated.scvi",
  dims = 1:30,
  reduction.name = "umap_scvi",
  a = 1, b = 1
)

DimPlot(mono_object, reduction = "umap_scvi", group.by = "tissue")
ggsave("plots/batch_effects/umap_scvi.jpg", units = "cm", height = 18, width = 18)

DimPlot(mono_object, reduction = "umap_scvi", split.by = "batch", ncol = 5) & NoLegend()
ggsave("plots/batch_effects/umap_scvi_split_by_batch.jpg", units = "cm", height = 16, width = 28)

# Save updated object
saveRDS(mono_object, "data/Mono_object.rds")