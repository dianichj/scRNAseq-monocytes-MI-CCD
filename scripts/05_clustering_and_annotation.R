# 05_clustering_and_annotation.R

# Description:
# Perform clustering and cluster annotation using Seurat. This script
# uses the integrated CCA reduction for clustering and exports the UMAP with cluster labels.

# Load libraries
library(Seurat)
library(ggplot2)

# Load object
mono_object <- readRDS("data/Mono_object.rds")

# Set default assay to SCT
DefaultAssay(mono_object) <- "SCT"

# Find neighbors and clusters using CCA reduction
mono_object <- FindNeighbors(mono_object, dims = 1:30, reduction = "integrated.cca")
mono_object <- FindClusters(mono_object, resolution = 0.8, algorithm = 1, graph.name = "SCT_snn")

# Plot clusters
DimPlot(mono_object, reduction = "umap_cca", label = TRUE, repel = TRUE, group.by = "seurat_clusters")
ggsave("plots/batch_effects/umap_clusters.jpg", height = 18, width = 18, units = "cm")

# Save clustered object
saveRDS(mono_object, "data/Mono_object.rds")