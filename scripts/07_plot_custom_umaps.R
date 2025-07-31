07_plot_custom_umaps.R

# Description:
# Generates custom UMAP plots for annotated clusters and Seurat clusters,
# placing bold labels at centroid positions.

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load object
mono_object <- readRDS("data/Mono_object.rds")

# Ensure UMAPs are computed for annotated view
mono_object <- RunUMAP(mono_object, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap_cca", a = 1, b = 0.7)

# Plot by cell type annotations
umap_data <- FetchData(mono_object, vars = c("umapcca_1", "umapcca_2", "annotated"))
centroids <- umap_data %>% group_by(annotated) %>% summarize(across(starts_with("umap"), mean))

umap_plot_annotated <- DimPlot(mono_object, group.by = "annotated", reduction = "umap_cca", label = FALSE, repel = TRUE) +
  geom_text(data = centroids, aes(x = umapcca_1, y = umapcca_2, label = annotated), size = 8, fontface = "bold") +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 23, face = "bold", hjust = 0.05),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    legend.text = element_text(size = 23),
    legend.title = element_text(size = 23, face = "bold")
  )

ggsave("plots/umaps/umap_annotated_high_res.png", plot = umap_plot_annotated, height = 18, width = 18, dpi = 600)

# Plot by Seurat clusters
umap_data2 <- FetchData(mono_object, vars = c("umapcca_1", "umapcca_2", "seurat_clusters"))
centroids2 <- umap_data2 %>% group_by(seurat_clusters) %>% summarize(across(starts_with("umap"), mean))

umap_plot_clusters <- DimPlot(mono_object, group.by = "seurat_clusters", reduction = "umap_cca", label = FALSE, repel = TRUE) +
  geom_text(data = centroids2, aes(x = umapcca_1, y = umapcca_2, label = seurat_clusters), size = 8, fontface = "bold") +
  ggtitle("Seurat Clusters") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  theme(legend.text = element_text(face = "bold", size = 10))

# Save cluster plot
ggsave("plots/umaps/umap_seurat_clusters_high_res.png", plot = umap_plot_clusters, height = 18, width = 18, dpi = 600)
