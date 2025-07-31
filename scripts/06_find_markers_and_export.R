# 06_find_markers_and_export.R

# Description:
# This script identifies marker genes for each Seurat cluster using the SCT assay
# and exports both the full marker table and top 10 genes per cluster.

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load object
mono_object <- readRDS("data/Mono_object.rds")

# Set identity and assay
DefaultAssay(mono_object) <- "SCT"
Idents(mono_object) <- "seurat_clusters"

# Prepare for marker discovery
mono_object <- PrepSCTFindMarkers(mono_object)

# Find markers
markers <- FindAllMarkers(
  mono_object,
  min.pct = 0.2,
  logfc.threshold = 0.25,
  only.pos = TRUE,
  assay = "SCT"
)

# Save full marker table
dir.create("tables", showWarnings = FALSE)
write.csv(markers, "tables/Mono_markers_seurat.csv", row.names = FALSE)

# Extract top 10 per cluster for heatmap
top10 <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

# Plot heatmap
p1 <- DoHeatmap(subset(mono_object, downsample = 50), features = top10$gene, size = 1) + NoLegend()
p1 + theme(axis.text.y = element_text(size = 20))

# Save heatmap
ggsave("plots/batch_effects/heatmap_seurat_markers.jpg", plot = p1, units = "cm", height = 80, width = 60)

# Save object (in case needed downstream)
saveRDS(mono_object, "data/Mono_object.rds")