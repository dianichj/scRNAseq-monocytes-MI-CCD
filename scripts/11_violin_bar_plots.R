# 11_violin_bar_plots.R
# This script creates violin plots and bar plots for quality metrics and cell-type distributions

library(Seurat)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggsignif)

# Load the object
mono_object <- readRDS("data/Mono_object.rds")

# Set identities by annotation
Idents(mono_object) <- mono_object$annotated

# Define genes of interest for violin plots
genes_of_interest <- c("CD14", "S100A12", "CSF3R", "VCAN", "MS4A6A", "FCGR3A", "CDKN1C", "HES4", "IFI44", "IFI44L", "IFI6", "CLEC4C")

# Violin Plot for genes
violin_genes <- VlnPlot(mono_object, features = genes_of_interest) + theme(legend.position = "none")
ggsave("plots/violin_gene_expression.png", violin_genes, width = 14, height = 14, dpi = 600)

# Quality metrics violin plots
v1 <- VlnPlot(mono_object, features = "nFeature_RNA") + ggtitle("nFeature_RNA by Cluster")
v2 <- VlnPlot(mono_object, features = "nCount_RNA", group.by = "condition", pt.size = 0) + ggtitle("nCount_RNA by Condition")
v3 <- VlnPlot(mono_object, features = "nFeature_RNA", group.by = "batch", pt.size = 0) + ggtitle("nFeature_RNA by Sample")
v4 <- VlnPlot(mono_object, features = "nCount_RNA", group.by = "batch", pt.size = 0) + ggtitle("nCount_RNA by Sample")

for (i in list(v1, v2, v3, v4)) {
  ggsave(paste0("plots/", gsub(" ", "_", i$labels$title), ".png"), plot = i, width = 10, height = 8, dpi = 300)
}

# Bar plot: percentage of annotated cells
cell_type_percentages <- mono_object@meta.data %>%
  count(annotated) %>%
  mutate(percentage = n / sum(n) * 100)

bar_plot1 <- ggplot(cell_type_percentages, aes(x = annotated, y = percentage, fill = annotated)) +
  geom_bar(stat = "identity") +
  labs(x = "Cell Type", y = "Percentage (%)", title = "Percentage of Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

ggsave("plots/barplot_celltype_percentage.png", plot = bar_plot1, width = 10, height = 6, dpi = 300)

# Bar plot: condition vs. cell type
cell_type_condition <- mono_object@meta.data %>%
  group_by(condition, annotated) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100)

bar_plot2 <- ggplot(cell_type_condition, aes(x = condition, y = Percentage, fill = annotated)) +
  geom_bar(stat = "identity", width = 0.3) +
  labs(x = "Condition", y = "Percentage (%)", title = "Cell Type Proportions per Condition") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Paired")

ggsave("plots/barplot_celltype_by_condition.png", plot = bar_plot2, width = 10, height = 6, dpi = 300)
