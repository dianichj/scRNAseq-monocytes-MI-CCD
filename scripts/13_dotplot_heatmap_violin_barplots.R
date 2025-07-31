# 13_dotplot_heatmap_violin_barplots.R

# Description:
# This script generates dot plots, heatmaps, violin plots, and bar plots for gene expression analysis across annotated monocyte subtypes and conditions.
# Subtypes include: Classical Monocytes, Non-Classical Monocytes, IFN Classical Monocytes
# Conditions include: MI and CCD (myocardial infarction and chronic coronary disease)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggsignif)
library(data.table)

# Load the Seurat object
mono_object <- readRDS("data/Mono_object.rds")

# Set identities
Idents(mono_object) <- mono_object$annotated

# --- Dot Plot ---
genes_of_interest <- c("CD14", "S100A12", "VCAN", "CD36", "S100A8", "S100A9", "FCGR3A", "CDKN1C", "BCL2A1", "CTSL", "IFI44L", "IFI44", "IFI6", "IFIT3", "CLEC4C", "C12orf75", "P2RY14", "NIBAN3", "COBLL1")

dot_plot <- DotPlot(mono_object, features = genes_of_interest, group.by = "annotated", dot.scale = 8) +
  scale_color_gradientn(colors = c("blue", "#FF0000")) +
  xlab("Features") + ylab("Cell Type") + ggtitle("Gene Expression by Monocyte Subtype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16))

ggsave("plots/plots_expression/dot_plot_genes.png", plot = dot_plot, width = 14, height = 8, dpi = 600)

# --- Violin Plot ---
genes_violin <- c("CD14", "S100A12", "CSF3R", "VCAN", "MS4A6A", "FCGR3A", "CDKN1C", "HES4", "IFI44", "IFI44L", "IFI6", "CLEC4C")

violin_plot <- VlnPlot(mono_object, features = genes_violin) +
  theme(legend.position = "none")

ggsave("plots/plots_expression/violin_plot_genes.png", plot = violin_plot, width = 14, height = 14, dpi = 600)

# --- Bar Plot: Percentage of Annotated Cells ---
cell_data <- mono_object@meta.data

cell_type_percentages <- cell_data %>%
  count(annotated) %>%
  mutate(percentage = n / sum(n) * 100)

bar_plot <- ggplot(cell_type_percentages, aes(x = annotated, y = percentage, fill = annotated)) +
  geom_bar(stat = "identity") +
  labs(x = "Cell Type", y = "Percentage (%)", title = "Percentage of Annotated Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

ggsave("plots/plots_expression/barplot_celltypes.png", plot = bar_plot, width = 10, height = 6, dpi = 300)

# --- Bar Plot: Cell Types by Condition ---
cell_type_percentages_condition <- cell_data %>%
  group_by(condition, annotated) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(condition) %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100)

bar_condition <- ggplot(cell_type_percentages_condition, aes(x = condition, y = Percentage, fill = annotated)) +
  geom_bar(stat = "identity", width = 0.3) +
  labs(x = "Condition", y = "Percentage (%)", title = "Percentage of Cell Types per Condition") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "gray98"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(palette = "Paired")

ggsave("plots/plots_expression/barplot_celltypes_per_condition.png", plot = bar_condition, width = 10, height = 6, dpi = 300)
