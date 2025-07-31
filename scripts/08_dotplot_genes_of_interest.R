# 08_dotplot_genes_of_interest.R

# Description:
# This script generates a DotPlot of selected genes across annotated monocyte subtypes.

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load Seurat object
mono_object <- readRDS("data/Mono_object.rds")

# Ensure consistent annotation labels
mono_object$annotated <- recode(mono_object$annotated,
                                "Classical_Monocyte" = "Classical Monocyte",
                                "Non_Classical_Monocyte" = "Non Classical Monocyte",
                                "IFN" = "IFN Classical Monocyte")

# Define genes of interest
genes_of_interest <- unique(c("CD14", "S100A12", "VCAN", "CD36", "S100A8", "S100A9", 
                              "MS4A6A", "FCGR3A", "CDKN1C", "HES4", "IFI44L", "IFI44", 
                              "IFI6", "IFIT3", "MX1", "EPSTI1", "CLEC4C", "LILRA4", "P2RY14"))

# Generate DotPlot
plot_dot <- DotPlot(mono_object, features = genes_of_interest, group.by = "annotated", dot.scale = 8) +
  xlab("Genes") + ylab("Cell Type") + ggtitle("Gene Expression by Cluster") +
  scale_color_gradientn(colors = c("blue", "#FF0000")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# Display and save plot
print(plot_dot)
dir.create("plots/", showWarnings = FALSE)
ggsave("plots/dot_plot_gene_expression.png", plot = plot_dot, width = 14, height = 8, dpi = 600)
