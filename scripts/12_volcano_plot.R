# 12_volcano_plot.R
# Generates volcano plots comparing blood vs. bone marrow for all monocyte subtypes

library(ggplot2)
library(ggrepel)
library(dplyr)

# Configuration
monocyte_types <- c("CM", "NCM", "IFNCM")
input_folder <- "data/"
output_folder <- "plots/volcano/"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Immune-related genes to highlight
immune_genes <- c("CD14", "CD16", "VCAN", "LYZ", "FCGR3A", "FCN1", "CDKN1C", "S100A8", "S100A9", 
                  "CSF3R", "MS4A6A", "CX3CR1", "ITGAM", "ITGAX", "TREM2", "MERTK", "CD163", "CCL3", 
                  "CCL4", "HLA-DRA", "HLA-DRB1", "ICAM1", "TNF", "IFNG", "IL6", "CXCL10", "CCL5", 
                  "BACH2", "MKI67")

# Parameters
logFC_threshold <- 0
FDR_threshold <- 0.05

# Loop through monocyte types
for (type in monocyte_types) {
  file_path <- paste0(input_folder, type, "_Blood-BM.tabular")
  output_file <- paste0(output_folder, "volcano_", type, "_blood_vs_bm.png")

  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    next
  }

  df <- read.table(file_path, header = TRUE, sep = "\t")

  # Annotate significance
  df <- df %>%
    mutate(Significance = case_when(
      FDR < FDR_threshold & logFC > logFC_threshold ~ "Upregulated",
      FDR < FDR_threshold & logFC < -logFC_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))

  # Add labels to immune genes and top hits
  df$Label <- NA
  df$Label[df$GeneID %in% immune_genes & df$Significance != "Not Significant"] <- df$GeneID[df$GeneID %in% immune_genes & df$Significance != "Not Significant"]

  top_up <- head(df %>% filter(Significance == "Upregulated") %>% arrange(FDR), 10)$GeneID
  top_down <- head(df %>% filter(Significance == "Downregulated") %>% arrange(FDR), 10)$GeneID
  df$Label[df$GeneID %in% c(top_up, top_down)] <- df$GeneID[df$GeneID %in% c(top_up, top_down)]

  # Volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = Significance)) +
    geom_point(data = df[df$Significance == "Not Significant", ], color = "gray80", size = 1.5, alpha = 0.5) +
    geom_point(data = df[df$Significance == "Upregulated", ], color = "#D7191C", size = 1.8) +
    geom_point(data = df[df$Significance == "Downregulated", ], color = "#2C7BB6", size = 1.8) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed") +
    geom_text_repel(data = df[!is.na(df$Label), ], aes(label = Label), size = 3.2, max.overlaps = 20) +
    scale_color_manual(values = c("Not Significant" = "gray80", "Upregulated" = "#D7191C", "Downregulated" = "#2C7BB6")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      legend.position = "top"
    ) +
    labs(
      title = paste("Volcano Plot:", type, "Monocytes"),
      subtitle = "Blood vs Bone Marrow",
      x = "Log2 Fold Change",
      y = "-Log10(FDR)",
      color = ""
    )

  ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
  message("Saved: ", output_file)
} 