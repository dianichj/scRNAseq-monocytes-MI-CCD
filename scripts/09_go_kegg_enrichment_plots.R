# 09_go_kegg_enrichment_plots.R

# Description:
# This script generates GO and KEGG enrichment plots for all monocyte subtypes (CM, NCM, IFNCM).

# Load required libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(readr)
library(stringr)

# Output directories
dir.create("plots/go", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/kegg", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# GO DotPlots
# -----------------------------
go_files <- list(
  CM = "data/GO_CM.csv",
  NCM = "data/GO_NCM.csv",
  IFNCM = "data/GO_IFNCM.csv"
)

for (name in names(go_files)) {
  df <- read_csv(go_files[[name]]) %>%
    mutate(log10p = -log10(over_represented_pvalue))

  top_terms <- df %>%
    group_by(Component, Tissue) %>%
    slice_min(order_by = over_represented_pvalue, n = 5) %>%
    ungroup() %>%
    mutate(
      term = str_trunc(term, width = 60),
      Tissue = factor(Tissue, levels = c("Blood", "BM")),
      term = fct_reorder(term, log10p)
    )

  p <- ggplot(top_terms, aes(x = Tissue, y = term)) +
    geom_point(aes(size = numDEInCat, color = log10p),
               position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3)) +
    facet_wrap(~ Component, scales = "free_y", nrow = 1) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(10, 20)) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Top 5 GO Terms per Category in", name, "Monocytes"),
      x = "Tissue",
      y = "GO Term",
      size = "Number of DEGs",
      color = expression(-log[10](p~value))
    ) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      strip.text = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    )

  ggsave(paste0("plots/go/go_dotplot_", tolower(name), ".png"), plot = p, width = 22, height = 14, dpi = 600)
}

# -----------------------------
# KEGG BarPlots
# -----------------------------
kegg_files <- list(
  CM = "data/KEGG_CM.csv",
  NCM = "data/KEGG_NCM.csv",
  IFNCM = "data/KEGG_IFNCM.csv"
)

for (name in names(kegg_files)) {
  df <- read_csv(kegg_files[[name]]) %>%
    mutate(log10_padj = -log10(p_adjust_over_represented)) %>%
    mutate(Pathway_Name = fct_reorder(Pathway_Name, log10_padj))

  p <- ggplot(df, aes(x = log10_padj, y = Pathway_Name, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Blood" = "#D7263D", "Bone Marrow" = "#1B9AAA")) +
    labs(
      title = paste("Top Enriched KEGG Pathways in", name, "Monocytes"),
      x = expression(-log[10]("adjusted p-value")),
      y = "KEGG Pathway",
      fill = "Tissue"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )

  ggsave(paste0("plots/kegg/kegg_barplot_", tolower(name), ".png"), plot = p, width = 15, height = 12, dpi = 600)
}