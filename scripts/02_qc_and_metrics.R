# 02_qc_and_metrics.R

# Description:
# This script performs quality control on the Seurat object, calculates
# common metrics (mitochondrial/ribosomal content), and generates violin plots.

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

# Load processed Seurat object
mono_object <- readRDS("data/Mono_object.rds")

# Calculate QC metrics
mono_object <- PercentageFeatureSet(mono_object, pattern = "^MT-", col.name = "percent.mt")
mono_object <- PercentageFeatureSet(mono_object, pattern = "^RPS", col.name = "percent.rps")
mono_object <- PercentageFeatureSet(mono_object, pattern = "^RPL", col.name = "percent.rpl")

# Extract metadata
md <- mono_object@meta.data %>% as.data.table()

# Calculate means per batch
mean_nCountRNA_per_sample <- md %>% group_by(batch) %>% summarize(mean_nCountRNA = mean(nCount_RNA))
mean_nFeatureRNA_per_sample <- md %>% group_by(batch) %>% summarize(mean_nFeatureRNA = mean(nFeature_RNA))
mean_percent.mt_per_sample <- md %>% group_by(batch) %>% summarize(mean_percent.mt = mean(percent.mt))

# Export summary tables
dir.create("tables", showWarnings = FALSE)
write.csv(mean_nCountRNA_per_sample, "tables/mean_nCountRNA_per_sample.csv", row.names = FALSE)
write.csv(mean_nFeatureRNA_per_sample, "tables/mean_nFeatureRNA_per_sample.csv", row.names = FALSE)
write.csv(mean_percent.mt_per_sample, "tables/mean_percent_mt_per_sample.csv", row.names = FALSE)

# Generate violin plots
dir.create("plots/qc", showWarnings = FALSE)

# Basic QC metrics
VlnPlot(mono_object, features = c("percent.mt", "percent.rps", "percent.rpl"), ncol = 3, pt.size = 0)
ggsave("plots/qc/features_before_subset.jpg", width = 8, height = 4)

VlnPlot(mono_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
ggsave("plots/qc/RNA_before_subset.jpg", width = 8, height = 4)

# Split by batch
VlnPlot(mono_object, features = c("percent.mt", "percent.rps", "percent.rpl"), ncol = 3, pt.size = 0, group.by = "batch")
ggsave("plots/qc/features_by_batch.jpg", width = 8, height = 4)

VlnPlot(mono_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 1, pt.size = 0, group.by = "batch")
ggsave("plots/qc/features_RNA_by_batch.jpg", width = 16, height = 10)

# Save updated Seurat object with QC metrics
saveRDS(mono_object, file = "data/Mono_object.rds")
