# 03_sctransform_and_integration.R

# Description:
# Applies SCTransform normalization and performs integration with Harmony, CCA, and scVI.

# Load libraries
library(Seurat)
library(SeuratWrappers)
library(reticulate)

# Load Seurat object
mono_object <- readRDS("data/Mono_object.rds")

# Optional: Split RNA layer by batch before SCTransform
mono_object[["RNA"]] <- split(mono_object[["RNA"]], f = mono_object$batch)

# Run SCTransform
mono_object <- SCTransform(
  mono_object,
  vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl")
)

# Run PCA on SCT assay
mono_object <- RunPCA(mono_object, npcs = 30, verbose = TRUE, assay = "SCT")

# Save intermediate object
saveRDS(mono_object, "data/Mono_object.rds")

# -----------------------------
# Harmony integration
# -----------------------------
mono_object <- IntegrateLayers(
  assay = "SCT",
  object = mono_object,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = TRUE
)

# -----------------------------
# CCA integration
# -----------------------------
mono_object <- IntegrateLayers(
  assay = "SCT",
  object = mono_object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  normalization.method = "SCT",
  verbose = TRUE
)

# -----------------------------
# scVI integration (requires Python env)
# -----------------------------
conda_env <- "~/miniforge3/envs/scvi_final"
use_condaenv(condaenv = conda_env, required = TRUE)

mono_object <- NormalizeData(mono_object)
mono_object <- FindVariableFeatures(mono_object)
mono_object <- ScaleData(mono_object)
mono_object <- RunPCA(mono_object)

mono_object <- IntegrateLayers(
  assay = "RNA",
  object = mono_object,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.scvi",
  verbose = TRUE,
  conda_env = conda_env
)

# Save final integrated object
saveRDS(mono_object, "data/Mono_object.rds")