# 01_convert_to_seurat.R

# Description:
# This script converts an .h5ad file into a Seurat RDS object using `sceasy`. 
# It assumes that the source file is located in `data/` and outputs the converted 
# Seurat object also into the `data/` folder.

# Load required libraries
library(Seurat)
library(sceasy)
library(reticulate)

# Define path to your conda environment with scanpy + anndata installed
conda_env <- "~/miniforge3/envs/sc_python"
use_condaenv(condaenv = conda_env, required = TRUE)

# Import scanpy (needed by sceasy)
sc <- import("scanpy", convert = FALSE)

# Convert .h5ad file to Seurat .rds object
sceasy::convertFormat(
  input = "data/Batched_Object.h5ad",
  from = "anndata",
  to = "seurat",
  outFile = "data/Mono_object.rds"
)

# Optional: Load and inspect metadata after conversion
mono_object <- readRDS("data/Mono_object.rds")
head(mono_object@meta.data)

# Save updated Seurat object (clean and minimal for downstream use)
mono_object <- UpdateSeuratObject(mono_object)
mono_object <- DietSeurat(mono_object, assays = "RNA", layers = "counts")
saveRDS(mono_object, "data/Mono_object.rds")