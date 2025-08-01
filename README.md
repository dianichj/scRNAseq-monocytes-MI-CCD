# scRNAseq-monocytes-MI-CCD

This repository contains the analysis code and workflows used to characterize human monocyte subsets in myocardial infarction (MI) and chronic coronary disease (CCD) patients using single-cell RNA sequencing (scRNA-seq). The analysis combines Galaxy-based preprocessing with downstream processing in R, and was developed as part of my doctoral research at the University of Freiburg.

## 🔬 Overview of Analysis Workflow

The scRNA-seq data processing and analysis were conducted in two main phases:

### 1. **Low-level preprocessing**
- Raw BCL files were demultiplexed into FASTQ using **Cell Ranger mkfastq** and **bcl2fastq**.
- Gene-barcode matrices were generated using **Cell Ranger count**.
- Initial quality control and preprocessing were performed on [Galaxy Europe](https://usegalaxy.eu):
  - Filtering of cells with high mitochondrial content or abnormal gene counts.
  - Merging of samples using **AnnData** tools.
  - Normalization based on the [Galaxy Scanpy tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html).

### 2. **Downstream analysis in R and Galaxy**
- AnnData files were converted to Seurat objects using **sceasy** and **reticulate**.
- Processed using **Seurat v5.0.3**, **SCTransform**, **UMAP**, **Louvain clustering**, and **CCA** for batch correction.
- Cell type annotation was done based on canonical markers.
- **Pseudobulk analysis** was carried out by:
  - Re-uploading the AnnData object to Galaxy.
  - Creating pseudobulk count matrices using **Decoupler**.
  - Running **edgeR** for differential gene expression (DEGs).
  - Performing GO/KEGG enrichment with **GOseq**, correcting for gene length bias.

## 🔄 Galaxy Workflows Used

The following Galaxy workflows were used as part of this analysis:

- **Pseudobulk differential expression analysis (Decoupler + edgeR)**  
  [Galaxy IWC: pseudobulk workflow](https://github.com/galaxyproject/iwc/tree/main/workflows/scRNAseq/pseudobulk-worflow-decoupler-edger)  
  *Developed by Diana Chiang Jurado (ORCID: 0000-0002-5857-1477), Pavankumar Videm (ORCID: 0000-0002-5192-126X), and Pablo Moreno (ORCID: 0000-0002-9856-1679).*

- **GOseq enrichment analysis**  
  [Galaxy IWC: goseq workflow](https://github.com/galaxyproject/iwc/tree/main/workflows/transcriptomics/goseq)  
  *Developed by Amirhossein Naghsh Nilchi.*

## 📁 Repository Structure

```
scRNAseq-monocytes-MI-CCD/
│
├── scripts/                  # Modular R scripts (cleaned, numbered)
│   ├── 01_load_and_prepare.R
│   ├── 02_qc_filtering.R
│   ├── 03_integration_harmony.R
│   ├── 04_integration_cca.R
│   ├── 05_integration_scvi.R
│   ├── 06_annotation.R
│   ├── 07_dimensionality_reduction.R
│   ├── 08_dotplots_markers.R
│   ├── 09_GO_KEGG_enrichment.R
│   ├── 10_heatmap_and_featureplots.R
│   ├── 11_violin_and_barplots.R
│   ├── 12_volcano_plots.R
│   ├── 13_combined_expression_plots.R
│
├── data/                    # Processed Seurat objects or sample data
│   └── Mono_object.rds      # Not publicly available until manuscript publication
│
├── README.md                # Main project documentation
├── LICENSE                  # MIT License
└── .gitignore               # Used to hide results and data files from GitHub
```

🔒 *Note: `figures/` and `results/` directories are excluded until manuscript can be made publicly available.*
