# Zika scRNA-seq antiviral transcriptional response analysis

This repository contains the R code used for the manuscript:

**Exploratory strain-associated patterns of antiviral transcriptional responses to Zika virus exposure in developing human neural tissue**

## Overview

The analysis is a secondary descriptive reanalysis of public transcriptomic datasets. The primary dataset is human fetal brain-derived neural tissue exposed ex vivo to Mock, ZIKV-BR, ZIKV-FSS/FSS13025, or IFNβ. Additional public datasets are used only as contextual comparisons of IFN/ISG-associated transcriptional activation.

The analysis should be interpreted as exploratory and descriptive. The primary fetal neural dataset contains one biological sample per condition, so cells are not treated as biological replicates and exploratory differential expression summaries are not interpreted as formal condition-level inference.

## Public datasets

- **GSE238140**: primary fetal neural scRNA-seq dataset
- **GSE230571**: contextual Vero/moDC scRNA-seq dataset
- **GSE97919**: contextual bulk RNA-seq dataset

Large raw GEO files are not included in this repository. Download the relevant GEO supplementary files and place them under the folders described in `data/README.md`.

## Repository structure

```text
.
├── scripts/
│   ├── 01_prepare_10x_data.R
│   ├── 02_qc_normalization_umap_clustering.R
│   ├── 03_celltype_annotation.R
│   ├── 04_exploratory_celltype_DE.R
│   ├── 05_progenitor_volcano_ISG_heatmap.R
│   ├── 06_ZIKV_IFNb_correlation_sensitivity.R
│   ├── 07_make_supplementary_tables.R
│   ├── 08_transcriptional_response_summary.R
│   ├── 09_bulk_RNAseq_contextual_GSE97919.R
│   └── run_all.R
├── data/
│   └── README.md
├── results/
├── figures/
├── tables/
├── supplementary_tables/
├── session_info/
├── install_packages.R
└── README.md
```

## Installation

Install R and the required packages, then run:

```r
source("install_packages.R")
```

Main R dependencies include Seurat, Matrix, tidyverse, dplyr, ggplot2, patchwork, readr, tibble, tidyr, purrr, scales, ragg, writexl, data.table, stringr, DESeq2, pheatmap, magick, and Cairo.

## How to run

From the repository root:

```r
source("scripts/run_all.R")
```

Alternatively, run the scripts sequentially from `01_` to `09_`.

## Expected outputs

The scripts generate processed R objects, exploratory differential expression tables, supplementary tables, session information, and publication-style figures in the following folders:

- `results/`
- `figures/`
- `tables/`
- `supplementary_tables/`
- `session_info/`

## Reproducibility notes

- Scripts are designed to be run from the repository root.
- Paths are relative to the repository root.
- Large raw data and generated outputs are intentionally excluded from version control.
- Session information is written by the analysis scripts where applicable.

## Citation

Please cite the associated manuscript if using this code.
