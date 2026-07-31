# GD3 Plasma and CSF Proteomics Analysis

## Overview
This repository contains R analysis scripts for the manuscript:
**"Plasma and CSF Proteomics: Neuroinflammatory and Lysosomal Dysregulation in Gaucher Disease Type 3"**

**Authors:** Pavithra Krishnaswami, Mikhail Levit, Luis Concha-Marambio, Carly M. Farris, Isabela Batsu, Katherine Klinger, Bailin Zhang, S. Pablo Sardi, Can Kayatekin

## Study Description
Comprehensive proteomic profiling of plasma and cerebrospinal fluid (CSF) from 9 GD3 patients and age-matched healthy controls using the Olink® Explore HT platform, quantifying 5,416 unique proteins.

## Script Contents
The main analysis script includes:
- Data loading and quality control (Olink NPX processing)
- Statistical analysis (two-sided t-tests, FDR and Bonferroni corrections)
- Pathway enrichment analysis (Gene Set Enrichment Analysis using clusterProfiler)
- Neuroinflammatory pathway prioritization
- Visualization (volcano plots, boxplots, pathway enrichment plots)

## Requirements

### R Version
R >= 4.4.12

### Required Packages
- OlinkAnalyze (v3.0)
- clusterProfiler (v4.8.0)
- ggplot2
- dplyr
- tidyr
- org.Hs.eg.db (v3.17.0)

### Installation
```R
# Install from CRAN
install.packages(c("ggplot2", "dplyr", "tidyr"))

# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

# Install OlinkAnalyze from GitHub
devtools::install_github("Olink-Proteomics/OlinkAnalyze")
