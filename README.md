# Bioinformatics Portfolio

This repository serves as a showcase for the bioinformatics applications and pipelines I have developed. It demonstrates my skills in data analysis, visualization,
and creating interactive tools for biological research.

## Featured Projects

### 1. Drug Sensitivity and Gene Expression Analysis for Senhwa Biosciences

### Overview

 This Shiny app provides a comprehensive analysis tool for exploring relationships between drug sensitivity and gene expression/copy number across various cancer
 types. It's designed to assist researchers at Senhwa Biosciences in drug discovery and development processes.

 ![Senhwa Biosciences Analysis App](shiny_app/img/app_screenshot.png)

### Key Features

- **Compound and Gene Selection**: Users can select specific compounds and genes of interest.
- **Multiple Data Sources**: Integrates data from PRISM, PRISM AUC, and GDSC AUC databases.
- **Interactive Visualizations**:
  - Scatter plots showing correlations between drug sensitivity and gene features.
  - Venn diagrams comparing cell lines used across different datasets.
- **Statistical Analysis**: Provides correlation coefficients for each cancer subtype.
- **Data Export**: Options to download figures and data tables for further analysis.

### Technical Details

- **Framework**: Built using R Shiny
- **Key Libraries**: tidyverse, ggplot2, ggpubr, ggvenn, shinydashboard
- **Data Handling**: Efficient data manipulation using dplyr and tidyr
- **Visualization**: Custom plots created with ggplot2 and its extensions

### How to Use

 1. Select a compound of interest from the dropdown menu.
 2. Choose a gene and specify whether to analyze its expression or copy number.
 3. Explore the correlation table to identify significant relationships.
 4. Click on a specific cancer subtype to view detailed scatter plots and Venn diagrams.
 5. Download figures or data tables for further analysis or reporting.

 This project showcases my ability to create user-friendly, interactive tools for complex biological data analysis, combining statistical methods with intuitive
 visualizations to aid in scientific research and drug discovery.

### 2. Single-Cell RNA Sequencing Analysis Pipeline

This project demonstrates my expertise in analyzing single-cell RNA sequencing data using R and Bioconductor packages.

#### Key Features

- **Quality Control**: Implemented in `scRNAseq_QC.R`
  - Filtering of low-quality cells based on library size, number of features, and mitochondrial percentage
  - Visualization of QC metrics across samples and experimental groups

- **Normalization and Plotting**: Implemented in `scRNAseq_normalization_plotting.R`
  - Normalization of UMI counts using deconvolution method
  - Variance stabilization using sctransform
  - Visualization of data before and after normalization

#### Technical Details

- **Key Libraries**: scater, scran, sctransform, DropletUtils, ggplot2
- **Data Handling**: Efficient manipulation of large-scale single-cell data
- **Visualization**: Custom plots for QC metrics and normalized data

#### How to Use

1. Run `scRNAseq_QC.R` to perform initial quality control on your single-cell RNA-seq data.
2. Use `scRNAseq_normalization_plotting.R` to normalize the filtered data and generate visualizations.

These scripts provide a basic pipeline for processing and analyzing single-cell RNA sequencing data, showcasing my experience to work with complex biological datasets.
