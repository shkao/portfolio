# Single-Cell RNA Sequencing Analysis Pipeline

This repository contains a comprehensive pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data, implemented in both R and Python. The pipeline covers quality control, filtering, normalization, and data transformation steps, providing a versatile toolkit for researchers working with scRNA-seq data.

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Components](#pipeline-components)
3. [Usage](#usage)
4. [File Descriptions](#file-descriptions)
5. [Dependencies](#dependencies)

## Overview

This scRNA-seq analysis pipeline is designed to process raw sequencing data, perform quality control, filter low-quality cells and genes, normalize the data, and prepare it for downstream analyses. The pipeline is implemented in both R and Python, allowing users to choose their preferred language or compare results between the two implementations.

## Pipeline Components

1. **Quality Control and Filtering**
   - Implemented in R (`01_scRNAseq_QC_and_filtering.R`) and Python (`01_scRNAseq_QC_and_filtering.py`)
   - Filters cells based on library size, number of features, and mitochondrial percentage
   - Provides visualizations of QC metrics (R version)
   - Allows flexible parameter setting for QC thresholds

2. **Normalization and Data Transformation**
   - Implemented in R (`02_scRNAseq_normalization.R`) and Python (`02_scRNAseq_normalization.py`)
   - Normalizes UMI counts using various methods:
     - Deconvolution method (R)
     - Log-normalization (Python)
   - Applies variance stabilization using sctransform (R version)
   - Visualizes data before and after normalization (R version)

## Usage

1. **Quality Control and Filtering:**
   - R: Run `01_scRNAseq_QC_and_filtering.R`
   - Python: Use functions in `01_scRNAseq_QC_and_filtering.py`

2. **Normalization and Data Transformation:**
   - R: Run `02_scRNAseq_normalization.R`
   - Python: Use `log_normalize` function from `02_scRNAseq_normalization.py`

3. Follow the workflow outlined in `scRNA_analysis/README.md` for a step-by-step guide on data processing.

## File Descriptions

- `scRNA_analysis/README.md`: Detailed workflow for the scRNA-seq analysis pipeline
- `01_scRNAseq_QC_and_filtering.R`: R script for quality control and filtering
- `01_scRNAseq_QC_and_filtering.py`: Python script for quality control and filtering
- `02_scRNAseq_normalization.R`: R script for normalization and data transformation
- `02_scRNAseq_normalization.py`: Python script for normalization

## Dependencies

### R Libraries
- scater
- scran
- sctransform
- DropletUtils
- tidyverse
- BiocParallel
- patchwork
- ggplot2
- ensembldb
- AnnotationHub

### Python Libraries
- scanpy
- anndata
- numpy

Please ensure all dependencies are installed before running the scripts. You can install R packages using `install.packages()` or `BiocManager::install()` for Bioconductor packages. For Python, you can use `pip install` to install the required libraries.

This pipeline provides a comprehensive approach to scRNA-seq data analysis, offering flexibility in language choice and methods. It's designed to be adaptable to various experimental designs and can be easily integrated into larger bioinformatics workflows.
