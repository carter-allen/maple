# maple <img src="https://carter-allen.github.io/maple.png" align="right" width="150"/>

## Bayesian spatial finite mixture models for identification of cell sub-populations in multi-sample spatial transcriptomics experiments

[![R-CMD-check](https://github.com/carter-allen/maple/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/carter-allen/maple/actions/workflows/R-CMD-check.yaml)

`maple` is a Bayesian statistical framework for modeling multi-sample spatial transcriptomics data sets. `maple` is designed to be used within standard [Seurat](https://satijalab.org/seurat/) workflows, and the user may specify to use principal components (PCs), highly variable genes (HVGs), spatially variable genes (SVGs), or custom cell/cell-spot embeddings such as those generated by [RESEPT](https://github.com/OSU-BMBL/RESEPT). A pre-print publication of `maple` can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.02.28.482296v2).

# Installation 

You may install the stable package version from [CRAN](https://CRAN.R-project.org/package=maple), or the development version as follows:

```
devtools::install_github("carter-allen/maple")
```

Installation time is under 10 minutes for all depedent R packages. `maple` has been tested on Windows, macOS, and Linux systems as detailed [here](https://cran.r-project.org/web/checks/check_results_maple.html)

# Vignettes

- [Multi-Sample Anterior Mouse Brain](https://carter-allen.github.io/stxBrain_multi_maple.html)
- [Multi-Sample Posterior Mouse Brain](https://carter-allen.github.io/stxBrain_posterior_maple.html)
- [Multi-Sample Anterior & Posterior Mouse Brain](https://carter-allen.github.io/stxBrain_all_maple.html)

# Usage

Below are command lines for analysis of a 2-sample anterior mouse brain data set sequenced with the 10X Visium platform. Each model takes < 10 minutes to run for 2,000 iterations on a M1 iMac.

```
# Single cell code
library(Seurat)
library(SeuratData)
library(maple)

# Load data
InstallData("stxBrain")
brain1 <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "anterior2")

brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

brain <- merge(brain1,brain2)
DefaultAssay(brain) <- "SCT"
VariableFeatures(brain) <- c(VariableFeatures(brain1),VariableFeatures(brain2))
brain <- RunPCA(brain)

save(brain, file = paste0(data_dir,"brain_anterior12_merged.RData"))
load(paste0(data_dir,"brain_anterior12_merged.RData"))

# fit maple model 
# PC features
brain_fit_PCs <- fit_maple(brain,K = 6,emb = "PCs")
# plot results
brain$maple_labels_PCs_K6 <- as.factor(brain_fit_PCs$z)
Idents(brain) <- "maple_labels_PCs_K6"
SpatialDimPlot(brain,label = TRUE)

# HVG features
brain_fit_HVGs <- fit_maple(brain,K = 6,emb = "HVGs")
# plot results
brain$maple_labels_HVGs_K6 <- as.factor(brain_fit_HVGs$z)
Idents(brain) <- "maple_labels_HVGs_K6"
SpatialDimPlot(brain,label = TRUE)

# visualize coefficients
spruce::plot_deltas(brain_fit_PCs)

# get scores
post_scores <- spruce::get_scores(brain_fit_PCs)
brain@meta.data <- cbind(brain@meta.data,post_scores)

# plot scores
SpatialFeaturePlot(brain,features = "u_score")
SpatialFeaturePlot(brain,features = "k1_score")
```

# Citation

Allen, C., Chang, Y., Ma, Q., & Chung, D. (2022). MAPLE: A Hybrid Framework for Multi-Sample Spatial Transcriptomics Data. bioRxiv.
