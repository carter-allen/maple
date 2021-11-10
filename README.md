# maple

Bayesian spatial finite mixture models for identification of cell sub-populations in multi-sample spatial transcriptomics experiments

# Installation 

```
devtools::install_github("carter-allen/maple")
```

# Requirements

`maple` depends on our existing package `spruce`, which contains a suite of functions for single sample spatial transcriptomics data analysis. `spruce` can be installed via:

```
devtools::install_github("carter-allen/spruce")
```

`maple` also depends heavily on `Seurat` for data pre-processing and visualization.

# Usage

Below are command lines for analysis of a 2-sample anterior mouse brain data set sequenced with the 10X Visium platform.

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