#' Fit Maple multi-sample Bayesian spatial mixture model
#'
#' This function allows you to detect sub-populations and explain membership with relevant covariates in multi-sample spatial transcriptomics experiments.
#' @param seurat_obj An integrated Seurat object 
#' @param K The number of sub-populations to infer. Each should be present in each sample.
#' @param emb Either one of "PCs", "HVGs", or "SVGs" OR a matrix with custom embeddings. If the latter, rows should be sorted as in meta data of Seurat object.
#' @param n_dim The number of dimensions to use if emb is specified as one of "PCs", "HVGs", or "SVGs". Ignored if emb is a matrix of custom embeddings.
#' @param covars Column names of Seurat meta data to use as covariates. If none specified, will fit a global intercept and sample-indicator model for cell type membership probabilities.
#' @param r Spatial smoothing parameter. Should be greater than 0 with larger values enforcing stronger prior spatial association.
#' @param nsim Number of total MCMC iterations to conduct. 
#' @param burn Number of initial MCMC iterations to discard as burn in. The number of saved iterations is nsim-burn
#' @param z_init Initialized cluster allocation vector to aid in MCMC convergence. If NULL z_init will be set using hierarchical clustering. 
#'
#' @keywords spatial transcriptomics Bayesian
#' @import Seurat
#' @importFrom spruce fit_mvn_PG_smooth
#' @importFrom dbarts makeModelMatrixFromDataFrame
#' @importFrom stats model.matrix
#' @export
#' @return A list of MCMC samples, including the MAP estimate of cluster indicators (z)
#' 
fit_maple <- function(seurat_obj,
                      K,
                      emb = "PCs",
                      n_dim = 8,
                      covars = NULL,
                      r = 3,
                      nsim = 2000,
                      burn = 1000,
                      z_init = NULL)
{
  if(emb == "PCs")
  {
    # check PCs are present
    if(!is.null(seurat_obj@reductions$pca))
    {
      emb <- seurat_obj@reductions$pca@cell.embeddings[,1:n_dim]
      rownames(emb) <- rownames(seurat_obj@reductions$pca@cell.embeddings)
    }
    else
    {
      message("Error: No PCA reductions found in the supplied Seurat object. Please add PCA embeddings to seurat_obj$reductions$pca@cell.embeddings.")
      return(NULL)
    }
  }
  else if(emb == "HVGs")
  {
    hvgs <- VariableFeatures(seurat_obj)[1:n_dim]
    emb <- t(seurat_obj@assays$SCT@scale.data[hvgs,])
  }
  else if(emb == "SVGs")
  {
    svgs <- SpatiallyVariableFeatures(seurat_obj)[1:n_dim]
    emb <- t(seurat_obj@assays$SCT@scale.data[svgs,])
  }
  else if(is.matrix(emb))
  {
    emb <- emb
  }
  else
  {
    message("Error: emb should be one of (PCs, HVGs, SVGs) or a matrix of custom embeddings with rows in same order as in seurat_obj@meta.data")
  }
  
  # compile coordinates
  L = length(seurat_obj@images)
  coords = NULL
  for(l in 1:L)
  {
    coords_x_l <- seurat_obj@images[[l]]@coordinates$col
    coords_y_l <- seurat_obj@images[[l]]@coordinates$row
    if(l > 1)
    {
      coords_x_l <- coords_x_l + max(seurat_obj@images[[l-1]]@coordinates$col) + 50
    }
    coords_l <- data.frame(x = coords_x_l,
                           y = coords_y_l)
    rownames(coords_l) <- rownames(seurat_obj@images[[l]]@coordinates)
    coords <- rbind(coords,coords_l)
  }
  
  # check dimensions match
  if(nrow(emb) != nrow(coords))
  {
    message("Number of rows (cells) of embedding must match that of spatial coordinates")
    return(NULL)
  }
  else
  {
    # order coords according to rownames of emb
    coords <- coords[rownames(emb),]
    N <- nrow(emb)
    meta <- seurat_obj@meta.data[rownames(emb),]
  }
  
  # check supplied covariates
  if(is.null(covars))
  {
    print("Note, no supplied covariates. Will use orig.ident as sample indicator.")
    sample <- as.factor(meta[,"orig.ident"])
    W <- as.matrix(stats::model.matrix(~ sample))
    rownames(W) <- rownames(meta)
  }
  else
  {
    if(all(covars %in% colnames(meta)))
    {
      Intercept <- 1
      W <- as.matrix(cbind(Intercept,dbarts::makeModelMatrixFromDataFrame(meta[,covars])))
      rownames(W) <- rownames(meta)
    }
    else
    {
      message("All covariate names must be column names of seurat_obj@meta.data")
      return(NULL)
    }
  }
  
  fit <- spruce::fit_mvn_PG_smooth(Y = emb,
                                   W = W,
                                   coords_df = coords,
                                   K = K,
                                   r = r,
                                   nsim = nsim, 
                                   burn = burn)
  fit$coords = coords
  return(fit)
}

