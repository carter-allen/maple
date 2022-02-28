#' Fit Maple multi-sample Bayesian spatial mixture model
#'
#' This function allows you to detect sub-populations and explain membership with relevant covariates in multi-sample spatial transcriptomics experiments.
#' @param seurat_obj An integrated Seurat object.
#' @param K The number of sub-populations to infer. Each should be present in each sample.
#' @param emb The cell spot embedding to use. Either one of "PCs", "scGNN", "harmony", "HVGs", or "SVGs".
#' @param n_dim The number of dimensions to use. 
#' @param covars Column names of Seurat meta data to use as covariates. If none specified, will fit a global intercept and sample-indicator model for cell type membership probabilities.
#' @param MCAR Logical. Include multivariate CAR random intercepts in gene expression model?
#' @param CAR Logical. Include univariate CAR random intercepts in multinomial gene expression model?
#' @param smooth Logical. Use manual spatial smoothing controlled by r parameter?
#' @param r Spatial smoothing parameter for if smooth == TRUE. Should be greater than 0 with larger values enforcing stronger prior spatial association.
#' @param nsim Number of total MCMC iterations to conduct. 
#' @param burn Number of initial MCMC iterations to discard as burn in. The number of saved iterations is nsim-burn.
#' @param z_init Initialized cluster allocation vector to aid in MCMC convergence. If NULL z_init will be set using hierarchical clustering. 
#'
#' @keywords spatial transcriptomics Bayesian
#' @import Seurat
#' @import spruce
#' @importFrom dbarts makeModelMatrixFromDataFrame
#' @importFrom stats model.matrix cutree
#' @export
#' @return A list of MCMC samples, including the MAP estimate of cluster indicators (z)
#' @examples 
#' \dontrun{
#' brain1 <- LoadData("stxBrain", type = "anterior1")
#' brain2 <- LoadData("stxBrain", type = "anterior2")
#' brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)
#' brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
#' brain <- merge(brain1,brain2)
#' DefaultAssay(brain) <- "SCT"
#' VariableFeatures(brain) <- c(VariableFeatures(brain1),VariableFeatures(brain2))
#' brain <- RunPCA(brain)
#' brain_fit_PCs <- fit_maple(brain,K = 6,emb = "PCs")
#' }
#' 
#' 
fit_maple <- function(seurat_obj,
                      K,
                      emb = "PCs",
                      n_dim = 8,
                      covars = NULL,
                      MCAR = FALSE,
                      CAR = FALSE,
                      smooth = TRUE,
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
      Y <- seurat_obj@reductions$pca@cell.embeddings[,1:n_dim]
      rownames(Y) <- rownames(seurat_obj@reductions$pca@cell.embeddings)
    }
    else
    {
      message("Error: No PCA reductions found in the supplied Seurat object. Please add PCA embeddings to seurat_obj$reductions$pca@cell.embeddings. See Seurat::RunPCA documentation for more details.")
      return(NULL)
    }
  }
  else if(emb == "scGNN")
  {
    # check harmony reductions are present
    if(!is.null(seurat_obj@reductions$scGNN))
    {
      Y <- seurat_obj@reductions$scGNN@cell.embeddings[,1:n_dim]
      rownames(Y) <- rownames(seurat_obj@reductions$scGNN@cell.embeddings)
    }
    else
    {
      message("Error: No scGNN reductions found in the supplied Seurat object. Please add scGNN embeddings to seurat_obj$reductions$scGNN@cell.embeddings. See Seurat::CreateDimReducObject() documentation for more details.")
      return(NULL)
    }
  }
  else if(emb == "harmony")
  {
    # check harmony reductions are present
    if(!is.null(seurat_obj@reductions$harmony))
    {
      Y <- seurat_obj@reductions$harmony@cell.embeddings[,1:n_dim]
      rownames(Y) <- rownames(seurat_obj@reductions$harmony@cell.embeddings)
    }
    else
    {
      message("Error: No harmony reductions found in the supplied Seurat object. Please add harmony embeddings to seurat_obj$reductions$harmony@cell.embeddings. See Seurat::CreateDimReducObject() documentation for more details.")
      return(NULL)
    }
  }
  else if(emb == "HVGs")
  {
    hvgs <- VariableFeatures(seurat_obj)[1:n_dim]
    Y <- t(seurat_obj@assays$SCT@scale.data[hvgs,])
  }
  else if(emb == "SVGs")
  {
    svgs <- SpatiallyVariableFeatures(seurat_obj)[1:n_dim]
    Y <- t(seurat_obj@assays$SCT@scale.data[svgs,])
  }
  else
  {
    message("Error: emb should be one of (PCs, scGNN, harmony, HVGs, or SVGs) or a matrix of custom embeddings with rows in same order as in seurat_obj@meta.data")
  }
  
  # offset images
  coords <- offset_images(seurat_obj)
  
  # check dimensions match
  if(nrow(Y) != nrow(coords))
  {
    message("Number of rows (cells) of embedding must match that of spatial coordinates")
    return(NULL)
  }
  else
  {
    # order coords according to rownames of Y
    coords <- coords[rownames(Y),]
    N <- nrow(Y)
    meta <- seurat_obj@meta.data[rownames(Y),]
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
      W <- as.matrix(cbind(Intercept,dbarts::makeModelMatrixFromDataFrame(as.data.frame(meta[,covars]))))
      rownames(W) <- rownames(meta)
    }
    else
    {
      message("All covariate names must be column names of seurat_obj@meta.data")
      return(NULL)
    }
  }
  
  # dispatch to spruce functions
  print("Dispatching to appropriate model fit function")
  if(emb %in% c("HVGs","SVGs"))
  {
    print("Fitting MSN model to account for skewness of HVGs/SVGs")
    fit <- spruce::fit_msn_PG_smooth(Y = Y,
                                     W = W,
                                     coords_df = coords,
                                     K = K,
                                     r = r,
                                     nsim = nsim, 
                                     burn = burn,
                                     z_init = z_init)
  }
  else
  {
    if((MCAR == TRUE) & (CAR == TRUE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_CAR_MCAR_smooth(Y = Y,
                                                W = W,
                                                coords_df = coords,
                                                K = K,
                                                r = r,
                                                nsim = nsim, 
                                                burn = burn,
                                                z_init = z_init)
    }
    else if((MCAR == TRUE) & (CAR == TRUE) & (smooth == FALSE))
    {
      fit <- spruce::fit_mvn_PG_CAR_MCAR(Y = Y,
                                         W = W,
                                         coords_df = coords,
                                         K = K,
                                         nsim = nsim, 
                                         burn = burn,
                                         z_init = z_init)
    }
    else if((MCAR == TRUE) & (CAR == FALSE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_MCAR_smooth(Y = Y,
                                            W = W,
                                            coords_df = coords,
                                            K = K,
                                            r = r,
                                            nsim = nsim, 
                                            burn = burn,
                                            z_init = z_init)
    }
    else if((MCAR == TRUE) & (CAR == FALSE) & (smooth == FALSE))
    {
      fit <- spruce::fit_mvn_PG_MCAR(Y = Y,
                                     W = W,
                                     coords_df = coords,
                                     K = K,
                                     nsim = nsim, 
                                     burn = burn,
                                     z_init = z_init)
    }
    else if((MCAR == FALSE) & (CAR == TRUE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_CAR_smooth(Y = Y,
                                           W = W,
                                           coords_df = coords,
                                           K = K,
                                           r = r,
                                           nsim = nsim, 
                                           burn = burn,
                                           z_init = z_init)
    }
    else if((MCAR == FALSE) & (CAR == TRUE) & (smooth == FALSE))
    {
      fit <- spruce::fit_mvn_PG_CAR(Y = Y,
                                    W = W,
                                    coords_df = coords,
                                    K = K,
                                    nsim = nsim, 
                                    burn = burn,
                                    z_init = z_init)
    }
    else if((MCAR == FALSE) & (CAR == FALSE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_smooth(Y = Y,
                                       W = W,
                                       coords_df = coords,
                                       K = K,
                                       r = r,
                                       nsim = nsim, 
                                       burn = burn,
                                       z_init = z_init)
    }
    else
    {
      fit <- spruce::fit_mvn_PG(Y = Y,
                                W = W,
                                K = K,
                                nsim = nsim, 
                                burn = burn,
                                z_init = z_init)
    }
  }
  
  fit$coords = coords
  return(fit)
}

