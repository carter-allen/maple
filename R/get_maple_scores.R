#' Get posterior probability scores
#'
#' This function allows you to compute posterior uncertainty and continuous phenotype scores
#' @param fit A list returned by fit_maple()
#'
#' @keywords spatial transcriptomics Bayesian
#' @importFrom spruce get_scores
#' @export
#' @return A fit object (list)
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
#' brain_fit_scores <- get_maple_scores(brain_fit_PCs)
#' }
get_maple_scores <- function(fit)
{
  post_scores <- spruce::get_scores(fit)
  fit$Uncertainty <- post_scores[,"u_score"]
  fit$CPhenotypes <- post_scores[,-ncol(post_scores)]
  return(fit)
}