#' Get posterior probability scores
#'
#' This function allows you to compute posterior uncertainty and continuous phenotype scores
#' @param fit A list returned by fit_maple()
#'
#' @keywords spatial transcriptomics Bayesian
#' @importFrom spruce get_scores
#' @export
#' @return A fit object (list)
#' 
get_maple_scores <- function(fit)
{
  post_scores <- spruce::get_scores(fit)
  fit$Uncertainty <- post_scores[,"u_score"]
  fit$CPhenotypes <- post_scores[,-ncol(post_scores)]
  return(fit)
}