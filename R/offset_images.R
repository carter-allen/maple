#' Offset multiple Seurat images
#'
#' Internal function for defining coordinates with from multi-sample Seurat objects
#' @param seurat_obj A Seurat object
#'
#' @keywords spatial transcriptomics Bayesian
#' @import Seurat
#' @export
#' @return A coordinate data frame
#' 
offset_images <- function(seurat_obj)
{
  L = length(seurat_obj@images)
  coords = NULL
  for(l in 1:L)
  {
    coords_x_l <- scale(seurat_obj@images[[l]]@coordinates$col,scale = F)
    coords_y_l <- scale(seurat_obj@images[[l]]@coordinates$row,scale = F)
    coords_x_l <- coords_x_l + abs(min(coords_x_l))
    if(l > 1)
    {
      coords_x_l <- coords_x_l + max_x + mean(coords_x_l)
    }
    coords_l <- data.frame(x = coords_x_l,
                           y = coords_y_l)
    rownames(coords_l) <- rownames(seurat_obj@images[[l]]@coordinates)
    coords <- rbind(coords,coords_l)
    max_x <- max(coords_x_l)
  }
  return(coords)
}
