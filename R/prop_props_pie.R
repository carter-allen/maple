#' Plot grouped pie charts of mixture component proportions
#'
#' This function allows you to visualize the relative abundance of sub-populations after running fit_maple()
#' @param fit A list returned by fit_maple()
#' @param group A column name of fit$W or a grouping vector of length N (nrow(fit$W))
#'
#' @keywords spatial transcriptomics Bayesian
#' @import dplyr
#' @import ggplot2
#' @export
#' @return A ggplot object
#' 
plot_props_pie <- function(fit, group)
{
  if(length(group) == 1)
  {
    al_df = data.frame(z = fit$z,group = fit$W[,group])
  }
  else if(length(group) == nrow(fit$W))
  {
    group = as.factor(group)
    al_df = data.frame(z = fit$z,group = group)
  }
  else
  {
    message("Error: group should be either a column name of fit$W, a column index of fit$W, or a length N categorical vector coercible to a factor.")
    return()
  }
  
  al_df_sum = al_df %>%
    group_by(z,group) %>%
    summarize(Freq = n())
  
  al_df_group = al_df %>%
    group_by(group) %>%
    summarize(n_group = n())
  
  al_df_sum = inner_join(al_df_sum,al_df_group,by = "group") %>%
    mutate(prop = Freq/n_group) %>%
    mutate(z = as.factor(z), group = as.factor(group))
  
  g <- ggplot(al_df_sum, aes(x = "", y = prop, fill = z)) + 
    geom_bar(stat = "identity", width = 1, color = "black") + 
    theme_void() + 
    coord_polar("y",start = 0) + 
    facet_wrap(~group)
  
  return(g)
}