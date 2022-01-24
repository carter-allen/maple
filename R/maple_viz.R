#' Plot tissue architecture labels
#'
#' This function allows you to plot (static or interactive) cell spot labels and uncertainty measures
#' @param fit A list returned by fit_maple()
#' @param pt.size The size of each cell spot point
#' @param interactive Logical parameter controlling static or interactive nature of plot
#' @param shade_uncertainty Logical parameter for shading of cell spots by posterior uncertainty. Must run get_maple_scores() first.
#' @param feature A user-provided feature (e.g., gene of interest) to visualize over tissue spaces instead of sub-population labels.
#'
#' @keywords spatial transcriptomics Bayesian
#' @import ggplot2
#' @import shiny
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom rlang .data
#' @export
#' @return A ggplot object or shiny app window
#' 

maple_viz <- function(fit,
                      pt.size = 1,
                      interactive = FALSE,
                      shade_uncertainty = FALSE,
                      feature = NULL)
{
  options(scipen = 999)
  if(!interactive)
  {
    coords_df = as.data.frame(fit$coords)
    coords_df$Label = factor(fit$z, levels = sort(unique(as.numeric(fit$z))), 
                             labels = paste("Sub-Population", sort(unique(as.numeric(fit$z)))))
    if(!shade_uncertainty)
    {
      if(!is.null(feature))
      {
        g = ggplot(data = coords_df, aes(x = .data$x, 
                                         y = .data$y, 
                                         color = feature)) + 
          geom_point(size = pt.size) + 
          theme_void() + 
          xlab(NULL) + 
          ylab(NULL)
      }
      else
      {
        g = ggplot(data = coords_df, aes(x = .data$x, 
                                         y = .data$y, 
                                         color = .data$Label)) + 
          geom_point(size = pt.size) + 
          theme_void() + 
          xlab(NULL) + 
          ylab(NULL)
      }
      
    }
    else
    {
      if(!is.null(feature))
      {
        coords_df$Uncertainty = fit$Uncertainty
        g = ggplot(data = coords_df, aes(x = .data$x, 
                                         y = .data$y, 
                                         color = feature, 
                                         alpha = .data$Uncertainty)) + 
          geom_point(size = pt.size) + 
          theme_void() + 
          xlab(NULL) + 
          ylab(NULL)
      }
      else
      {
        coords_df$Uncertainty = fit$Uncertainty
        g = ggplot(data = coords_df, aes(x = .data$x, 
                                         y = .data$y, 
                                         color = .data$Label, 
                                         alpha = .data$Uncertainty)) + 
          geom_point(size = pt.size) + 
          theme_void() + 
          xlab(NULL) + 
          ylab(NULL)
      }
    }
    return(g)
  }
  else
  {
    shinyApp(
      ui <- fluidPage(
        plotlyOutput("distPlot")
      ), 
      server <- function(input, output) {
        output$distPlot <- renderPlotly({
          coords_df = as.data.frame(fit$coords)
          coords_df$Label = factor(fit$z, levels = sort(unique(as.numeric(fit$z))), 
                                   labels = paste("Sub-Population", sort(unique(as.numeric(fit$z)))))
          if(!shade_uncertainty)
          {
            if(!is.null(feature))
            {
              g = ggplot(data = coords_df, aes(x = .data$x, 
                                               y = .data$y, 
                                               color = feature)) + 
                geom_point(size = pt.size) + 
                theme_void() + 
                xlab(NULL) + 
                ylab(NULL)
            }
            else
            {
              g = ggplot(data = coords_df, aes(x = .data$x, 
                                               y = .data$y, 
                                               color = .data$Label)) + 
                geom_point(size = pt.size) + 
                theme_void() + 
                xlab(NULL) + 
                ylab(NULL)
            }
            
          }
          else
          {
            if(!is.null(feature))
            {
              coords_df$Uncertainty = fit$Uncertainty
              g = ggplot(data = coords_df, aes(x = .data$x, 
                                               y = .data$y, 
                                               color = feature, 
                                               alpha = .data$Uncertainty)) + 
                geom_point(size = pt.size) + 
                theme_void() + 
                xlab(NULL) + 
                ylab(NULL)
            }
            else
            {
              coords_df$Uncertainty = fit$Uncertainty
              g = ggplot(data = coords_df, aes(x = .data$x, 
                                               y = .data$y, 
                                               color = .data$Label, 
                                               alpha = .data$Uncertainty)) + 
                geom_point(size = pt.size) + 
                theme_void() + 
                xlab(NULL) + 
                ylab(NULL)
            }
          }
        })
      }
    )
  }
}