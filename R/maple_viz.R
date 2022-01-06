#' Plot tissue architecture labels
#'
#' This function allows you to plot (static or interactive) cell spot labels and uncertainty measures
#' @param fit A list returned by fit_maple()
#' @param pt.size The size of each cell spot point
#' @param interactive Logical parameter controlling static or interactive nature of plot
#' @param shade_uncertainty Logical parameter for shading of cell spots by posterior uncertainty. Must run get_maple_scores() first.
#'
#' @keywords spatial transcriptomics Bayesian
#' @import ggplot2
#' @import shiny
#' @import plotly
#' @export
#' @return A ggplot object or shiny app window
#' 

maple_viz <- function(fit,
                      pt.size = 1,
                      interactive = FALSE,
                      shade_uncertainty = FALSE)
{
  options(scipen = 999)
  if(!interactive)
  {
    coords_df = as.data.frame(fit$coords)
    coords_df$Label = as.factor(fit$z)
    if(!shade_uncertainty)
    {
      g = ggplot(data = coords_df, aes(x = .data$x, 
                                       y = .data$y, 
                                       color = .data$Label)) + 
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
    return(g)
  }
  else
  {
    require(shiny)
    shinyApp(
      ui <- fluidPage(
        plotlyOutput("distPlot")
      ), 
      server <- function(input, output) {
        output$distPlot <- renderPlotly({
          coords_df = as.data.frame(fit$coords)
          coords_df$Label = as.factor(fit$z)
          if(!shade_uncertainty)
          {
            ggplot(data = coords_df, aes(x = .data$x, 
                                         y = .data$y, 
                                         color = .data$Label)) + 
              geom_point(size = pt.size) + 
              theme_void() + 
              xlab(NULL) + 
              ylab(NULL)
          }
          else
          {
            coords_df$Uncertainty = round(fit$Uncertainty,2)
            ggplot(data = coords_df, aes(x = .data$x, 
                                         y = .data$y, 
                                         color = .data$Label, 
                                         alpha = .data$Uncertainty)) + 
              geom_point(size = pt.size) + 
              theme_void() + 
              xlab(NULL) + 
              ylab(NULL)
          }
        })
      }
    )
  }
}