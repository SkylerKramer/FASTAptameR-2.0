# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates a line plot of cluster abundance across populations
fa_clusterEnrichTracker <- function(clusterEnrichDF = NULL){
  
  # make line plot
  p <- clusterEnrichDF %>%
    ggplot2::ggplot(., ggplot2::aes(
      Population, TotalRPM,
      colour = Seeds,
      text = glue::glue(
        "
        Population: {Population}
        File Name: {FileName}
        Seed: {Seeds}
        Cluster: {Cluster}
        Unique Sequences: {TotalSequences}
        Total Reads: {TotalReads}
        Total RPM: {TotalRPM}
        "
      )
    )) +
    ggplot2::geom_line(group = 1, size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_continuous(breaks = unique(clusterEnrichDF$Population)) +
    ggplot2::labs(x = "Population", y = "Total RPM", title = "Seed Tracker") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text") %>% plotly::layout(legend = list(orientation = 'h', y = -0.2))
  
  # return interactive figure
  return(fig)
}
