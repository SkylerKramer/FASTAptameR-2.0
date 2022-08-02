# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function plots the results of the motif discovery performed by FSBC
fa_fsbc_motifDiscoveryPlot <- function(motifDF = NULL){
  
  # generate bubble plot where x-axis is the "Rank by Normalized Z-Score", and y-axis is the -log10(p-value)
  p <- ggplot(motifDF) +
    geom_point(aes(x = Rank, y = -log10(P), colour = Motif_Length, size = Motif_Length)) +
    geom_point(aes(x = Rank, y = -log10(P), size = Motif_Length, text = glue::glue("Motif: {Motif}")), shape = 21, colour = "black", stroke = 0.15) +
    scale_colour_viridis_c() +
    scale_size(guide = "none") +
    labs(x = "Rank by Normalized Z-Score", y = "-log10(p)", colour = "Length", title = "Over-Enriched Strings") +
    theme_classic()
  
  # Note, the second geom_point adds the bubble outlines; plotly returns an error when "fill = Motif_Length" (when fill scales when length)
  # The current code does not have this error
  
  # make plot bold
  p <- boldPlots(p)
  
  # make figure
  fig <- plotly::ggplotly(p, tooltip = "text")
  
  # return figure
  return(fig)
}