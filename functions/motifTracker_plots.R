# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates a line plot of motif abundance across populations
fa_motif_motifTrackerPlot <- function(targetDF = NULL){
  # make base plot
  if("Alias" %in% colnames(targetDF)){
    p <- ggplot2::ggplot(targetDF, ggplot2::aes(
      Population, TotalRPM,
      colour = Alias,
      text = glue::glue(
        "
      Population: {Population}
      File Name: {FileName}
      Alias: {Alias}
      Motif: {Motif}
      Total Reads: {TotalReads}
      Total RPM: {TotalRPM}
      "
      )
    ))
  } else{
    p <- ggplot2::ggplot(targetDF, ggplot2::aes(
      Population, TotalRPM,
      colour = Motif,
      text = glue::glue(
        "
      Population: {Population}
      File Name: {FileName}
      Motif: {Motif}
      Total Reads: {TotalReads}
      Total RPM: {TotalRPM}
      "
      )
    ))
  }
  
  # create line plot
  p <- p +
    ggplot2::geom_line(group = 1, size = 1) +
    ggplot2::geom_point(size = 5) +
    ggplot2::scale_x_continuous(breaks = unique(targetDF$Population)) +
    ggplot2::labs(x = "Population", y = "Total RPM", title = "Motif Tracker") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text") %>% plotly::layout(legend = list(orientation = 'h', y = -0.2))
  
  # return interactive figure
  return(fig)
}

#' This function creates a line plot of sequence abundance across populations
fa_motif_sequenceTrackerPlot <- function(targetDF = NULL){
  # make base plot
  if("Alias" %in% colnames(targetDF)){
    p <- ggplot2::ggplot(targetDF, ggplot2::aes(
      Population, RPM,
      colour = Alias,
      text = glue::glue(
        "
      Population: {Population}
      File Name: {FileName}
      Alias: {Alias}
      Sequence: {seqs}
      Rank: {Rank}
      Reads: {Reads}
      RPM: {RPM}
      "
      )
    ))
  } else{
    p <- ggplot2::ggplot(targetDF, ggplot2::aes(
      Population, RPM,
      colour = seqs,
      text = glue::glue(
        "
      Population: {Population}
      File Name: {FileName}
      Sequence: {seqs}
      Rank: {Rank}
      Reads: {Reads}
      RPM: {RPM}
      "
      )
    ))
  }
  
  # make line plot
  p <- p +
    ggplot2::geom_line(group = 1, size = 1) +
    ggplot2::geom_point(size = 5) +
    ggplot2::scale_x_continuous(breaks = unique(targetDF$Population)) +
    ggplot2::labs(x = "Population", y = "RPM", title = "Sequence Tracker") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text") %>% plotly::layout(legend = list(orientation = 'h', y = -0.2))
  
  # return interactive figure
  return(fig)
}
