# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates two histograms of distance: 1) for unique sequences and 2) for all reads
fa_distance_histogram <- function(distanceData = NULL, querySequence = NULL){
  
  # distance histogram with unique sequences
  p1 <- ggplot2::ggplot(distanceData, ggplot2::aes(Distance)) +
    ggplot2::geom_bar(fill = "skyblue", colour = "black") +
    ggplot2::xlab("Distance from Query") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_classic()
  
  # make bold text
  p1 <- boldPlots(p = p1) %>%
    plotly::ggplotly() %>%
    plotly::add_text(
      x = mean(distanceData$Distance),
      y = max(ggplot2::ggplot_build(p1)$data[[1]]$ymax),
      text = querySequence, textfont = list(size = 12)
    )
  
  if(!("Reads" %in% colnames(distanceData))){
    return(p1)
  } else{
    # distance histogram with total reads
    p2 <- distanceData %>%
      dplyr::group_by(Distance) %>%
      dplyr::summarise(TotalReads = sum(Reads)) %>%
      
      ggplot2::ggplot(ggplot2::aes(x = Distance, y = TotalReads)) +
      ggplot2::geom_bar(stat = "identity", fill = "skyblue", colour = "black") +
      ggplot2::xlab("Distance from Query") + ggplot2::ylab("Total number of reads") +
      ggplot2::theme_classic()
    
    # make bold text
    p2 <- boldPlots(p = p2) %>%
      plotly::ggplotly() %>%
      plotly::add_text(
        x = mean(distanceData$Distance),
        y = max(ggplot2::ggplot_build(p2)$data[[1]]$ymax),
        text = querySequence, textfont = list(size = 12)
      )
    
    # make interactive figure
    fig <- plotly::subplot(p1, p2, titleY = T, titleX = T, nrows = 2, margin = 0.1, shareX = TRUE) %>%
      plotly::layout(title = "<b>Distance Histograms</b>", showlegend = F, margin = list(t = 50))
    
    # return interactive histogram
    return(fig)
  }
}
