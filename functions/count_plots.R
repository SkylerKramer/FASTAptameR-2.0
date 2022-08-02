# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function makes a line plot showing the relationship between rank and reads
fa_count_rpr <- function(countData = NULL, minReads = NULL, maxRanks = NULL){
  
  # make data.frame by breaking 'id' column into 'Reads' and 'Ranks'; filter by minReads and maxRanks
  seqCounts <- data.frame(
    Reads = countData$id %>% strsplit(., "-") %>% lapply(., `[[`, 2) %>% unlist(.) %>% as.numeric(.),
    SequenceRank = countData$id %>% strsplit(., "-") %>% lapply(., `[[`, 1) %>% unlist(.) %>% gsub(">", "", .) %>% as.numeric(.)
  ) %>%
    dplyr::filter(Reads > minReads & SequenceRank < maxRanks)
  
  # make plot
  p <- ggplot2::ggplot(seqCounts, ggplot2::aes(x = SequenceRank, y = Reads, group = 1)) +
    ggplot2::geom_line(colour = "skyblue", size = 2) +
    ggplot2::labs(x = "Ranks of unique sequences", y = "Total reads per unique sequence", title = "Read count for each rank") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return figure
  return(fig)
}

#' This function makes two histograms of sequence lengths: 1) for unique sequences and 2) for all reads
fa_count_histogram <- function(countData = NULL){
  # make histogram of sequence lengths
  p1 <- ggplot2::ggplot(countData, ggplot2::aes(Length)) +
    ggplot2::geom_bar(fill = "skyblue", colour = "black") +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    ggplot2::labs(x = "Sequence length", y = "Number of unique sequences") +
    ggplot2::theme_classic()
  
  p2 <- countData %>%
    dplyr::group_by(Length) %>%
    dplyr::summarise(TotalReads = sum(Reads)) %>%
    
    ggplot2::ggplot(ggplot2::aes(x = Length, y = TotalReads)) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", colour = "black") +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    ggplot2::labs(x = "Sequence length", y = "Total number of reads") +
    ggplot2::theme_classic()
  
  # make bold plots
  p1 <- boldPlots(p = p1)
  p2 <- boldPlots(p = p2)
  
  # make interactive figure
  fig <- plotly::subplot(plotly::ggplotly(p1), plotly::ggplotly(p2), titleY = T, titleX = T, nrows = 2, margin = 0.1) %>%
    plotly::layout(title = "<b>Sequence-Length Histogram</b>", showlegend = F, margin = list(t = 50))
  
  # return interactive plot
  return(fig)
}

#' This function bins sequences by read counts
#' The bar plot shows bins (x-axis) against fraction of the population (y-axis), colored by number of unique reads
fa_count_binnedAbundance <- function(countData = NULL, useSingleton = TRUE, breaks = c(10,100,1000)){
  
  # initialize breaks if not specified by user
  if(length(breaks) == 1){
    if(breaks == ""){
      breaks <- c(10,100,1000)
    }
  }
  
  # remove breaks outside of the range of the data
  breaks <- breaks[breaks > min(countData$Reads) & breaks < max(countData$Reads)]
  
  # add min and max read values to breaks vector, sort, filter for unique values
  breaks_mod <- c(0, breaks, max(countData$Reads)) %>% sort() %>% unique()
  
  # get vector of labels for cut
  labelNames <- vector()
  for(i in 2:length(breaks_mod)){
    if(i == 2 & useSingleton){
      labelNames[i-1] <- paste0("1 < Reads < ", breaks_mod[i])
    } else{
      labelNames[i-1] <- paste0(breaks_mod[i-1], " <= Reads < ", breaks_mod[i])
    }
  }
  
  # get vector for factor labels
  if(useSingleton){
    factorNames <- c("Singleton", labelNames)
  } else{
    factorNames <- labelNames
  }
  
  # label reads by break points
  countData_formatted <- countData %>%
    dplyr::mutate(BinnedReads = as.character(cut(Reads, breaks = breaks_mod, labels = labelNames, right = TRUE))) %>%
    dplyr::mutate(
      BinnedReads = dplyr::case_when(
        useSingleton & Reads == 1 ~ "Singleton",
        TRUE ~ BinnedReads
      )
    ) %>%
    
    # group by bins, count number of unique sequences per bin, and count number of reads per bin
    dplyr::group_by(BinnedReads) %>%
    dplyr::summarise(
      n = dplyr::n(),
      TotalReads = sum(Reads)
    ) %>%
    
    # get fraction of reads in each bin; sort by fraction
    dplyr::mutate(Fraction = round(TotalReads / sum(TotalReads), 2)) %>%
    dplyr::arrange(desc(Fraction)) %>%
    
    # reorder bin levels
    dplyr::mutate(BinnedReads = factor(BinnedReads, levels = factorNames))
  
  # make barplot
  p <- ggplot2::ggplot(countData_formatted, ggplot2::aes(BinnedReads, Fraction, fill = n)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_gradient(low = "skyblue", high = "salmon") +
    ggplot2::labs(
      x = "No. Reads", y = "Fraction of Population",
      fill = "No. Unique\nSequences",
      title = "Binned Sequence Abundance"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}
