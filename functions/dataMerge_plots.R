# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates a bar plot to show in how many rounds sequences persist (e.g., 1 round, 2 rounds, etc.)
fa_dataMerge_seqPersistence <- function(fastaInputs = NULL){
  
  # initialize empty data frame
  seqPersist <- data.frame()
  
  # iterate through all inputs
  for(i in 1:length(fastaInputs)){
    
    # bind main df with newly read sequences
    seqPersist <- rbind(
      seqPersist, 
      fa_formatInput(fastaInput = fastaInputs[i]) %>% dplyr::select(seqs, Reads)
    )
  }
  
  # count number of sequences in 1 round, 2 rounds, etc.
  seqPersist <- seqPersist %>%
    dplyr::select(seqs) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::group_by(Freq) %>%
    dplyr::summarise(., seqCount = dplyr::n_distinct(.))
  
  # make bar plot
  p <- ggplot2::ggplot(
    seqPersist,
    ggplot2::aes(
      Freq, seqCount,
      text = glue::glue(
        "
        No. Rounds Detected: {Freq}
        No. Unique Reads: {seqCount}
        "
      )
    )
  ) +
    ggplot2::geom_bar(stat = "identity", colour = "black", fill = "skyblue") +
    ggplot2::labs(x = "No. Populations Detected", y = "No. Unique Reads", title = "Sequence Persistence Analysis") +
    ggplot2::scale_x_continuous(breaks = unique(seqPersist$Freq)) +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text")
  
  # return interactive figure
  return(fig)
}

#' This function creates an UpSetR plot
fa_dataMerge_UpSetR <- function(fastaInputs = NULL, fastaNames = NULL){
  
  # initialize empty list
  sequenceList <- list()
  
  # create list of lists for sequences
  for(i in 1:length(fastaInputs)){
    sequenceList[[i]] <- fa_formatInput(fastaInput = fastaInputs[i]) %>% dplyr::pull(seqs)
    names(sequenceList)[i] <- fastaNames[i]
  }
  
  # make UpSetR plot (not a ggplot2 object, so it cannot be made interactive)
  p <- UpSetR::upset(
    UpSetR::fromList(sequenceList),
    order.by = "freq",
    main.bar.color = "skyblue", sets.bar.color = "skyblue",
    mainbar.y.label = "Sequence Intersections", sets.x.label = "Sequences Per Set", 
    text.scale = 2, point.size = 3, line.size = 1,
    mb.ratio = c(0.65,0.35)
  )
  
  # return plot
  return(p)
}
