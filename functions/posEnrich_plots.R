# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates a bar plot to show the average enrichment of non-reference residues at each position of a reference sequence
fa_enrich_avgSequenceBar <- function(
    dataPath = NULL, refSeq = NULL, enrichRange = c(0,5), seqType = "Nucleotide", modList = "",
    lowCol = "red2", midCol = "gold", highCol = "yellow", breakpoints = c()
){
  
  # nucleotides or amino acids
  if(seqType == "Nucleotide"){
    charList <- c("A", "C", "G", "T", "U")
  } else if(seqType == "AminoAcid"){
    charList <- c(
      "*",
      "A", "C", "F", "G", "I", "L", "M", "P", "V",
      "W", "N", "Q", "S", "T", "Y",
      "H", "K", "R", "D", "E"
    )
  } else{
    
    message("seqType must be one of Nucleotide or AminoAcid!")
    return(NULL)
  }
  
  # check if any modifications are requested
  if(modList != ""){
    
    # split modifications by newline
    mods_formatted <- strsplit(modList, split = "\n") %>% unlist()
    
    # check if any additions
    if(sum(nchar(mods_formatted) == 1) > 0){
      additions <- mods_formatted[nchar(mods_formatted) == 1]
      
      # add new characters
      charList <- charList %>%
        append(additions)
    }
    
    # check if any replacements
    if(sum(nchar(mods_formatted) > 1) > 0){
      replacements <- mods_formatted[nchar(mods_formatted) > 1]
      
      # make replacements
      charList <- charList %>%
        replace(
          .,
          . == replacements %>% strsplit(., ",") %>% lapply(function(x) x[1]) %>% unlist(),
          replacements %>% strsplit(., ",") %>% lapply(function(x) x[2]) %>% unlist()
        )
    }
  }
  
  # only keep unique characters (in case duplicates are added)
  charList <- unique(charList)
  
  # length of reference sequence
  seqLength <- nchar(refSeq)
  
  # replace ambiguities in reference sequence with X
  refSeq <- gsub(
    paste0("[^", paste(charList, collapse = ""), "]") %>% gsub("\\*", "\\\\*", .),
    "X",
    refSeq
  )
  
  # read data and only keep sequences with the same length as the reference sequence; also convert ambiguities to X
  hmData <- read.csv(dataPath) %>%
    dplyr::mutate(., seqs = as.character(seqs)) %>%
    dplyr::filter(., nchar(seqs) == seqLength) %>%
    dplyr::mutate(
      seqs = gsub(
        paste0("[^", paste(charList, collapse = ""), "]") %>% gsub("\\*", "\\\\*", .),
        "X",
        seqs
      )
    )
  
  # initialize matrix; rows are charLists and cols are the positions of the WT seq
  hmDF <- matrix(
    data = NA,
    nrow = length(charList), ncol = seqLength,
    dimnames = list(charList, unlist(stringr::str_split(refSeq, pattern = "")))
  )
  
  # iterate through each position of the WT seq
  for(i in 1:seqLength){
    # extract all charLists that are not WT
    muts <- data.frame(table(substr(hmData$seqs, i, i))) %>%
      dplyr::filter(Var1 != substr(refSeq, i, i) & Var1 != "X") %>%
      dplyr::pull(Var1) %>%
      as.character()
    
    # iterate through each of the substitutions
    for(mut in muts) {
      # get mean enrichment score from each entry
      hmDF[mut, i] <- hmData %>%
        dplyr::filter(substr(seqs, i, i) == mut) %>%
        dplyr::summarize(meanEnrichment = mean(enrichment_ba)) %>%
        dplyr::pull(meanEnrichment)
    }
  }
  
  # split reference seq. into single characters appended with "_{position}"
  refSeq_formatted <- paste0(unlist(stringr::str_split(refSeq, pattern = "")), "_", 1:seqLength)
  
  # make df that repeats charList (amino acids or nucleotides) for every position of the reference sequence; add column for enrichment
  dta <- expand.grid(Characters = charList, refSeq = refSeq_formatted) %>%
    dplyr::mutate(Enrichment = as.vector(unlist(hmDF)))
  
  # floor and ceiling with enrichment values, depending on user preference
  dta$Enrichment[dta$Enrichment < min(enrichRange)] <- min(enrichRange)
  dta$Enrichment[dta$Enrichment > max(enrichRange)] <- max(enrichRange)
  
  # group by refSeq, which is the column with formatted positions from the reference sequence; get mean enrichment per character
  dta_average <- dta %>%
    dplyr::group_by(refSeq) %>%
    dplyr::summarize(AvEnrich = mean(Enrichment, na.rm = T))
  
  # add min and max possible values, sort, and omit duplicates
  breakpoints_mod <- c(1, as.numeric(breakpoints), length(refSeq_formatted)) %>% sort() %>% unique()
  
  # get mean of AvEnrich between each pair of breakpoints
  bpMeans <- lapply(
    2:length(breakpoints_mod),
    function(x) dta_average$AvEnrich[breakpoints_mod[x-1]:(breakpoints_mod[x] - 1)] %>% mean()
  ) %>% unlist()
  
  # make bar plot of average enrichment per position
  p <- ggplot2::ggplot(
    dta_average, ggplot2::aes(refSeq, AvEnrich, fill = AvEnrich, text = glue::glue(
      "
    WT charList: {refSeq}
    Average Enrichment: {AvEnrich}
    "
    ))
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_discrete(labels = gsub("_.*", "", refSeq_formatted)) +
    ggplot2::scale_fill_gradient2(
      midpoint = quantile(dta_average$AvEnrich, 0.65, na.rm = T),
      low = lowCol, mid = midCol, high = highCol
    ) +
    ggplot2::geom_vline(xintercept = breakpoints_mod, linetype = "dotted") +
    ggplot2::labs(
      x = "Reference Sequence",
      y = "Avg. Enrichment of Non-Reference Residues",
      title = "Average Enrichment per Position",
      fill = "Average\nEnrichment"
    ) +
    ggplot2::theme_classic()
  
  # iterate through breakpoints to add average bar to each segment
  for(i in 2:length(breakpoints_mod)){
    p <- p + geom_segment(
      x = breakpoints_mod[i-1], xend = breakpoints_mod[i],
      y = bpMeans[i-1], yend = bpMeans[i-1]
    )
  }
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text")
  
  return(fig)
}

#' This function creates a heat map to show the average enrichment of non-reference residues at each position of a reference sequence,
#' further resolved by all possible mutated residues
fa_enrich_heatMap <- function(
    dataPath = NULL, refSeq = NULL, enrichRange = c(0,5), seqType = "Nucleotide", modList = "",
    lowCol = "red2", midCol = "gold", highCol = "yellow"
){
  
  # nucleotides or amino acids
  if(seqType == "Nucleotide"){
    
    charList <- c("A", "C", "G", "T", "U")
  } else if(seqType == "AminoAcid"){
    
    charList <- c(
      "*",
      "A", "C", "F", "G", "I", "L", "M", "P", "V",
      "W", "N", "Q", "S", "T", "Y",
      "H", "K", "R", "D", "E"
    )
  } else{
    
    message("seqType must be one of Nucleotide or AminoAcid!")
    return(NULL)
  }
  
  # check if any modifications are requested
  if(modList != ""){
    
    # split modifications by newline
    mods_formatted <- strsplit(modList, split = "\n") %>% unlist()
    
    # check if any additions
    if(sum(nchar(mods_formatted) == 1) > 0){
      additions <- mods_formatted[nchar(mods_formatted) == 1]
      
      # add new characters
      charList <- charList %>%
        append(additions)
    }
    
    # check if any replacements
    if(sum(nchar(mods_formatted) > 1) > 0){
      replacements <- mods_formatted[nchar(mods_formatted) > 1]
      
      # make replacements
      charList <- charList %>%
        replace(
          .,
          . == replacements %>% strsplit(., ",") %>% lapply(function(x) x[1]) %>% unlist(),
          replacements %>% strsplit(., ",") %>% lapply(function(x) x[2]) %>% unlist()
        )
    }
  }
  
  # only keep unique characters (in case duplicates are added)
  charList <- unique(charList)
  
  # length of reference sequence
  seqLength <- nchar(refSeq)
  
  # replace ambiguities in reference sequence with X
  refSeq <- gsub(
    paste0("[^", paste(charList, collapse = ""), "]") %>% gsub("\\*", "\\\\*", .),
    "X",
    refSeq
  )
  
  # read data and only keep sequences with the same length as the reference sequence; also convert ambiguities to X
  hmData <- read.csv(dataPath) %>%
    dplyr::mutate(., seqs = as.character(seqs)) %>%
    dplyr::filter(., nchar(seqs) == seqLength) %>%
    dplyr::mutate(
      seqs = gsub(
        paste0("[^", paste(charList, collapse = ""), "]") %>% gsub("\\*", "\\\\*", .),
        "X",
        seqs
      )
    )
  
  # initialize matrix; rows are charLists and cols are the positions of the WT seq
  hmDF <- matrix(
    data = NA,
    nrow = length(charList), ncol = seqLength,
    dimnames = list(charList, unlist(stringr::str_split(refSeq, pattern = "")))
  )
  
  # iterate through each position of the WT seq
  for(i in 1:seqLength){
    # extract all charLists that are not WT
    muts <- data.frame(table(substr(hmData$seqs, i, i))) %>%
      dplyr::filter(Var1 != substr(refSeq, i, i) & Var1 != "X") %>%
      dplyr::pull(Var1) %>%
      as.character()
    
    # iterate through each of the substitutions
    for(mut in muts) {
      # get mean enrichment score from each entry
      hmDF[mut, i] <- hmData %>%
        dplyr::filter(substr(seqs, i, i) == mut) %>%
        dplyr::summarize(meanEnrichment = mean(enrichment_ba)) %>%
        dplyr::pull(meanEnrichment)
    }
  }
  
  # split reference seq. into single characters appended with "_{position}"
  refSeq_formatted <- paste0(unlist(stringr::str_split(refSeq, pattern = "")), "_", 1:seqLength)
  
  # make df that repeats charList (amino acids or nucleotides) for every position of the reference sequence; add column for enrichment
  dta <- expand.grid(Characters = charList, refSeq = refSeq_formatted) %>%
    dplyr::mutate(Enrichment = as.vector(unlist(hmDF)))
  
  # floor and ceiling with enrichment values, depending on user preference
  dta$Enrichment[dta$Enrichment < min(enrichRange)] <- min(enrichRange)
  dta$Enrichment[dta$Enrichment > max(enrichRange)] <- max(enrichRange)
  
  # make heatmap as raster
  p <- ggplot2::ggplot(dta, ggplot2::aes(refSeq, Characters, fill = Enrichment)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(midpoint = quantile(dta$Enrichment, 0.65, na.rm = T),
                                  low = lowCol, mid = midCol, high = highCol) +
    ggplot2::scale_x_discrete(labels = gsub("_.*", "", refSeq_formatted)) +
    ggplot2::labs(
      x = "Reference Sequence",
      y = "Residues",
      title = "Enrichment Heat Map",
      fill = "Average\nEnrichment"
    ) +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}
