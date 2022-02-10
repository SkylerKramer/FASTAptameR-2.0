# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function makes text in ggplot2 objects bold
boldPlots <- function(p = NULL){
  # adjust theme
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size=14, face= "bold", colour= "black" ),
      axis.title.x = ggplot2::element_text(size=14, face="bold", colour = "black"),    
      axis.title.y = ggplot2::element_text(size=14, face="bold", colour = "black"),    
      axis.text.x = ggplot2::element_text(size=12, face="bold", colour = "black"), 
      axis.text.y = ggplot2::element_text(size=12, face="bold", colour = "black"),
      strip.text.x = ggplot2::element_text(size = 10, face="bold", colour = "black" ),
      strip.text.y = ggplot2::element_text(size = 10, face="bold", colour = "black"),
      axis.line.x = ggplot2::element_line(color="black", size = 0.3),
      axis.line.y = ggplot2::element_line(color="black", size = 0.3),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=0.3)
    )
  
  # return modified plot
  return(p)
}

#' This function makes a line plot showing the relationship between rank and reads
fa_count_rpr <- function(countData = NULL, minReads = NULL, maxRanks = NULL){
  # make data.frame by breaking 'id' column into 'Reads' and 'Ranks'; filter by minReads and maxRanks
  seqCounts <- data.frame(Reads = countData$id %>%
                            strsplit(., "-") %>% lapply(., `[[`, 2) %>%
                            unlist(.) %>% as.numeric(.),
                          SequenceRank = countData$id %>%
                            strsplit(., "-") %>% lapply(., `[[`, 1) %>%
                            unlist(.) %>% gsub(">", "", .) %>% as.numeric(.)) %>%
    .[which(.$Reads > minReads & .$SequenceRank < maxRanks),]
  
  # make plot
  p <- ggplot2::ggplot(seqCounts, ggplot2::aes(x = SequenceRank, y = Reads, group = 1)) +
    ggplot2::geom_line(colour = "skyblue", size = 2) +
    ggplot2::labs(x = "Ranks of unique sequences", y = "Total reads per unique sequence") +
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
    ggplot2::xlab("Sequence length") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_classic()
  
  p2 <- countData %>%
    dplyr::group_by(Length) %>%
    dplyr::summarise(TotalReads = sum(Reads)) %>%
    
    ggplot2::ggplot(ggplot2::aes(x = Length, y = TotalReads)) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", colour = "black") +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    ggplot2::xlab("Sequence length") + ggplot2::ylab("Total number of reads") +
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
    ggplot2::labs(x = "No. Reads", y = "Fraction of Population",
                  fill = "No. Unique\nSequences",
                  title = "Binned Sequence Abundance") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

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

#' This function creates two histograms of distance: 1) for unique sequences and 2) for all reads
fa_distance_histogram <- function(distanceData = NULL, querySequence = NULL){
  # distance histogram with unique sequences
  p1 <- ggplot2::ggplot(distanceData, ggplot2::aes(Distance)) +
    ggplot2::geom_bar(fill = "skyblue", colour = "black") +
    ggplot2::xlab("Distance from Query") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_classic()
  
  # make bold text
  p1 <- boldPlots(p = p1)
  
  p1 <- p1 %>%
    plotly::ggplotly() %>%
    plotly::add_text(x = mean(distanceData$Distance), y = max(ggplot2::ggplot_build(p1)$data[[1]]$ymax),
                     text = querySequence, textfont = list(size = 12))
  
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
    p2 <- boldPlots(p = p2)
    
    p2 <- p2 %>%
      plotly::ggplotly() %>%
      plotly::add_text(x = mean(distanceData$Distance), y = max(ggplot2::ggplot_build(p2)$data[[1]]$ymax),
                       text = querySequence, textfont = list(size = 12))
    
    # make interactive figure
    fig <- plotly::subplot(p1, p2, titleY = T, titleX = T, nrows = 2, margin = 0.1, shareX = TRUE) %>%
      plotly::layout(title = "<b>Distance Histograms</b>", showlegend = F, margin = list(t = 50))
    
    # return interactive histogram
    return(fig)
  }
}

#' This function creates a bar plot to show in how many rounds sequences persist (e.g., 1 round, 2 rounds, etc.)
fa_enrich_seqPersistence <- function(fastaInputs = NULL, minReads = 0){
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
    dplyr::filter(Reads > minReads) %>%
    dplyr::select(seqs) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::group_by(Freq) %>%
    dplyr::summarise(., seqCount = dplyr::n_distinct(.))
  
  # make bar plot
  p <- ggplot2::ggplot(seqPersist, ggplot2::aes(Freq, seqCount, text = glue::glue(
    "
  No. Rounds Detected: {Freq}
  No. Unique Reads: {seqCount}
  "
  ))) +
    ggplot2::geom_bar(stat = "identity", colour = "black", fill = "skyblue") +
    ggplot2::labs(x = "No. Rounds Detected", y = "No. Unique Reads", title = "Sequence Persistence Analysis") +
    ggplot2::scale_x_continuous(breaks = unique(seqPersist$Freq)) +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text")
  
  # return interactive figure
  return(fig)
}

#' This function creates a histogram of fold changes
fa_enrich_histogram <- function(df = NULL){
  # get population information from log2E column
  popInfo <- substr(colnames(df)[1], 7, 8)
  
  # rename column to appease ggplot2 grammar
  colnames(df) <- "log2E"
  
  # make plot
  p <- ggplot2::ggplot(df, ggplot2::aes(log2E)) +
    ggplot2::geom_histogram(color = "black", fill = "skyblue", bins = 30) +
    ggplot2::labs(x = "log2(Enrichment)", 
                  y = "Number of Unique Sequences",
                  title = paste0("log2E Histogram - ", substr(popInfo, 1, 1), ":", substr(popInfo, 2, 2))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45))
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates a scatter plot of RPMs in two populations
fa_enrich_scatter <- function(df = NULL){
  # save RPM names
  rpmNames <- colnames(df)[1:2]
  
  # rename RPM columns to appease plotly grammar
  colnames(df)[1:2] <- c("RPM1", "RPM2")
  
  # add epsilon factor if RPM is equal to 0
  df <- df %>%
    dplyr::mutate(
      RPM1 = ifelse(RPM1 == 0, RPM1 + 0.01, RPM1),
      RPM2 = ifelse(RPM2 == 0, RPM2 + 0.01, RPM2)
    )
  
  # make plot
  p <- ggplot2::ggplot(df, ggplot2::aes(RPM1, RPM2, text = seqs)) +
    ggplot2::geom_point(colour = "skyblue", alpha = 0.5) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = rpmNames[1], y = rpmNames[2], title = paste0(rpmNames[1], " vs ", rpmNames[2])) +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates an RA plot from RPMs in two populations; epsilon = 0.1
fa_enrich_ra <- function(df = NULL){
  # save RPM names
  rpmNames <- colnames(df)[1:2]
  
  # rename RPM columns to appease plotly grammar
  colnames(df)[1:2] <- c("RPM1", "RPM2")
  
  # add epsilon factor if RPM is equal to 0
  df <- df %>%
    dplyr::mutate(
      RPM1 = ifelse(RPM1 == 0, RPM1 + 0.01, RPM1),
      RPM2 = ifelse(RPM2 == 0, RPM2 + 0.01, RPM2)
    )
  
  # add fold change and average log RPM
  df <- df %>%
    dplyr::mutate(
      R = log2(RPM2 / RPM1),
      A = 0.5 * log2(RPM2 * RPM1)
    )
  
  # make plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = A, y = R, text = seqs)) +
    ggplot2::geom_point(colour = "skyblue", alpha = 0.5) +
    ggplot2::labs(x = "Average log2(RPM)", y = "Fold Change", title = "Enrichment RA Plot") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates box plots of enrichment per cluster
fa_enrich_clusterBoxplots <- function(df = NULL){
  # save population information
  populations <- gsub(".*_", "", colnames(df)[4]) %>% toupper() %>% strsplit("") %>% unlist()
  
  # rename columns for plotting functionality
  colnames(df) <- c("Seqs", "Cluster", "SeedDistance", "Enrichment")
  
  # get longest cluster number string
  padLength <- df$Cluster %>% as.character() %>% nchar() %>% max()
  
  # left-pad cluster number and sort
  df_format <- df %>%
    dplyr::mutate(Cluster = stringr::str_pad(Cluster, padLength, pad = "0", side = "left")) %>%
    dplyr::arrange(Cluster) %>%
    dplyr::mutate(Cluster = factor(Cluster, levels = unique(Cluster)))
  
  # filter for seed sequences
  df_seeds <- df_format %>% filter(SeedDistance == 0)
  
  # make box plots
  p <- ggplot2::ggplot() +
    ggplot2::geom_boxplot(
      data = df_format,
      ggplot2::aes(Cluster, Enrichment),
      fill = "skyblue", outlier.shape = NA
    ) +
    ggplot2::geom_point(
      data = df_seeds,
      ggplot2::aes(Cluster, Enrichment, text = glue::glue(
        "
      Cluster no.: {Cluster}
      Enrichment: {Enrichment}
      Seed: {Seqs}
      "
      )),
      colour = "red", pch = 18
    ) +
    ggplot2::labs(x = "Cluster No.", y = "Enrichment",
                  title = paste0("Sequence Enrichment by Cluster - ", populations[1], ":", populations[2])) +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text", dynamicTicks = TRUE)
  
  # return interactive figure
  return(fig)
}

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
    charList <- c("*",
                  "A", "C", "F", "G", "I", "L", "M", "P", "V",
                  "W", "N", "Q", "S", "T", "Y",
                  "H", "K", "R", "D", "E")
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
  hmDF <- matrix(data = NA, nrow = length(charList), ncol = seqLength,
                 dimnames = list(charList,
                                 unlist(stringr::str_split(refSeq, pattern = ""))))
  
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
    ggplot2::labs(x = "Reference Sequence",
                  y = "Residues",
                  title = "Enrichment Heat Map",
                  fill = "Average\nEnrichment") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates a PCA scatter plot (first two components) of a kmer matrix
fa_clusterDiversity_kmerPCA <- function(clusterFile = NULL, kmerSize = 3, topClusters = 10, keepNC = T){
  # format clustered file
  clusterDF <- fa_formatInput(fastaInput = clusterFile, population = NULL)
  
  # convert ambiguous bases to X
  clusterDF$seqs <- gsub("[^ACGTU]", "X", clusterDF$seqs)
  
  # convert NC values to NA
  clusterDF[clusterDF == "NC"] <- NA
  clusterDF$Cluster <- as.numeric(as.character(clusterDF$Cluster))
  
  # only keep the top clusters
  if(keepNC){
    clusterDF <- clusterDF[which(clusterDF$Cluster <= topClusters | is.na(clusterDF$Cluster)),]
  } else{
    clusterDF <- clusterDF[which(clusterDF$Cluster <= topClusters),]
  }
  
  # compute kmer matrix
  kmerMatrix <- clusterDF$seqs %>%
    as.character(.) %>%
    strsplit(., split = "") %>%
    kmer::kcount(., k = kmerSize, named = F)
  rownames(kmerMatrix) <- NULL
  
  # pca with kmer matrix
  kmerPCA <- prcomp(kmerMatrix, scale. = F)
  
  # make plot
  kmerPCA_plot <- factoextra::fviz_pca_ind(
    kmerPCA,
    addEllipses = F, geom = c("point"),
    col.ind = as.factor(clusterDF$Cluster), legend.title = "Cluster"
  ) +
    ggplot2::theme_classic()
  
  # make bold text
  kmerPCA_plot <- boldPlots(p = kmerPCA_plot)
  
  # make interactive figure
  fig <- plotly::ggplotly(kmerPCA_plot)
  
  # return interactive PCA plot
  return(fig)
}

#' This function creates line plots of metadata per cluster
fa_clusterDiversity_metaplot <- function(diversityDF = NULL){
  # remove non-clustered sequences and convert cluster to character
  diversityDF <- diversityDF %>%
    dplyr::filter(Cluster != "NC") %>%
    dplyr::mutate(Cluster = as.character(Cluster))
  
  # line plot of unique sequences per cluster
  seqsPlot <- plotly::plot_ly(diversityDF, type = "scatter", mode = "lines", x = ~Cluster, y = ~TotalSequences) %>%
    plotly::layout(xaxis = list(title = "<b>Cluster Number</b>", categoryorder = "array", categoryarray = ~Cluster, zeroline = T),
                   yaxis = list(title = "<b>Unique Sequences</b>"))
  
  # line plot of total sequences per cluster
  readsPlot <- plotly::plot_ly(diversityDF, type = "scatter", mode = "lines", x = ~Cluster, y = ~TotalReads) %>%
    plotly::layout(xaxis = list(title = "<b>Cluster Number</b>", categoryorder = "array", categoryarray = ~Cluster, zeroline = T),
                   yaxis = list(title = "<b>Total Reads</b>"))

  # line plot of average LED to seed per cluster
  ledPlot <- plotly::plot_ly(diversityDF, type = "scatter", mode = "lines", x = ~Cluster, y = ~AverageLED) %>%
    plotly::layout(xaxis = list(title = "<b>Cluster Number</b>", categoryorder = "array", categoryarray = ~Cluster, zeroline = T),
                   yaxis = list(title = "<b>Average LED to Cluster Seed</b>"))

  # make interactive figure from subplots
  fig <- plotly::subplot(seqsPlot, readsPlot, ledPlot, titleY = TRUE, shareX = TRUE, nrows = 3) %>%
    plotly::layout(title = "<b>Cluster Metaplots</b>", showlegend = FALSE)
  
  # return interactive figure
  return(fig)
}

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