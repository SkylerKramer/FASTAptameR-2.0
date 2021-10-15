### Skyler Kramer - FASTAptameR 2.0

##### imports ####
library(dplyr)
library(purrr)
library(ggplot2)
library(plotly)

##### support functions #####
fa_formatInput <- function(fastaInput = NULL, population = NULL){
  # read data, save IDs, convert sequences to data.frame
  countData <- readLines(fastaInput)
  
  seqIDs <- as.character(countData[seq(1, length(countData), 2)])
  countData <- data.frame(seqs = countData[seq(2, length(countData), 2)], stringsAsFactors = F)
  
  # process seqIDs to make named columns
  idDF <- gsub('>', '', seqIDs) %>%
    strsplit(., '-') %>%
    unlist(.) %>%
    matrix(., nrow = length(seqIDs), byrow = T) %>%
    as.data.frame(.)
  
  # rename columns and optionally append the population value
  colnames(idDF) <- c('Rank', 'Reads', 'RPM', 'Cluster', 'ClusterRank', 'SeedDistance')[1:ncol(idDF)]
  if(!is.null(population)){
    colnames(idDF) <- paste0(colnames(idDF), '.', population)
  }
  
  # treat values as numerics
  idDF[,1] <- as.numeric(as.character(idDF[,1]))
  idDF[,2] <- as.numeric(as.character(idDF[,2]))
  idDF[,3] <- as.numeric(as.character(idDF[,3]))
  
  # return the data.frame
  return(cbind(countData, idDF))
}
fa_formatOutput <- function(outputData = NULL){
  # create shape for FASTA output
  fastaOutput <- rep(NA, nrow(outputData)*2)
  
  # add sequence IDs and sequences
  fastaOutput[seq(1,length(fastaOutput),2)] <- outputData$id
  fastaOutput[seq(2,length(fastaOutput),2)] <- as.character(outputData$seqs)
  
  # return the data structure
  return(fastaOutput)
}
fa_translate_converge <- function(translateDF = NULL){
  # split 'id' column from data.frame input
  seqIDs <- translateDF$id
  
  idDF <- gsub('>', '', seqIDs) %>%
    strsplit(., '-') %>%
    unlist(.) %>%
    matrix(., nrow = length(seqIDs), byrow = T) %>%
    as.data.frame(.)
  colnames(idDF) <- c("Rank", "Reads", "RPM")
  
  # merge data.frame input with 'id' data.frame; remove Rank and id, which will be remade
  translateDF <- cbind(idDF, translateDF) %>% .[,-c(1,4)]
  
  # group by 'seqs' and sum up 'Reads' and 'RPM'
  translateDF <- as.data.frame(translateDF %>%
                                 dplyr::group_by(seqs) %>%
                                 dplyr::summarise(sumReads = sum(as.numeric(Reads)), sumRPM = sum(as.numeric(RPM))) %>%
                                 dplyr::arrange(-sumReads))
  
  # add 'Rank' and 'id' columns
  translateDF$Rank <- row.names(translateDF)
  translateDF$id <- paste0(">", translateDF$Rank, "-", translateDF$sumReads, "-", translateDF$sumRPM)
  
  # return 'id' and 'seqs'
  return(translateDF[,c(5,1)])
}
fa_translate_mapping <- function(inputChanges = "", translateSelection = "Standard"){
  # convert user input into a usable format
  translationDF_changed <- inputChanges %>% # input string should have 1 codon / translation pair per newline, separated by commas
    gsub("[[:punct:]]", ",", .) %>% # replace all special characters with comma
    gsub(",{2,}", ",", .) %>% # replace consecutive commas with a single one
    gsub(" | \t", "", .) %>% # omit spaces and tabs
    gsub("\\(|\\)", "", .) %>% # omit parentheses
    toupper() %>% # convert all alphabetic characters to uppercase
    strsplit("\n") %>% # split by newline
    unlist() %>% # unlist
    as.data.frame() %>% # convert to data.frame; 1 pair per row
    
    tidyr::separate(., ., c("Codon", "Translation"), ",") %>% # split into multiple columns
    tidyr::drop_na() %>% # drop missing values
    dplyr::mutate(., Codon = gsub("T", "U", Codon)) %>% # replace T with U in Codon column
    dplyr::filter(., nchar(Codon) == 3) %>% # remove codons that are not 3 characters long
    dplyr::filter(., nchar(Translation) == 1) %>% # remove translations that are not 1 characters long
    
    dplyr::mutate(., Codon = as.character(Codon), Translation = as.character(Translation)) # confirm both columns are characters
  
  # select appropriate filepath for JSON file
  filepath <- switch(
    selection,
    "Standard" = "trans/translations_standard.json",
    "Vertebrate mitochondrial" = "trans/translations_vertMito.json",
    "Yeast mitochondrial" = "trans/translations_yeastMito.json",
    "Mold, protozoan, and coelenterate mitochondrial + Mycoplasma / Spiroplasma" = "trans/translations_mold.json",
    "Invertebrate mitochondrial" = "trans/translations_invertMito.json",
    "Ciliate, dasycladacean and Hexamita nuclear" = "trans/translations_cilNucl.json",
    "Echinoderm and flatworm mitochondrial" = "trans/translations_echinoMito.json",
    "Euplotid nuclear" = "trans/translations_euplotidNucl.json",
    "Alternative yeast nuclear" = "trans/translations_altYeastNucl.json",
    "Ascidian mitochondrial" = "trans/translations_ascMito.json",
    "Alternative flatworm mitochondrial" = "trans/translations_altWormMito.json",
    "Blepharisma nuclear" = "trans/translations_blephNucl.json",
    "Chlorophycean mitochondrial" = "trans/translations_chloroMito.json",
    "Trematode mitochondrial" = "trans/translations_tremMito.json",
    "Scenedesmus obliquus mitochondrial" = "trans/translations_sceneMito.json",
    "Pterobranchia mitochondrial" = "trans/translations_pteroMito.json"
  )
  
  # read JSON file with standard genetic code
  translationDF <- jsonlite::fromJSON(filepath) %>% # read JSON file with translation mappings
    dplyr::mutate(., Codon = as.character(Codon), Translation = as.character(Translation)) # confirm both columns are characters
  
  # if input codons and translations are each of length 0, return unmodified data.frame
  if(nrow(translationDF_changed) == 0){
    message("No changes to translation mapping!")
    return(translationDF)
  }
  
  # iterate through dataframe of changes
  for(i in 1:nrow(translationDF_changed)){
    
    # if codon of interest already exists in dataframe, replace the translation
    if(translationDF_changed$Codon[i] %in% translationDF$Codon){
      translationDF[which(translationDF$Codon == translationDF_changed$Codon[i]),2] <- translationDF_changed$Translation[i]
    } else{
      
      # else, append to the data.frame
      translationDF <- rbind(
        data.frame(Codon = translationDF_changed$Codon[i], Translation = translationDF_changed$Translation[i], stringsAsFactors = FALSE),
        translationDF
      )
    }
  }
  
  # return modified mappings
  return(translationDF)
}
fa_motif_format <- function(motifList = NULL, motifType = "Nucleotide"){
  # format motif(s) as a Perl query
  # if "Nucleotide" is specified, consider the degenerate codes
  if(motifType == "Nucleotide"){
    motifList <- motifList %>% toupper(.) %>%
      gsub(" ", "", .) %>%
      gsub('U', 'T', .) %>%
      gsub('R', '[AG]', .) %>%
      gsub('Y', '[CT]', .) %>%
      gsub('W', '[AT]', .) %>%
      gsub('S', '[GC]', .) %>%
      gsub('M', '[AC]', .) %>%
      gsub('K', '[GT]', .) %>%
      gsub('B', '[^A]', .) %>%
      gsub('D', '[^C]', .) %>%
      gsub('H', '[^G]', .) %>%
      gsub('V', '[^T]', .) %>%
      gsub('N', '[ACGT]', .) %>%
      gsub(",", "|", .)
  } else{
    motifList <- motifList %>% toupper(.) %>%
      gsub(" ", "", .) %>%
      gsub(",", "|", .)
  }
  
  # return formatted motif(s)
  return(motifList)
}

##### plot functions #####
fa_count_rpr <- function(countData = NULL, minReads = NULL, maxRanks = NULL){
  # make data.frame by breaking 'id' column into 'Reads' and 'Ranks'; filter by minReads and maxRanks
  seqCounts <- data.frame(Reads = countData$id %>%
                            strsplit(., "-") %>% lapply(., `[[`, 2) %>%
                            unlist(.) %>% as.numeric(.),
                          SequenceRank = countData$id %>%
                            strsplit(., "-") %>% lapply(., `[[`, 1) %>%
                            unlist(.) %>% gsub(">", "", .) %>% as.numeric(.)) %>%
    .[which(.$Reads > minReads & .$SequenceRank < maxRanks),]
  
  # make interactive plot
  plotly::plot_ly(seqCounts, x = ~SequenceRank, y = ~Reads, type = "scatter", mode = "lines") %>%
    plotly::layout(xaxis = list(title = "Ranks of unique sequences", tickfont = list(size = 15), titlefont = list(size = 20)),
                   yaxis = list(title = "Total reads per unique sequence", tickfont = list(size = 15), titlefont = list(size = 20)))
}
fa_count_histogram <- function(countData = NULL){
  # make histogram of sequence lengths
  p1 <- ggplot2::ggplot(countData, ggplot2::aes(Length)) +
    ggplot2::geom_bar(fill = "skyblue", colour = "black") +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    ggplot2::xlab("Sequence length") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
  p2 <- countData %>%
    dplyr::group_by(Length) %>%
    dplyr::summarise(TotalReads = sum(Reads)) %>%
    
    ggplot2::ggplot(ggplot2::aes(x = Length, y = TotalReads)) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", colour = "black") +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    ggplot2::xlab("Sequence length") + ggplot2::ylab("Total number of reads") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
  # return interactive plot
  return(plotly::subplot(plotly::ggplotly(p1), plotly::ggplotly(p2), titleY = T, titleX = T, nrows = 2, margin = 0.1) %>%
    plotly::layout(title = "Sequence-Length Histogram", showlegend = F, height = 600, margin = list(t = 50)))
}
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
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # return interactive figure
  return(plotly::ggplotly(p, height = 600))
}
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
  
  # return interactive figure
  return(plotly::ggplotly(p, height = 500, tooltip = "text") %>% plotly::layout(legend = list(orientation = 'h', y = -0.2)))
}
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
  
  # return interactive figure
  return(plotly::ggplotly(p, height = 500, tooltip = "text") %>% plotly::layout(legend = list(orientation = 'h', y = -0.2)))
}
fa_distance_histogram <- function(distanceData = NULL, querySequence = NULL){
  # distance histogram with unique sequences
  p1 <- ggplot2::ggplot(distanceData, ggplot2::aes(Distance)) +
    ggplot2::geom_bar(fill = "skyblue", colour = "black") +
    ggplot2::xlab("Distance") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
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
      ggplot2::xlab("Distance") + ggplot2::ylab("Total number of reads") +
      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
    
    p2 <- p2 %>%
      plotly::ggplotly() %>%
      plotly::add_text(x = mean(distanceData$Distance), y = max(ggplot2::ggplot_build(p2)$data[[1]]$ymax),
                       text = querySequence, textfont = list(size = 12))
    
    # return interactive histogram
    return(plotly::subplot(p1, p2, titleY = T, titleX = T, nrows = 2, margin = 0.1) %>%
             plotly::layout(title = "Distance Histogram", showlegend = F, height = 600, margin = list(t = 50)))
  }
}
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
  
  # return interactive figure
  return(plotly::ggplotly(p, tooltip = "text"))
}
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
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45))
  
  # return interactive figure
  return(plotly::ggplotly())
}
fa_enrich_scatter <- function(df = NULL){
  # save RPM names
  rpmNames <- colnames(df)[1:2]
  
  # rename RPM columns to appease plotly grammar
  colnames(df)[1:2] <- c("RPM1", "RPM2")
  
  # make interactive figure
  fig <- plotly::plot_ly(df, type = "scatter", mode = "markers", x = ~RPM1, y = ~RPM2, text = ~seqs, color = I("skyblue"), alpha = 0.5) %>%
    plotly::layout(title = paste0(rpmNames[1], " vs ", rpmNames[2]),
                   xaxis = list(type = "log", title = rpmNames[1]), yaxis = list(type = "log", title = rpmNames[2]))
  
  # return interactive figure
  return(fig)
}
fa_enrich_volcano <- function(df = NULL){
  # get population information from log2E column
  popInfo <- substr(colnames(df)[1], 7, 8)
  
  # rename 1st column to appease plotly grammar
  colnames(df)[1] <- "log2E"
  
  # make interactive figure
  fig <- plotly::plot_ly(df, type = "scatter", mode = "markers",
                         x = ~log2E, y = ~statStrength, text = ~seqs) %>%
    plotly::layout(title = paste0("RPM Volcano Plot - ", substr(popInfo, 1, 1), ":", substr(popInfo, 2, 2)),
                   xaxis = list(title = "log2(Enrichment)"),
                   yaxis = list(title = "Statistical Strength"))
  
  # return interactive figure
  return(fig)
}
fa_enrich_ra <- function(df = NULL){
  # save RPM names
  rpmNames <- colnames(df)[1:2]
  
  # rename RPM columns to appease plotly grammar
  colnames(df)[1:2] <- c("RPM1", "RPM2")
  
  # add epsilon factor if RPM is equal to 0
  df <- df %>%
    dplyr::mutate(
      RPM1 = ifelse(RPM1 == 0, RPM1 + 0.1, RPM1),
      RPM2 = ifelse(RPM2 == 0, RPM2 + 0.1, RPM2)
    )
  
  # add fold change and average log RPM
  df <- df %>%
    dplyr::mutate(
      R = log2(RPM2 / RPM1),
      A = 0.5 * log2(RPM2 * RPM1)
    )
  
  # make plot
  p <- ggplot(df, aes(x = A, y = R, text = seqs)) +
    geom_point(colour = "skyblue", alpha = 0.5) +
    labs(x = "Average log2(RPM)", y = "Fold Change", title = "Enrichment RA Plot") +
    theme_classic()
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}
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
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text", dynamicTicks = TRUE)
  
  # return interactive figure
  return(fig)
}
fa_enrich_avgSequenceBar <- function(dataPath = NULL, refSeq = NULL, enrichRange = c(0,5), seqType = "Nucleotide", modList = "",
                                     lowCol = "red2", midCol = "gold", highCol = "yellow"){
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
  
  # group by refSeq, which is the column with formatted positions from the reference sequence; get mean enrichment per character
  dta_average <- dta %>%
    dplyr::group_by(refSeq) %>%
    dplyr::summarize(AvEnrich = mean(Enrichment, na.rm = T))
  
  p <- ggplot2::ggplot(dta_average, ggplot2::aes(refSeq, AvEnrich, fill = AvEnrich, text = glue::glue(
    "
    WT charList: {refSeq}
    Average Enrichment: {AvEnrich}
    "
  ))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_discrete(labels = gsub("_.*", "", refSeq_formatted)) +
    ggplot2::scale_fill_gradient2(midpoint = quantile(dta_average$AvEnrich, 0.65, na.rm = T),
                                  low = lowCol, mid = midCol, high = highCol) +
    ggplot2::labs(x = "Reference Sequence",
                  y = "Avg. Enrichment of Non-Reference Residues",
                  title = "Average Enrichment per Position",
                  fill = "Average\nEnrichment") +
    ggplot2::theme_classic()
  
  return(plotly::ggplotly(p, tooltip = "text"))
}
fa_enrich_heatMap <- function(dataPath = NULL, refSeq = NULL, enrichRange = c(0,5), seqType = "Nucleotide", modList = "",
                              lowCol = "red2", midCol = "gold", highCol = "yellow"){
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
  
  # return interactive figure
  return(plotly::ggplotly(p))
}
fa_clusterDiversity_kmerPCA <- function(clusterFile = NULL, kmerSize = 3, topClusters = 10, keepNC = T){
  # format clustered file
  clusterDF <- fa_formatInput(fastaInput = clusterFile, population = NULL)
  
  # return NULL if the sequences are made of a non-nucleotide alphabet
  # if(sum(grepl("[^ACGTU]", clusterDF$seqs)) != 0){
  #   return(NULL)
  # }
  
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
  kmerPCA_plot <- factoextra::fviz_pca_ind(kmerPCA, addEllipses = F, geom = c("point"), col.ind = as.factor(clusterDF$Cluster)) +
    ggplot2::labs(color="Cluster")
  
  # return interactive PCA plot
  return(plotly::ggplotly(kmerPCA_plot))
}
fa_clusterDiversity_metaplot <- function(diversityDF = NULL){
  # convert "NC" to NA and make Cluster column from factors to numerics
  diversityDF$Cluster <- as.character(diversityDF$Cluster)
  
  # line plot of unique sequences per cluster
  seqsPlot <- plotly::plot_ly(diversityDF, type = "scatter", mode = "lines", x = ~Cluster, y = ~TotalSequences) %>%
    plotly::layout(xaxis = list(title = "Cluster Number", categoryorder = "array", categoryarray = ~Cluster, zeroline = T),
                   yaxis = list(title = "Unique Sequences"))
  
  # line plot of total sequences per cluster
  readsPlot <- plotly::plot_ly(diversityDF, type = "scatter", mode = "lines", x = ~Cluster, y = ~TotalReads) %>%
    plotly::layout(xaxis = list(title = "Cluster Number", categoryorder = "array", categoryarray = ~Cluster, zeroline = T),
                   yaxis = list(title = "Total Reads"))
  
  # line plot of average LED to seed per cluster
  ledPlot <- plotly::plot_ly(diversityDF, type = "scatter", mode = "lines", x = ~Cluster, y = ~AverageLED) %>%
    plotly::layout(xaxis = list(title = "Cluster Number", categoryorder = "array", categoryarray = ~Cluster, zeroline = T),
                   yaxis = list(title = "Average LED to Cluster Seed"))
  
  # return plots
  return(plotly::subplot(seqsPlot, readsPlot, ledPlot, titleY = T, shareX = T, nrows = 3, margin = 0.05) %>%
           plotly::layout(title = "Cluster Metaplots", showlegend = F, height = 600))
}
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
  
  # return interactive figure
  return(plotly::ggplotly(p, height = 800, tooltip = "text") %>% plotly::layout(legend = list(orientation = 'h', y = -0.2)))
}

##### user functions #####
fa_count <- function(dataInput = NULL){
  
  # get start time
  startTime <- Sys.time()
  
  # read FASTQ/A file and keep the sequences in a data.frame
  if(toupper(sub(".*(?=.$)", "", dataInput, perl=T)) == "Q"){
    allSeqs <- data.frame(seqs = toupper(readLines(dataInput) %>% .[seq(from = 2, to = length(.), by = 4)]), stringsAsFactors = F)
  }else{
    allSeqs <- readLines(dataInput) %>%
      .[!grepl('>', .)] %>%
      data.frame(seqs = toupper(.), stringsAsFactors = F)
  }
  
  # count/sort unique sequences; add rank, rpm, and sequence ID to data.frame
  seqCounts <- as.data.frame(allSeqs %>% dplyr::count(seqs, sort = T, name = "Reads"))
  seqCounts$Rank <- as.numeric(row.names(seqCounts))
  seqCounts$RPM <- round(as.numeric(seqCounts$Reads / (nrow(allSeqs) / 1e6)), 2)
  seqCounts$Length <- as.numeric(nchar(seqCounts$seqs))
  seqCounts$id <- paste0(">", seqCounts$Rank, "-", seqCounts$Reads, "-",
                         round(seqCounts$RPM, digits = 2))
  
  # message elapsed time
  message(gsub("Time difference of", "Elapsed time:", capture.output(round(Sys.time() - startTime, 2))))
  
  # return data.frame
  return(seqCounts[,c(6,3,2,4:5,1)])
}
fa_count2 <- function(dataInput = NULL){
  
  # get start time
  startTime <- Sys.time()
  
  # read FASTQ/A file and keep the sequences in a data.frame
  if(toupper(sub(".*(?=.$)", "", dataInput, perl=T)) == "Q"){
    allSeqs <- data.frame(seqs = LaF::get_lines(dataInput, line_numbers = seq(2, LaF::determine_nlines(dataInput), by = 4)))
  }else{
    
    # check if sequences have IDs
    if(readLines(dataInput, n = 1) %>% substr(1, 1) == ">"){
      allSeqs <- data.frame(seqs = LaF::get_lines(dataInput, line_numbers = seq(2, LaF::determine_nlines(dataInput), by = 2)))
    } else{
      allSeqs <- data.frame(seqs = LaF::get_lines(dataInput, line_numbers = seq(1, LaF::determine_nlines(dataInput), by = 1)))
    }
  }
  
  # count/sort unique sequences; add rank, rpm, and sequence ID to data.frame
  seqCounts <- allSeqs %>%
    dplyr::count(seqs, sort = T, name = "Reads") %>%
    tibble::rowid_to_column(var = "Rank") %>%
    dplyr::mutate(
      seqs = as.character(seqs) %>% gsub("\r", "", .),
      RPM = round(Reads / (sum(Reads) / 1e6), 2),
      Length = as.numeric(nchar(seqs)),
      id = paste0(">", Rank, "-", Reads, "-", RPM)
    ) %>%
    dplyr::relocate(id, Rank, Reads, RPM, Length, seqs)
  
  # message elapsed time
  message(gsub("Time difference of", "Elapsed time:", capture.output(round(Sys.time() - startTime, 2))))
  
  # return data.frame
  return(seqCounts)
}
fa_count3 <- function(dataInput = NULL){
  
  # get start time
  startTime <- Sys.time()
  
  # read FASTQ/A file and keep the sequences in a data.frame
  if(toupper(sub(".*(?=.$)", "", dataInput, perl=T)) == "Q"){
    allSeqs <- data.frame(seqs = LaF::get_lines(dataInput, line_numbers = seq(2, LaF::determine_nlines(dataInput), by = 4)))
  }else{
    
    # check if sequences have IDs
    if(readLines(dataInput, n = 1) %>% substr(1, 1) == ">"){
      allSeqs <- data.frame(seqs = LaF::get_lines(dataInput, line_numbers = seq(2, LaF::determine_nlines(dataInput), by = 2)))
    } else{
      allSeqs <- data.frame(seqs = LaF::get_lines(dataInput, line_numbers = seq(1, LaF::determine_nlines(dataInput), by = 1)))
    }
  }
  
  # define function for counting
  countFxn <- function(x){
    data.table::data.table(x)[, .N, keyby = x]
  }
  
  # count/sort unique sequences; add rank, rpm, and sequence ID to data.frame
  seqCounts <- countFxn(allSeqs$seqs) %>% as.data.frame() %>%
    dplyr::rename(seqs = x, Reads = N) %>%
    dplyr::arrange(dplyr::desc(Reads)) %>%
    tibble::rowid_to_column(var = "Rank") %>%
    dplyr::mutate(
      seqs = as.character(seqs) %>% gsub("\r", "", .),
      RPM = round(Reads / (sum(Reads) / 1e6), 2),
      Length = as.numeric(nchar(seqs)),
      id = paste0(">", Rank, "-", Reads, "-", RPM)
    ) %>%
    dplyr::relocate(id, Rank, Reads, RPM, Length, seqs)
  
  # message elapsed time
  message(gsub("Time difference of", "Elapsed time:", capture.output(round(Sys.time() - startTime, 2))))
  
  # return data.frame
  return(seqCounts)
}
fa_translate <- function(fastaInput = NULL, orf = 1, converge = T, inputChanges = "", translateSelection = "Standard"){
  # read FASTA file, save sequence IDs and sequences
  fastaData <- readLines(fastaInput)
  
  # return NULL if not nucleotide sequence
  if(sum(grepl("[^ACGTU]", fastaData[seq(2, length(fastaData), 2)])) != 0){
    return(NA)
  }
  
  # get sequence IDs
  seqIDs <- as.character(fastaData[seq(1, length(fastaData), 2)])
  
  # make all sequences capital, modify by ORF, convert U to T,
  # remove end characters if length is not divisible by 3 (would not be translated),
  # split into groups of three, remove white space at end of split strings
  seqs <- trimws(toupper(as.character(fastaData[seq(2, length(fastaData), 2)])) %>%
                   substring(., orf) %>%
                   gsub('T', 'U', .) %>%
                   ifelse(nchar(.) %% 3 == 1, gsub(".{1}$", "", .), .) %>%
                   ifelse(nchar(.) %% 3 == 2, gsub(".{2}$", "", .), .) %>%
                   gsub("(.{3})", "\\1 ", .), which = "right")
  
  # translate sequences: ambiguous codons to X, remove white spaces after translation
  translationMapping <- fa_translate_mapping(inputChanges = inputChanges, translateSelection = translateSelection)
  
  # translate sequences according to user modifications (if any)
  for(i in 1:nrow(translationMapping)){
    seqs <- gsub(translationMapping[i,1], translationMapping[i,2], seqs)
  }
  
  # finally, remove white spaces between amino acids
  seqs <- gsub(" ", "", seqs)
  
  # optionally merge non-unique amino acid sequences
  if(converge){
    translateDF <- fa_translate_converge(translateDF = data.frame(id = seqIDs, seqs = seqs))
  }else{
    translateDF <- data.frame(id = seqIDs, seqs = seqs)
  }
  
  # breakdown 'id' column into 'Rank', 'Reads', and 'RPM'
  translateDF$Rank <- translateDF$id %>%
    gsub(">", "", .) %>%
    strsplit(., split = "-") %>%
    lapply(., `[[`, 1) %>%
    unlist() %>%
    as.numeric()
  translateDF$Reads <- translateDF$id %>%
    strsplit(., split = "-") %>%
    lapply(., `[[`, 2) %>%
    unlist() %>%
    as.numeric()
  translateDF$RPM <- translateDF$id %>%
    strsplit(., split = "-") %>%
    lapply(., `[[`, 3) %>%
    unlist() %>%
    as.numeric()
  
  # return translated sequences and corresponding data.frame
  # optionally return the number of unique nucleotide sequences for each amino acid sequence if convergence == T
  if(converge){
    translateDF <- merge(translateDF, as.data.frame(table(seqs)), by = "seqs")
    names(translateDF)[names(translateDF) == "Freq"] <- "Unique.Nt.Count"
    translateDF <- translateDF[order(translateDF$Rank),c(2:6, 1)]
  }else{
    return(translateDF[,c(1,3:5,2)])
  }
}
fa_motifSearch <- function(fastaInput = NULL, motif = NULL, highlight = F, partial = F, motifType = "Nucleotide"){
  # read file
  countData <- readLines(fastaInput)
  
  # make countData into data.frame and handle U/T conversion
  countData <- data.frame(id = as.character(countData[seq(1, length(countData), 2)]),
                          seqs = toupper(countData[seq(2, length(countData), 2)]))
  countData$seqs <- gsub('U', 'T', countData$seqs)
  
  # check format of IDs; should have either 3 or 6 components
  if(sum(lengths(strsplit(countData$id, "-")) %% 3) != 0){
    message("Incorrect formatting for sequence IDs!")
    return(NULL)
  }

  # convert case of motif; convert U to T, handle degenerate code in motif
  motif <- fa_motif_format(motifList = motif, motifType = motifType)
  
  # filter countData for sequences with the motif
  if(partial){
    # partial filter uses the OR operator
    countData <- countData[grepl(motif, countData$seqs),]
  }else{
    # full filter repeatedly filters sequences
    fullMatchMotifs <- motif %>%
      paste0("(?=.*", ., ")") %>%
      gsub("\\|", ")(?=.*", .)
    
    countData <- countData[grepl(fullMatchMotifs, countData$seqs, perl = T),]
  }

  # put parentheses around motif
  if(highlight){
    motif <- unlist(strsplit(motif, "\\|"))
    for(i in 1:length(motif)){
      countData$seqs <- gsub(motif[i], paste0('(', motif[i], ')'), countData$seqs)
    }
  }
  
  # decompose the id
  countData$Rank <- countData$id %>%
    gsub(">", "", .) %>%
    strsplit(., split = "-") %>%
    lapply(., `[[`, 1) %>%
    unlist() %>%
    as.numeric()
  countData$Reads <- countData$id %>%
    strsplit(., split = "-") %>%
    lapply(., `[[`, 2) %>%
    unlist() %>%
    as.numeric()
  countData$RPM <- countData$id %>%
    strsplit(., split = "-") %>%
    lapply(., `[[`, 3) %>%
    unlist() %>%
    as.numeric()
  
  # handle clustered input
  if(lengths(strsplit(countData$id[1], split = "-")) == 6){
    countData$cluster <- countData$id %>%
      strsplit(., split = "-") %>%
      lapply(., `[[`, 4) %>%
      unlist() %>%
      as.numeric()
    countData$rankInCluster <- countData$id %>%
      strsplit(., split = "-") %>%
      lapply(., `[[`, 5) %>%
      unlist() %>%
      as.numeric()
    countData$LED <- countData$id %>%
      strsplit(., split = "-") %>%
      lapply(., `[[`, 6) %>%
      unlist() %>%
      as.numeric()
    
    return(return(countData[,c(1,3:8,2)]))
  } else{
    return(countData[,c(1,3:5,2)])
  }
}
fa_motif_motifTracker <- function(fastaInputs = NULL, fileNames = NULL, queryList = NULL, queryAliases = NULL, motifType = "Nucleotide"){
  
  # one motif query per line
  queryList <- strsplit(queryList, "\n") %>% unlist()
  
  # one alias per line
  if(!is.null(queryAliases)){
    queryAliases <- strsplit(queryAliases, "\n") %>% unlist()
    
    # confirm one alias per query
    if(length(queryList) != length(queryAliases)){
      message("Number of aliases should equal number of queries!")
      return(NULL)
    }
  }
  
  # for each query, format the motif
  queryList_mod <- queryList %>%
    sapply(., function(x) fa_motif_format(motifList = x, motifType = motifType)) %>%
    as.vector() %>%
    paste0("(?=.*", ., ")") %>%
    gsub("\\|", ")(?=.*", .)
  
  # initialize target data.frame
  targetDF <- data.frame()
  
  # initialize list of sequence counts
  lengthList <- list()
  
  # iterate through each fasta file
  for(i in 1:length(fastaInputs)){
    
    # read and format FASTA input
    formatDF <- fa_formatInput(fastaInput = fastaInputs[i])
    
    # get number of sequences in FASTA file
    lengthList[[i]] <- sum(formatDF$Reads) / 1e6
    
    # iterate through each formatted query
    for(j in 1:length(queryList_mod)){
      
      # row bind targetDF with formatted fasta input after searching for query
      targetDF <- rbind(
        targetDF,
        formatDF %>%
          dplyr::filter(., grepl(queryList_mod[j], seqs, perl = T)) %>%
          dplyr::mutate(., Motif = queryList[j], Population = i, FileName = fileNames[i])
      )
    }
  }
  
  # unlist length list
  lengthList <- lengthList %>% unlist()
  
  # group by population and motif, get total reads and RPM for each motif
  targetDF <- targetDF %>%
    dplyr::group_by(., Population, FileName, Motif) %>%
    dplyr::summarise(., TotalReads = sum(Reads)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(., TotalRPM = round(TotalReads / lengthList[Population], 2)) %>%
    dplyr::arrange(Motif)
  
  # add Alias column if they were supplied
  if(!is.null(queryAliases)){
    # initialize column
    targetDF <- targetDF %>% dplyr::mutate(Alias = NA, .after = Motif)
    
    # iterate through IDs
    for(i in 1:length(queryAliases)){
      targetDF <- targetDF %>%
        dplyr::mutate(Alias = ifelse(Motif == queryList[i], queryAliases[i], Alias))
    }
  }
  
  # return targetDF
  return(targetDF)
}
fa_motif_sequenceTracker <- function(fastaInputs = NULL, fileNames = NULL, queryList = NULL, queryAliases = NULL){
  
  # one sequence query per line
  queryList <- strsplit(queryList, "\n") %>% unlist()
  
  # one alias per line
  if(!is.null(queryAliases)){
    queryAliases <- strsplit(queryAliases, "\n") %>% unlist()
    
    # confirm one alias per query
    if(length(queryList) != length(queryAliases)){
      message("Number of aliases should equal number of queries!")
      return(NULL)
    }
  }
  
  # initialize targetDF
  targetDF <- data.frame()
  
  # iterate through all input files; format file, pull sequences from query list, add population; bind to previous df
  for(i in 1:length(fastaInputs)){
    targetDF <- rbind(
      targetDF,
      fa_formatInput(fastaInput = fastaInputs[i]) %>%
        dplyr::filter(., seqs %in% queryList) %>%
        dplyr::mutate(., Population = i, FileName = fileNames[i])
    )
  }
  
  # rearrange columns and sort by seqs
  targetDF <- targetDF %>%
    dplyr::relocate(., c(Population, FileName, seqs, Rank, Reads, RPM)) %>%
    dplyr::arrange(seqs)
  
  # add Alias column if they were supplied
  if(!is.null(queryAliases)){
    # initialize column
    targetDF <- targetDF %>% dplyr::mutate(Alias = NA, .after = seqs)
    
    # iterate through IDs
    for(i in 1:length(queryAliases)){
      targetDF <- targetDF %>%
        dplyr::mutate(Alias = ifelse(seqs == queryList[i], queryAliases[i], Alias))
    }
  }
  
  # return targetDF
  return(targetDF)
}
fa_distance <- function(dataInput = NULL, querySequence = NULL, seqRange = c(1, 1000)){
  # read file
  if(grepl("CSV", toupper(dataInput))){
    inputDF <- read.csv(dataInput)
  } else{
    inputDF <- fa_formatInput(fastaInput = dataInput)
  }
  
  # add ID field
  inputDF$id <- paste0(">", inputDF$Rank, "-", inputDF$Reads, "-", inputDF$RPM)
  
  # most inputs have a seqs column, but cluster analysis gives a seeds column; standardize the naming
  if("Seeds" %in% colnames(inputDF)){
    inputDF <- inputDF %>%
      dplyr::rename(., seqs = Seeds)
  }
  
  # get length of query sequence
  seqLength <- nchar(querySequence)
  
  # truncate sequences if 1) min(seqRange) > 1 and/or 2) max(seqRange) < seqLength
  querySequence_trunc <- substr(querySequence,
                                start = min(seqRange),
                                stop = min(seqLength, max(seqRange)))
  
  # check if the truncated sequence is the same as the input sequence
  if(querySequence == querySequence_trunc){
    # compute LED between each target sequence and the query sequence
    inputDF$Distance <- drop(adist(x = toupper(inputDF$seqs), y = toupper(querySequence)))
    
    # rearrange columns
    inputDF <- inputDF[,c(5, 2:4, 6, 1)]
  } else{
    # truncate target sequences
    inputDF$TruncSeqs <- substr(inputDF$seqs,
                                start = min(seqRange),
                                stop = min(seqLength, max(seqRange)))
    
    # compute LED between each truncated target sequence and the truncated query sequence
    inputDF$Distance <- drop(adist(x = toupper(inputDF$TruncSeqs), y = toupper(querySequence_trunc)))
    
    # rearrange columns
    inputDF <- inputDF[,c(5, 2:4, 7, 1, 6)]
  }
  
  # order by distance and rearrange columns
  inputDF <- inputDF[order(as.numeric(inputDF$Distance)),]
  
  # return df after ordering by distance
  return(inputDF)
}
fa_enrich <- function(fastaInputs = NULL, keepNA = FALSE){
  # initialize list of formatted files
  inputList <- list()
  
  # read and format all input files; store them as a list
  for(i in 1:length(fastaInputs)){
    inputList[[i]] <- fa_formatInput(fastaInput = fastaInputs[i], population = letters[i])
  }
  
  # merge all files together
  for(i in 2:length(fastaInputs)){
    # this 1st merge initializes the df
    if(i == 2){
      mergeDF <- merge(inputList[[i-1]], inputList[[i]], by = "seqs", all = ifelse(keepNA, TRUE, FALSE))
    } else{
      # merge every other file into this merged df
      mergeDF <- merge(mergeDF, inputList[[i]], by = "seqs", all = ifelse(keepNA, TRUE, FALSE))
    }
  }
  
  # calculate enrichment and log2(enrichment) between all consecutive files
  for(i in 2:length(fastaInputs)){
    # enrichment
    mergeDF[[paste0("enrichment_", letters[i], letters[i-1])]] <- mergeDF[[paste0("RPM.", letters[i])]] / mergeDF[[paste0("RPM.", letters[i-1])]]
    
    # log2(enrichment)
    mergeDF[[paste0("log2E_", letters[i], letters[i-1])]] <- log2(mergeDF[[paste0("enrichment_", letters[i], letters[i-1])]])
    
    # finally, round both values to 3 decimal points
    mergeDF[[paste0("enrichment_", letters[i], letters[i-1])]] <- round(mergeDF[[paste0("enrichment_", letters[i], letters[i-1])]], 3)
    mergeDF[[paste0("log2E_", letters[i], letters[i-1])]] <- round(mergeDF[[paste0("log2E_", letters[i], letters[i-1])]], 3)
  }
  
  # replace NAs with 0
  mergeDF <- mergeDF %>% replace(is.na(.), 0)
  
  # return merged data.frame after ordering by 'Rank' in the first population
  return(mergeDF[order(as.numeric(mergeDF$Rank.a)),])
}
fa_clusterLED <- function(fastaInput = NULL, minReads = 10, maxLED = 7, totalClusters = 30, multipleOutputs = F, outputDirectory = NULL, keepNC = T){
  # format the output directory path for R
  if(outputDirectory != "" & dir.exists(outputDirectory)){
    outputDirectory <- outputDirectory
  }
  
  # read FASTA file and convert to data.frame; only keep sequences in which 'Reads' is greater than minReads
  countData <- fa_formatInput(fastaInput = fastaInput)
  countData <- countData[which(countData$Reads > minReads),]
  
  # initialize extra columns in data.frame
  countData$cluster <- NA
  countData$rankInCluster <- NA
  countData$LED <- NA
  
  # sample sequences; assume they are sorted by number of reads
  seqs <- as.list(countData$seqs)
  
  # create cluster iterator
  clusterNumber <- 1
  
  # only create a user-specified number of clusters and while you still have sequences
  while(clusterNumber <= totalClusters & length(seqs) != 0){
    
    # get start time
    startTime <- Sys.time()
    
    # reset iterator (counter)
    counter <- 1
    
    # first element of seq. list (most abundant seq.) is seed; remove from seq. list
    clusterList <- seqs[1]
    seqs <- seqs[-1]
    
    # make list for string distances
    seqLED <- list(0)
    
    # only iterate while you still have sequences and while you are not at the end of the list
    while(length(seqs) != 0 & counter < length(seqs)){
      
      # check LED between seed and current sequence
      if(as.numeric(adist(clusterList[[1]], seqs[[counter]])) <= maxLED){
        
        # add LED to seqLED list
        seqLED <- c(seqLED, as.numeric(adist(clusterList[[1]], seqs[[counter]])))
        
        # if LED is < maxED, add sequence to cluster list; set sequence from original list to NA
        clusterList <- c(clusterList, seqs[[counter]])
        seqs[[counter]] <- NA
      }
      
      # check next sequence in list
      counter <- counter + 1
    }
    
    # remove NAs from sequence list; these seq.s are now clustered
    seqs <- seqs[!is.na(seqs)]
    
    # turn clustered sequences in data frame and merge back into countData
    seqDF <- data.frame(seqs = unlist(clusterList), cluster = clusterNumber,
                        rankInCluster = 1:length(clusterList), LED = unlist(seqLED),
                        stringsAsFactors = F)
    countData[which(countData$seqs %in% seqDF$seqs),5:7] <- seqDF[,2:4]
    
    # check if user wants multiple outputs
    if(multipleOutputs){
      clustDF <- countData[which(countData$cluster == clusterNumber),]
      clustDF$id <- paste0('>', clustDF[,2], '-', clustDF[,3], '-', clustDF[,4], '-',
                           clustDF[,5], '-', clustDF[,6], '-', clustDF[,7])
      write.table(fa_formatOutput(outputData = clustDF[,c(1,8)]), file = paste0(outputDirectory, clusterNumber, ".fasta"),
                  quote = F, row.names = F, col.names = F)
    }
    
    # message with number of unique sequences in cluster
    message(paste(paste0("Finished cluster ", clusterNumber, ": ", nrow(seqDF), " unique sequences"),
                  paste0(gsub("Time difference of", "Elapsed time:", capture.output(round(Sys.time() - startTime, 2)))),
                  sep = "<br/>"))

    # bump the cluster number
    clusterNumber <- clusterNumber + 1
  }
  
  # optionally replace the NA values of non-clustered sequences with 'NC'
  if(keepNC){
    countData[is.na(countData)] <- 'NC'
  } else{
    countData <- na.omit(countData)
  }
  
  # last output if user desires multiple outputs
  if(multipleOutputs){
    if(keepNC){
      clustDF <- countData[which(countData$cluster == "NC"),]
      clustDF$id <- paste0('>', clustDF[,2], '-', clustDF[,3], '-', clustDF[,4], '-',
                           clustDF[,5], '-', clustDF[,6], '-', clustDF[,7])
      write.table(fa_formatOutput(outputData = clustDF[,c(1,8)]), file = paste0(outputDirectory, "NC.fasta"),
                  quote = F, row.names = F, col.names = F)
    }
  }else{
    # make id
    countData$id <- paste0('>', countData[,2], '-', countData[,3], '-', countData[,4], '-',
                           countData[,5], '-', countData[,6], '-', countData[,7])
    
    return(countData[,c(8,2:7,1)])
  }
}
fa_clusterDiversity <- function(clusterFASTA = NULL){
  # read and format output from cluster function
  countData <- fa_formatInput(fastaInput = clusterFASTA, population = NULL)
  countData[countData == 'NC'] <- NA
  
  clDF_format <- countData[which(countData$ClusterRank == 1), c(5,1)]
  colnames(clDF_format)[2] <- 'Seeds'
  
  # get additional features per cluster
  clusterStats <- countData %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(TotalSequences = dplyr::n(), TotalReads = sum(Reads), TotalRPM = sum(RPM),
                     AverageLED = round(mean(as.numeric(as.character(SeedDistance))),2)) %>%
    as.data.frame(.)

  clusterStats <- clusterStats[order(as.numeric(clusterStats$Cluster)),]
  clusterStats[is.na(clusterStats)] <- 'NC'
  
  # return data.frame with clustering metadata
  return(merge(clDF_format, clusterStats, by = 'Cluster', all = T, sort = F))
}
fa_clusterEnrich <- function(clusterCSVs = NULL, fileNames = NULL){
  # initialize list of formatted files
  stackedDF <- data.frame()
  
  # read and format all input files; store them as a list
  for(i in 1:length(clusterCSVs)){
    stackedDF <- rbind(
      stackedDF,
      read.csv(clusterCSVs[i]) %>% dplyr::mutate(., Population = i, FileName = fileNames[i])
    )
  }
  
  # return data.frame after merging by cluster number in first population
  return(stackedDF)
}

##### UI functions #####
fa_count_outputName <- function(inputFile = NULL, outputType = NULL){
  # if output extension is 'FASTA', replace previous extension with '-count.fasta'
  # else, replace previous extension with '-count.csv'
  if(outputType == "FASTA"){
    fileName <- ifelse(sum(endsWith(toupper(inputFile), c("FASTA", "FASTQ"))) > 0,
                       sub(rightSubstr(inputFile, 6), "-count.fasta", inputFile),
                       sub(rightSubstr(inputFile, 3), "-count.fasta", inputFile))
  }else{
    fileName <- ifelse(sum(endsWith(toupper(inputFile), c("FASTA", "FASTQ"))) > 0,
                       sub(rightSubstr(inputFile, 6), "-count.csv", inputFile),
                       sub(rightSubstr(inputFile, 3), "-count.csv", inputFile))
  }
  
  # return output file name
  return(fileName)
}
fa_count_metadata <- function(countData = NULL){
  # get number of (non-)redundant sequences
  uniqueSeqs <- paste0("Unique sequences: ", nrow(countData))
  totalSeqs <- paste0("Total sequences: ", countData$id %>%
                        strsplit(., "-") %>% lapply(., `[[`, 2) %>%
                        unlist(.) %>% as.numeric(.) %>% sum(.))
  
  # return number of (non-)redundant sequences
  return(c(totalSeqs, uniqueSeqs))
}
fa_motifTracker_scores <- function(RPM = NULL){
  # compute enrichment scores
  enrichment <- list()
  for(i in 2:length(RPM)){
    for(j in i:1){
      if(i != j){
        enrichment <- append(enrichment, paste0("Enrichment (populations ", i, ":", j, "): ", round(RPM[i] / RPM[j], 3)))
      }
    }
  }
  
  # return list of enrichment scores
  return(unlist(enrichment))
}
rightSubstr <- function(x, n){
  # return the last n characters of string
  return(substr(x, nchar(x) - n + 1, nchar(x)))
}
patternPerms <- function(patterns = NULL){
  if(length(patterns) == 1){
    return(patterns)
  } else{
    patternMatrix <- matrix(nrow = 0, ncol = length(patterns))
    for(i in seq_along(patterns)){
      patternMatrix <- rbind(patternMatrix, cbind(patterns[i], Recall(patterns[-i])))
    }
    return(patternMatrix)
  }
}