### Skyler Kramer - FASTAptameR 2.0

##### imports ####
library(dplyr)

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
fa_enrich_merge <- function(countData1 = NULL, countData2 = NULL, countData3 = NULL, removeNA = T){
  # merge at least two counted files - and optionally a third counted file - by 'seqs'
  countMerge <- merge(countData1, countData2, by = 'seqs', all = T)
  if(!is.null(countData3)){
    countMerge <- merge(countMerge, countData3, by = 'seqs', all = T)
  }
  
  # optionally keep na values
  if(removeNA){
    countMerge <- na.omit(countMerge)
  }

  # enrichment and log2(E) (the fold change): y to x
  countMerge$enrichment_yx <- round(countMerge$RPM.y / countMerge$RPM.x, 3)
  countMerge$log2_E_yx <- round(log2(countMerge$RPM.y / countMerge$RPM.x), 3)
  
  #  enrichment and log2(E) (the fold change): z to y, z to x
  if(!is.null(countData3)){
    countMerge$enrichment_zy <- round(countMerge$RPM.z / countMerge$RPM.y, 3)
    countMerge$log2_E_zy <- round(log2(countMerge$RPM.z / countMerge$RPM.y), 3)
    
    countMerge$enrichment_zx <- round(countMerge$RPM.z / countMerge$RPM.x, 3)
    countMerge$log2_E_zx <- round(log2(countMerge$RPM.z / countMerge$RPM.x), 3)
  }
  
  # return merged file with fold changes and enrichment scores
  return(countMerge)
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
  p1 <- ggplot2::ggplot(countData, ggplot2::aes(Length)) + ggplot2::geom_bar(fill = "blue") +
    ggplot2::xlab("Sequence length") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
  p2 <- countData %>%
    dplyr::group_by(Length) %>%
    dplyr::summarise(TotalReads = sum(Reads)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Length, y = TotalReads)) +
    ggplot2::geom_bar(stat = "identity", fill = "blue") +
    ggplot2::xlab("Sequence length") + ggplot2::ylab("Total number of reads") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
  # return interactive plot
  return(plotly::subplot(plotly::ggplotly(p1), plotly::ggplotly(p2), titleY = T, titleX = T, nrows = 2, margin = 0.1) %>%
    plotly::layout(title = "Sequence-Length Histogram", showlegend = F, height = 600, margin = list(t = 50)))
}
fa_distance_histogram <- function(distanceData = NULL, querySequence = NULL){
  # distance histogram with unique sequences
  p1 <- ggplot2::ggplot(distanceData, ggplot2::aes(Distance)) + ggplot2::geom_bar(fill = "blue") +
    ggplot2::xlab("Distance") + ggplot2::ylab("Number of unique sequences") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
  p1 <- p1 %>%
    plotly::ggplotly() %>%
    plotly::add_text(x = mean(distanceData$Distance), y = max(ggplot2::ggplot_build(p1)$data[[1]]$ymax),
                     text = querySequence, textfont = list(size = 12))
  
  # distance histogram with total reads
  p2 <- distanceData %>%
    dplyr::group_by(Distance) %>%
    dplyr::summarise(TotalReads = sum(Reads)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Distance, y = TotalReads)) +
    ggplot2::geom_bar(stat = "identity", fill = "blue") +
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
fa_enrich_scatter <- function(df = NULL){
  # if only two populations (two occurrences of 'RPM' substring), create interactive 2d scatter plot
  # else, create interactive 3d scatter plot
  if(sum(grepl("RPM", colnames(df))) == 2){
    plotly::plot_ly(df, type = "scatter", mode = "markers", x = ~RPM.x, y = ~RPM.y, text = ~seqs) %>%
      plotly::layout(title = "RPM.x vs RPM.y",
                     xaxis = list(type = "log"), yaxis = list(type = "log"))
  }else{
    plotly::plot_ly(df, type = "scatter3d", mode = "markers", x = ~RPM.x, y = ~RPM.y, z = ~RPM.z, text = ~seqs) %>%
      plotly::layout(title = "RPM.x vs RPM.y vs RPM.z",
                     xaxis = list(type = "log"), yaxis = list(type = "log", zaxis = list(type = "log")))
  }
}
fa_enrich_histogram <- function(df = NULL){
  # if only two populations (one occurence of the 'foldChange' substring), create a single interactive histogram
  # else, create one interactive histogram per comparison (three histograms total): y-x, z-y, z-x
  if(sum(grepl("log2_E", colnames(df))) == 1){
    fig <- plotly::ggplotly(ggplot2::ggplot(df, ggplot2::aes(log2_E_yx)) +
                              ggplot2::geom_histogram(color = "black", fill = "blue", bins = 30) + 
                              ggplot2::ggtitle("Population y vs Population x") + ggplot2::xlab("log2(Enrichment)") +
                              ggplot2::ylab("Number of Unique Sequences") +
                              ggplot2::theme_bw() + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45)))
  } else{
    fig1 <- plotly::ggplotly(ggplot2::ggplot(df, ggplot2::aes(log2_E_yx)) +
                               ggplot2::geom_histogram(color = "black", fill = "blue", bins = 30) + 
                               ggplot2::ggtitle("Population y vs Population x") + ggplot2::xlab("log2(Enrichment)") +
                               ggplot2::ylab("Number of Unique Sequences") +
                               ggplot2::theme_bw() + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45)))
    
    fig2 <- plotly::ggplotly(ggplot2::ggplot(df, ggplot2::aes(log2_E_zy)) +
                               ggplot2::geom_histogram(color = "black", fill = "orange", bins = 30) + 
                               ggplot2::ggtitle("Population z vs Population y") + ggplot2::xlab("log2(Enrichment)") +
                               ggplot2::ylab("Number of Unique Sequences") +
                               ggplot2::theme_bw() + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45)))
    
    fig3 <- plotly::ggplotly(ggplot2::ggplot(df, ggplot2::aes(log2_E_zx)) +
                               ggplot2::geom_histogram(color = "black", fill = "green", bins = 30) + 
                               ggplot2::ggtitle("Population z vs Population x") + ggplot2::xlab("log2(Enrichment)") +
                               ggplot2::ylab("Number of Unique Sequences") +
                               ggplot2::theme_bw() + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45)))
    
    fig <- plotly::subplot(fig1, fig2, fig3) %>%
      plotly::layout(title = "log2(Enrichment) Histogram", showlegend = F)
  }
  
  # return histogram(s)
  return(fig)
}
fa_enrich_volcano <- function(df = NULL){
  # if only two populations (two occurrences of the 'RPM' substring), create a single interactive volcano plot
  # else, create one interactive volcano plot per comparison (three volcano plots total): y-x, z-y, z-x
  if(sum(grepl("RPM", colnames(df))) == 2){
    fig <- plotly::plot_ly(df, type = "scatter", mode = "markers",
                   x = ~log2_E_yx, y = ~sqrt(log2(Reads.y) / log2(sum(Reads.y))), text = ~seqs) %>%
      plotly::layout(title = "RPM Volcano Plot",
             xaxis = list(title = "log2(Enrichment)"),
             yaxis = list(title = "Sqrt(log2(Reads.y) / log2(sum(Reads.y)))"))
  } else{
    fig1 <- plotly::plot_ly(df, type = "scatter", mode = "markers",
                    x = ~log2_E_yx, y = ~sqrt(log2(Reads.y) / log2(sum(Reads.y))), text = ~seqs) %>%
      plotly::layout(xaxis = list(title = "log2(Enrichment) - (y:x)"),
             yaxis = list(title = "Sqrt(log2(Reads.y) / log2(sum(Reads.y)))"))
    
    fig2 <- plotly::plot_ly(df, type = "scatter", mode = "markers",
                    x = ~log2_E_zy, y = ~sqrt(log2(Reads.z) / log2(sum(Reads.z))), text = ~seqs) %>%
      plotly::layout(xaxis = list(title = "log2(Enrichment) - (z:y)"),
             yaxis = list(title = "Sqrt(log2(Reads.z) / log2(sum(Reads.z)))"))
    
    fig3 <- plotly::plot_ly(df, type = "scatter", mode = "markers",
                    x = ~log2_E_zx, y = ~sqrt(log2(Reads.z) / log2(sum(Reads.z))), text = ~seqs) %>%
      plotly::layout(xaxis = list(title = "log2(Enrichment) - (z:x)"),
             yaxis = list(title = "Sqrt(log2(Reads.z) / log2(sum(Reads.z)))"))
    
    fig <- plotly::subplot(fig1, fig2, fig3, titleY = T, titleX = T) %>%
      plotly::layout(title = "log2(Enrichment) Volcano Plot", showlegend = F)
  }
  
  # return volcano plot(s)
  return(fig)
}
fa_clusterDiversity_kmerPCA <- function(clusterFile = NULL, kmerSize = 3, topClusters = 10, keepNC = T){
  # format clustered file
  clusterDF <- fa_formatInput(fastaInput = clusterFile, population = NULL)
  
  # return NULL if the sequences are made of a non-nucleotide alphabet
  if(sum(grepl("[^ACGTU]", clusterDF$seqs)) != 0){
    return(NULL)
  }
  
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
fa_translate <- function(fastaInput = NULL, orf = 1, converge = T){
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
  
  # translate sequences, remove white spaces
  seqs <- gsub("UUU|UUC", "F", seqs) %>%
    gsub("UUA|UUG|CUU|CUC|CUA|CUG", "L", .) %>%
    gsub("AUU|AUC|AUA", "I", .) %>%
    gsub("AUG", "M", .) %>%
    gsub("GUU|GUC|GUA|GUG", "V", .) %>%
    gsub("UCU|UCC|UCA|UCG|AGU|AGC", "S", .) %>%
    gsub("CCU|CCC|CCA|CCG", "P", .) %>%
    gsub("ACU|ACC|ACA|ACG", "T", .) %>%
    gsub("GCU|GCC|GCA|GCG", "A", .) %>%
    gsub("UAU|UAC", "Y", .) %>%
    gsub("UAA|UAG|UGA", "-", .) %>%
    gsub("CAU|CAC", "H", .) %>%
    gsub("CAA|CAG", "Q", .) %>%
    gsub("AAU|AAC", "N", .) %>%
    gsub("AAA|AAG", "K", .) %>%
    gsub("GAU|GAC", "D", .) %>%
    gsub("GAA|GAG", "E", .) %>%
    gsub("UGU|UGC", "C", .) %>%
    gsub("UGG", "W", .) %>%
    gsub("CGU|CGC|CGA|CGG|AGA|AGG", "R", .) %>%
    gsub("GGU|GGC|GGA|GGG", "G", .) %>%
    gsub(" ", "", .)
  
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
    translateDF <- translateDF[order(translateDF$Rank),]
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
fa_motifTracker <- function(fastaInputs = NULL, motif = NULL, fileNames = NULL, motifType = "Nucleotide"){
  # read files and get total reads per population
  countData1 <- fa_formatInput(fastaInput = fastaInputs[1], population = 'x')
  totalReads.x <- sum(countData1$Reads.x)
  
  countData2 <- fa_formatInput(fastaInput = fastaInputs[2], population = 'y')
  totalReads.y <- sum(countData2$Reads.y)
  
  if(length(fastaInputs) == 3){
    countData3 <- fa_formatInput(fastaInput = fastaInputs[3], population = 'z')
    totalReads.z <- sum(countData3$Reads.z)
  }
  
  # format motif, convert Boolean OR ("|") into character vector, make into matrix of all pattern permutations
  motif <- fa_motif_format(motifList = motif, motifType = "Nucleotide") %>%
    strsplit(., "\\|") %>% unlist(.) %>%
    patternPerms(patterns = .) %>% as.matrix(.)
  
  # make one regex pattern from all permutations
  motif[,1] <- paste0("(?=", motif[,1])
  motif[,ncol(motif)] <- paste0(motif[,ncol(motif)], ")")
  motif <- apply(motif, 1, paste0, collapse = ".*") %>%
    paste0(., collapse = "|")
  
  # count number of pattern occurrences
  motifCount.x <- lengths(regmatches(countData1$seqs, gregexpr(motif, countData1$seqs, perl = T)))
  motifCount.y <- lengths(regmatches(countData2$seqs, gregexpr(motif, countData2$seqs, perl = T)))
  
  if(length(fastaInputs) == 3){
    motifCount.z <- lengths(regmatches(countData3$seqs, gregexpr(motif, countData3$seqs)))
  }
  
  # format output
  motifEnrich <- data.frame(File.Names = fileNames[1:2],
                            Unique.Reads = c(nrow(countData1),
                                             nrow(countData2)),
                            Total.Reads = c(totalReads.x,
                                            totalReads.y),
                            Seqs.With.Motif = c(sum(countData1[motifCount.x != 0, "Reads.x"]),
                                                sum(countData2[motifCount.y != 0, "Reads.y"])),
                            Seqs.RPM = c(round(sum(countData1[motifCount.x != 0, "Reads.x"]) / (totalReads.x / 1e6), 2),
                                         round(sum(countData2[motifCount.y != 0, "Reads.y"]) / (totalReads.y / 1e6), 2)),
                            Motif.Occurrences = c(sum(motifCount.x * countData1$Reads.x),
                                                  sum(motifCount.y * countData2$Reads.y)),
                            Motif.RPM = c(round(sum(motifCount.x * countData1$Reads.x) / (totalReads.x / 1e6), 2),
                                          round(sum(motifCount.y * countData2$Reads.y) / (totalReads.y / 1e6), 2)))
  
  if(length(fastaInputs) == 3){
    motifEnrich[3,] <- c(fileNames[3],
                         nrow(countData3),
                         totalReads.z,
                         sum(countData3[motifCount.z != 0, "Reads.z"]),
                         round(sum(countData1[motifCount.z != 0, "Reads.z"]) / (totalReads.z / 1e6), 2),
                         sum(motifCount.z * countData3$Reads.z),
                         round(sum(motifCount.z * countData3$Reads.z) / (totalReads.z / 1e6), 2))
  }
  
  # make RPM column a numeric
  motifEnrich$Motif.RPM <- as.numeric(motifEnrich$Motif.RPM)  
  
  # message with enrichment scores
  motifEnrich_messages <- rep(NA, ifelse(length(fastaInputs) == 3, 3, 1))
  motifEnrich_messages[1] <- paste0('Enrichment (population 2 to population 1): ', round(motifEnrich$Motif.RPM[2] / motifEnrich$Motif.RPM[1], 3))
  
  if(length(fastaInputs) == 3){
    motifEnrich_messages[2] <- paste0('Enrichment (population 3 to population 2): ', round(motifEnrich$Motif.RPM[3] / motifEnrich$Motif.RPM[2], 3))
    motifEnrich_messages[3] <- paste0('Enrichment (population 3 to population 1): ', round(motifEnrich$Motif.RPM[3] / motifEnrich$Motif.RPM[1], 3))
  }
  
  # return data.frame
  return(motifEnrich)
}
fa_distance <- function(dataInput = NULL, querySequence = NULL){
  # read file
  if(grepl("CSV", toupper(dataInput))){
    inputDF <- read.csv(dataInput)
  } else{
    inputDF <- fa_formatInput(fastaInput = dataInput)
  }
  
  # compute LED between each sequence and the query sequence
  inputDF$Distance <- drop(adist(x = toupper(inputDF$seqs), y = toupper(querySequence)))
  
  # return df after ordering by distance
  return(inputDF[order(as.numeric(inputDF$Distance)),])
}
fa_enrich <- function(fastaInputs = NULL, removeNA = T){
  # read and process FASTA files; merge by sequences and optionally remove sequences missing in multiple populations
  countData1 <- fa_formatInput(fastaInput = fastaInputs[1], population = 'x')
  countData2 <- fa_formatInput(fastaInput = fastaInputs[2], population = 'y')
  
  if(length(fastaInputs) == 3){
    countData3 <- fa_formatInput(fastaInput = fastaInputs[3], population = 'z')
    countMerge <- fa_enrich_merge(countData1 = countData1, countData2 = countData2, countData3 = countData3, removeNA = removeNA)
  }else{
    countMerge <- fa_enrich_merge(countData1 = countData1, countData2 = countData2, removeNA = removeNA)
  }
  
  # return merged data.frame after ordering by 'Rank' in the first population
  return(countMerge[order(as.numeric(countMerge$Rank.x)),])
}
fa_clusterLED <- function(fastaInput = NULL, minReads = 10, maxLED = 7, totalClusters = 30, multipleOutputs = F, outputDirectory = NULL,
                          keepNC = T){
  # parameter check
  if(multipleOutputs == T & is.null(outputDirectory)){
    stop("Must provide directory to receive multiple outputs!")
  }
  if(multipleOutputs == T & !dir.exists(outputDirectory)){
    stop("Must provide valid directory path to receive multiple outputs!")
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
fa_clusterEnrich <- function(clusterCSVs = NULL){
  # read clustered CSVs
  cluster1 <- read.csv(clusterCSVs[1])
  cluster2 <- read.csv(clusterCSVs[2])
  
  if(length(clusterCSVs) == 3){
    cluster3 <- read.csv(clusterCSVs[3]) %>%
      rename(Cluster.z = Cluster, TotalSequences.z = TotalSequences, TotalReads.z = TotalReads, TotalRPM.z = TotalRPM)
  }
  
  # merge by seed sequences
  clusterMerge <- merge(cluster1, cluster2, by = 'Seeds', all = T)
  if(length(clusterCSVs) == 3){
    clusterMerge <- merge(clusterMerge, cluster3, by = 'Seeds', all = T)
  }
  clusterMerge <- na.omit(clusterMerge)
  
  
  # calculate enrichment
  clusterMerge$Enrichment_yx <- clusterMerge$TotalRPM.y / clusterMerge$TotalRPM.x
  if(length(clusterCSVs) == 3){
    clusterMerge$Enrichment_zy <- clusterMerge$TotalRPM.z / clusterMerge$TotalRPM.y
    clusterMerge$Enrichment_zx <- clusterMerge$TotalRPM.z / clusterMerge$TotalRPM.x
  }
  
  # return data.frame after merging by cluster number in first population
  return(clusterMerge[order(as.numeric(clusterMerge$Cluster.x)),])
}
fa_clusterEnrich_byCluster <- function(clusterCSV1 = NULL, clusterCSV2 = NULL, clusterCSV3 = NULL, keepNA = F,
                                        csvOutput = NULL){
  ## reads CSVs
  cluster1 <- read.csv(clusterCSV1)
  cluster2 <- read.csv(clusterCSV2)
  if(!is.null(clusterCSV3)){
    cluster3 <- read.csv(clusterCSV3)
  }
  
  ## merge
  clusterMerge <- merge(cluster1, cluster2, by = 'Cluster', all = T)
  if(!is.null(clusterCSV3)){
    clusterMerge <- merge(clusterMerge, cluster3, by = 'Cluster', all = T)
  }
  if(!keepNA){
    clusterMerge <- na.omit(clusterMerge)
  }
  
  ## calculate enrichment
  clusterMerge$Enrichment_yx <- clusterMerge$TotalRPM.y / clusterMerge$TotalRPM.x
  if(!is.null(clusterCSV3)){
    clusterMerge$Enrichment_zy <- clusterMerge$TotalRPM.z / clusterMerge$TotalRPM.y
    clusterMerge$Enrichment_zx <- clusterMerge$TotalRPM.z / clusterMerge$TotalRPM.x
  }
  
  ## order by Cluster
  clusterMerge[clusterMerge == 'NC'] <- NA
  clusterMerge <- clusterMerge[order(as.numeric(clusterMerge$Cluster)),]
  clusterMerge[is.na(clusterMerge)] <- 'NC'
  
  ## calculate distances between seeds
  clusterMerge$SeedDistance_yx <- diag(adist(clusterMerge$Seeds.y, clusterMerge$Seeds.x))
  if(!is.null(clusterCSV3)){
    clusterMerge$SeedDistance_zy <- diag(adist(clusterMerge$Seeds.z, clusterMerge$Seeds.y))
  }
  
  ## write the output
  write.csv(clusterMerge, csvOutput, row.names = F, quote = F)
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
fa_distance_outputName <- function(inputFile = NULL){
  # add "-distance" to filename
  fileName <- ifelse(sum(endsWith(toupper(inputFile), c("FASTA"))) > 0,
                     sub(rightSubstr(inputFile, 6), "-distance.csv", inputFile),
                     sub(rightSubstr(inputFile, 4), "-distance.csv", inputFile))
  
  # return output file name
  return(fileName)
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