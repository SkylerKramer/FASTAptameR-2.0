# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function clusters sequences according to a set of user-defined parameters
fa_clusterLED <- function(fastaInput = NULL, minReads = 10, maxLED = 7, totalClusters = 30, multipleOutputs = FALSE, outputDirectory = NULL, keepNC = TRUE){
  
  # format the output directory path for R
  if(outputDirectory != "" & dir.exists(outputDirectory)){
    outputDirectory <- outputDirectory
  }
  
  # read FASTA file and convert to data.frame; only keep sequences in which 'Reads' is greater than minReads; initialize extra columns
  countData <- fa_formatInput(fastaInput = fastaInput) %>%
    dplyr::filter(Reads >= minReads) %>%
    dplyr::mutate(cluster = NA, rankInCluster = NA, LED = NA)
  
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
    seqDF <- data.frame(
      seqs = unlist(clusterList),
      cluster = clusterNumber,
      rankInCluster = 1:length(clusterList),
      LED = unlist(seqLED),
      stringsAsFactors = FALSE
    )
    countData[which(countData$seqs %in% seqDF$seqs),5:7] <- seqDF[,2:4]
    
    # check if user wants multiple outputs
    if(multipleOutputs){
      clustDF <- countData[which(countData$cluster == clusterNumber),]
      clustDF$id <- paste0('>', clustDF[,2], '-', clustDF[,3], '-', clustDF[,4], '-', clustDF[,5], '-', clustDF[,6], '-', clustDF[,7])
      write.table(
        fa_formatOutput(outputData = clustDF[,c(1,8)]),
        file = paste0(outputDirectory, clusterNumber, ".fasta"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }
    
    # message with number of unique sequences in cluster
    message(
      paste(
        paste0("Finished cluster ", clusterNumber, ": ", nrow(seqDF), " unique sequences"),
        paste0(gsub("Time difference of", "Elapsed time:", capture.output(round(Sys.time() - startTime, 2)))),
        sep = "<br/>"
      )
    )
    
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
      clustDF$id <- paste0('>', clustDF[,2], '-', clustDF[,3], '-', clustDF[,4], '-', clustDF[,5], '-', clustDF[,6], '-', clustDF[,7])
      write.table(
        fa_formatOutput(outputData = clustDF[,c(1,8)]),
        file = paste0(outputDirectory, "NC.fasta"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }
  }else{
    
    # make id
    countData$id <- paste0('>', countData[,2], '-', countData[,3], '-', countData[,4], '-', countData[,5], '-', countData[,6], '-', countData[,7])
    
    return(countData[,c(8,2:7,1)])
  }
}
