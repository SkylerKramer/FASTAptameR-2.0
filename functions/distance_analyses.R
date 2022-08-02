# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function computes the LED between a query sequence and every other sequence in the target file
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
  querySequence_trunc <- substr(
    querySequence,
    start = min(seqRange),
    stop = min(seqLength, max(seqRange))
  )
  
  # check if the truncated sequence is the same as the input sequence
  if(querySequence == querySequence_trunc){
    
    # compute LED between each target sequence and the query sequence
    inputDF$Distance <- drop(adist(x = toupper(inputDF$seqs), y = toupper(querySequence)))
    
    # rearrange columns
    inputDF <- inputDF[,c(5, 2:4, 6, 1)]
  } else{
    
    # truncate target sequences
    inputDF$TruncSeqs <- substr(
      inputDF$seqs,
      start = min(seqRange),
      stop = min(seqLength, max(seqRange))
    )
    
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