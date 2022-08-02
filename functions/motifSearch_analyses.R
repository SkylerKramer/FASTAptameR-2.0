# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function filters sequences for the presence of user-defined motif(s)
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
