# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function takes a set of input FASTA files and merges them according to user-defined criteria.
#' Currently, the only supported merges are full, inner, and left.
fa_dataMerge <- function(fastaInputs = NULL, mergeType = "Intersection"){
  
  # initialize empty list
  sequenceDF_list <- list()
  
  # add each DF-converted FASTA file to the list
  for(i in 1:length(fastaInputs)){
    sequenceDF_list[[i]] <- fa_formatInput(fastaInput = fastaInputs[i], population = letters[i])
  }
  
  # do the merging
  if(mergeType == "Union"){
    
    mergeDF <- plyr::join_all(sequenceDF_list, by = "seqs", type = "full")
  }
  
  if(mergeType == "Intersection"){
    
    mergeDF <- plyr::join_all(sequenceDF_list, by = "seqs", type = "inner")
  }
  
  if(mergeType == "Left"){
    
    mergeDF <- plyr::join_all(sequenceDF_list, by = "seqs", type = "left")
  }
  
  # return the merged data.frame
  return(mergeDF)
}
