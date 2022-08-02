# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function computes the enrichment of sequences across multiple rounds
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
  mergeDF <- mergeDF %>% replace(is.na(.), 0.01)
  
  # return merged data.frame after ordering by 'Rank' in the first population
  return(mergeDF[order(as.numeric(mergeDF$Rank.a)),])
}
