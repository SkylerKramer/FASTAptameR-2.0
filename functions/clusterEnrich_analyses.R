# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function looks at the enrichment of clusters across multiple rounds
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

#' This function computes the enrichment of clusters across multiple rounds
fa_clusterEnrich_tracker <- function(trackerDF = NULL){
  
  # rename 2nd column (either seqs or Motif) to Query
  colnames(trackerDF)[2] <- "Query"
  
  # ungroup data.frame, omit NAs, sort by Query, and omit any queries that are not duplicated
  trackerDF <- trackerDF %>%
    dplyr::ungroup() %>%
    na.omit() %>%
    dplyr::arrange(Query) %>%
    subset(., duplicated(Query) | duplicated(Query, fromLast = TRUE))
  
  # initialize empty data.frame to hold enrichment scores
  enrichmentDF <- data.frame()
  
  # fill enrichmentDF with one query at a time
  for(query in unique(trackerDF$Query)){
    
    # split tibble by query
    temp <- trackerDF[which(trackerDF$Query == query),] %>% as.data.frame()
    
    # iterate through each consecutive population pair
    for(i in 2:nrow(temp)){
      
      # fill enrichmentDF with comparisons and corresponding enrichments; column 4 is TotalRPM
      enrichmentDF <- rbind(
        enrichmentDF, data.frame(
          Comparison = paste0(temp$FileName[i], " : ", temp$FileName[i-1]),
          Clusters = paste0(temp$Cluster[i], " : ", temp$Cluster[i-1]),
          Query = query,
          Enrichment = round(temp$TotalRPM[i] / temp$TotalRPM[i-1], 3)
        )
      )
    }
  }
  
  # return DF
  return(enrichmentDF)
}
