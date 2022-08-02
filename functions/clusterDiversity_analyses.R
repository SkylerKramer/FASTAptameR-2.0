# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function returns summary statistics for a clustered FASTA file
fa_clusterDiversity <- function(clusterFASTA = NULL){
  # read and format output from cluster function
  countData <- fa_formatInput(fastaInput = clusterFASTA, population = NULL)
  countData[countData == 'NC'] <- NA
  
  clDF_format <- countData[which(countData$ClusterRank == 1), c(5,1)]
  colnames(clDF_format)[2] <- 'Seeds'
  
  # get additional features per cluster
  clusterStats <- countData %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      TotalSequences = dplyr::n(),
      TotalReads = sum(Reads),
      TotalRPM = sum(RPM),
      AverageLED = round(mean(as.numeric(as.character(SeedDistance))),2)
    ) %>%
    as.data.frame(.)
  
  clusterStats <- clusterStats[order(as.numeric(clusterStats$Cluster)),]
  clusterStats[is.na(clusterStats)] <- 'NC'
  
  # return data.frame with clustering metadata
  return(merge(clDF_format, clusterStats, by = 'Cluster', all = TRUE, sort = FALSE))
}
