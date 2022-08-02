# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function tracks how motif counts change between populations
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
          dplyr::filter(., grepl(queryList_mod[j], seqs, perl = TRUE)) %>%
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

#' This function tracks how sequence counts change between populations
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

#' This function reports the query enrichment from motifTracker or sequenceTracker
fa_motif_trackerEnrichment <- function(trackerDF = NULL){
  # ungroup data.frame
  trackerDF <- trackerDF %>% dplyr::ungroup()
  
  # rename 3rd column (either seqs or Motif) to Query
  colnames(trackerDF)[3] <- "Query"
  
  # initialize empty data.frame to hold enrichment scores
  enrichmentDF <- data.frame()
  
  # fill enrichmentDF with one query at a time
  for(query in unique(trackerDF$Query)){
    
    # split tibble by query
    temp <- trackerDF[which(trackerDF$Query == query),] %>% as.data.frame()
    
    # iterate through each consecutive population pair
    for(i in 2:nrow(temp)){
      
      # fill enrichmentDF with comparisons and corresponding enrichments
      enrichmentDF <- rbind(
        enrichmentDF, data.frame(
          Comparison = paste0(temp$FileName[i], " : ", temp$FileName[i-1]),
          Query = query,
          Alias = temp$Alias[i],
          Enrichment = temp[i, ncol(temp)] / temp[i-1, ncol(temp)]
        )
      )
    }
  }
  
  return(enrichmentDF)
}