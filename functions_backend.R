#' This function takes a FASTQ/A file and counts the number of unique sequences
#' These counts and other summary statistics are returned as a dataframe
fa_count <- function(dataInput = NULL){
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

#' This function translates a FA2-formated FASTA file, according to a user-defined genetic code (many options provided)
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
  mergeDF <- mergeDF %>% replace(is.na(.), 0)
  
  # return merged data.frame after ordering by 'Rank' in the first population
  return(mergeDF[order(as.numeric(mergeDF$Rank.a)),])
}

#' This function clusters sequences according to a set of user-defined parameters
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
    dplyr::summarise(TotalSequences = dplyr::n(), TotalReads = sum(Reads), TotalRPM = sum(RPM),
                     AverageLED = round(mean(as.numeric(as.character(SeedDistance))),2)) %>%
    as.data.frame(.)
  
  clusterStats <- clusterStats[order(as.numeric(clusterStats$Cluster)),]
  clusterStats[is.na(clusterStats)] <- 'NC'
  
  # return data.frame with clustering metadata
  return(merge(clDF_format, clusterStats, by = 'Cluster', all = T, sort = F))
}

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