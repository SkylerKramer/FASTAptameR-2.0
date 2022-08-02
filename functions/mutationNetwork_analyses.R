#' This function 
fa_mutationalIntermediates <- function(fastaInput = NULL, startNode = NULL, endNode = NULL, maxCost = NULL){
  
  # pull sequences from FASTA file
  seqList <- fa_formatInput(fastaInput = fastaInput) %>% dplyr::pull(seqs)
  
  # add start and end nodes
  seqList <- unique(c(seqList, startNode, endNode))
  
  # initialize data.frame of nodes and connections
  seqConnections <- data.frame(from_vertex = NA, to_vertex = NA, cost = NA)
  rowCount <- 1
  
  # add connections if their edit distances are within the user-defined threshold
  for(i in 1:(length(seqList) - 1)){
    
    j <- i + 1
    
    while(j <= length(seqList)){
      
      editDistance <- as.vector(adist(seqList[i], seqList[j]))
      
      if(editDistance <= maxCost){
        
        seqConnections[rowCount,] <- c(seqList[i], seqList[j], editDistance)
        rowCount <- rowCount + 1
      }
      
      j <- j + 1
    }
  }
  
  # make graph and only keep connections when 
  seqGraph <- cppRouting::makegraph(seqConnections)
  
  # get path from starting node to ending node, allowing only one hop per 
  mutationPath <- tryCatch(
    expr = cppRouting::get_path_pair(seqGraph, from = startNode, to = endNode)[[1]],
    error = function(cond) return(glue::glue("The minimal path for the given constraints is not represented in the dataset!"))
  )
  
  # check if path was not found
  if(length(mutationPath) == 1){
    
    return(mutationPath)
  }
  
  # initialize data.frame for results
  mutationNetwork <- data.frame(From_Sequence = rep(NA, length(mutationPath) - 1), To_Sequence = rep(NA, length(mutationPath) - 1))
  
  # format results into data.frame
  for(i in 2:length(mutationPath)){
    mutationNetwork$From_Sequence[i-1] <- mutationPath[i-1]
    mutationNetwork$To_Sequence[i-1] <- mutationPath[i]
  }
  
  # add LED between sequence pairs
  mutationNetwork$Transition_Cost <- sapply(1:nrow(mutationNetwork), function(i) adist(mutationNetwork$From_Sequence[i], mutationNetwork$To_Sequence[i]) %>% as.vector())
  
  # return the network
  return(mutationNetwork)
}
