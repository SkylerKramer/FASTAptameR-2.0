# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function appends -count to the file name supplied to FASTAptameR-Count
fa_count_outputName <- function(inputFile = NULL, outputType = NULL){
  # this small function extracts the file extension, which later has -count added to it
  rightSubstr <- function(x, n){
    # return the last n characters of string
    return(substr(x, nchar(x) - n + 1, nchar(x)))
  }
  
  # if output extension is 'FASTA', replace previous extension with '-count.fasta'
  # else, replace previous extension with '-count.csv'
  if(outputType == "FASTA"){
    fileName <- ifelse(sum(endsWith(toupper(inputFile), c("FASTA", "FASTQ"))) > 0,
                       sub(rightSubstr(inputFile, 6), "-count.fasta", inputFile),
                       sub(rightSubstr(inputFile, 3), "-count.fasta", inputFile))
  }else{
    fileName <- ifelse(sum(endsWith(toupper(inputFile), c("FASTA", "FASTQ"))) > 0,
                       sub(rightSubstr(inputFile, 6), "-count.csv", inputFile),
                       sub(rightSubstr(inputFile, 3), "-count.csv", inputFile))
  }
  
  # return output file name
  return(fileName)
}

#' This function returns the total number of reads and unique sequences after counting
fa_count_metadata <- function(countData = NULL){
  # get number of (non-)redundant sequences
  uniqueSeqs <- paste0("Unique sequences: ", nrow(countData))
  totalSeqs <- paste0("Total sequences: ", countData$id %>%
                        strsplit(., "-") %>% lapply(., `[[`, 2) %>%
                        unlist(.) %>% as.numeric(.) %>% sum(.))
  
  # return number of (non-)redundant sequences
  return(c(totalSeqs, uniqueSeqs))
}

# patternPerms <- function(patterns = NULL){
#   if(length(patterns) == 1){
#     return(patterns)
#   } else{
#     patternMatrix <- matrix(nrow = 0, ncol = length(patterns))
#     for(i in seq_along(patterns)){
#       patternMatrix <- rbind(patternMatrix, cbind(patterns[i], Recall(patterns[-i])))
#     }
#     return(patternMatrix)
#   }
# }