# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function takes a FASTQ/A file and counts the number of unique sequences
#' These counts and other summary statistics are returned as a dataframe
fa_count <- function(dataInput = NULL, reverseComplement = NULL){
  # get start time
  startTime <- Sys.time()
  
  # read FASTQ/A file and keep the sequences in a data.frame
  if(toupper(sub(".*(?=.$)", "", dataInput, perl=TRUE)) == "Q"){
    
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
  
  # optionally make reverse complement
  if(reverseComplement){
    seqCounts$seqs <- Biostrings::DNAStringSet(seqCounts$seqs) %>%
      Biostrings::reverseComplement() %>%
      as.character()
  }
  
  # message elapsed time
  message(gsub("Time difference of", "Elapsed time:", capture.output(round(Sys.time() - startTime, 2))))
  
  # return data.frame
  return(seqCounts)
}
