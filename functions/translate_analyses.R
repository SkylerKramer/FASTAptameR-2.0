# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

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
  seqs <- trimws(
    toupper(as.character(fastaData[seq(2, length(fastaData), 2)])) %>%
      substring(., orf) %>%
      gsub('T', 'U', .) %>%
      ifelse(nchar(.) %% 3 == 1, gsub(".{1}$", "", .), .) %>%
      ifelse(nchar(.) %% 3 == 2, gsub(".{2}$", "", .), .) %>%
      gsub("(.{3})", "\\1 ", .), which = "right"
  )
  
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
