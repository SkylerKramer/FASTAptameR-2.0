#' This function takes a FASTA file (WITH THE FA2 ID) and returns the same information as a dataframe
fa_formatInput <- function(fastaInput = NULL, population = NULL){
  # read data, save IDs, convert sequences to data.frame
  countData <- readLines(fastaInput)
  
  seqIDs <- as.character(countData[seq(1, length(countData), 2)])
  countData <- data.frame(seqs = countData[seq(2, length(countData), 2)], stringsAsFactors = F)
  
  # process seqIDs to make named columns
  idDF <- gsub('>', '', seqIDs) %>%
    strsplit(., '-') %>%
    unlist(.) %>%
    matrix(., nrow = length(seqIDs), byrow = T) %>%
    as.data.frame(.)
  
  # rename columns and optionally append the population value
  colnames(idDF) <- c('Rank', 'Reads', 'RPM', 'Cluster', 'ClusterRank', 'SeedDistance')[1:ncol(idDF)]
  if(!is.null(population)){
    colnames(idDF) <- paste0(colnames(idDF), '.', population)
  }
  
  # treat values as numerics
  idDF[,1] <- as.numeric(as.character(idDF[,1]))
  idDF[,2] <- as.numeric(as.character(idDF[,2]))
  idDF[,3] <- as.numeric(as.character(idDF[,3]))
  
  # return the data.frame
  return(cbind(countData, idDF))
}

#' This function takes a dataframe generated in FA2 and formats it as a FASTA
fa_formatOutput <- function(outputData = NULL){
  # create shape for FASTA output
  fastaOutput <- rep(NA, nrow(outputData)*2)
  
  # add sequence IDs and sequences
  fastaOutput[seq(1,length(fastaOutput),2)] <- outputData$id
  fastaOutput[seq(2,length(fastaOutput),2)] <- as.character(outputData$seqs)
  
  # return the data structure
  return(fastaOutput)
}

#' This function merges redundant amino acid sequences that arose from unique nucleotide sequences
fa_translate_converge <- function(translateDF = NULL){
  # split 'id' column from data.frame input
  seqIDs <- translateDF$id
  
  idDF <- gsub('>', '', seqIDs) %>%
    strsplit(., '-') %>%
    unlist(.) %>%
    matrix(., nrow = length(seqIDs), byrow = T) %>%
    as.data.frame(.)
  colnames(idDF) <- c("Rank", "Reads", "RPM")
  
  # merge data.frame input with 'id' data.frame; remove Rank and id, which will be remade
  translateDF <- cbind(idDF, translateDF) %>% .[,-c(1,4)]
  
  # group by 'seqs' and sum up 'Reads' and 'RPM'
  translateDF <- as.data.frame(translateDF %>%
                                 dplyr::group_by(seqs) %>%
                                 dplyr::summarise(sumReads = sum(as.numeric(Reads)), sumRPM = sum(as.numeric(RPM))) %>%
                                 dplyr::arrange(-sumReads))
  
  # add 'Rank' and 'id' columns
  translateDF$Rank <- row.names(translateDF)
  translateDF$id <- paste0(">", translateDF$Rank, "-", translateDF$sumReads, "-", translateDF$sumRPM)
  
  # return 'id' and 'seqs'
  return(translateDF[,c(5,1)])
}

#' This function determines the genetic code to use for translating and also formats / applies user-defined alterations to that code
fa_translate_mapping <- function(inputChanges = "", translateSelection = "Standard"){
  # convert user input into a usable format
  translationDF_changed <- inputChanges %>% # input string should have 1 codon / translation pair per newline, separated by commas
    gsub("[[:punct:]]", ",", .) %>% # replace all special characters with comma
    gsub(",{2,}", ",", .) %>% # replace consecutive commas with a single one
    gsub(" | \t", "", .) %>% # omit spaces and tabs
    gsub("\\(|\\)", "", .) %>% # omit parentheses
    toupper() %>% # convert all alphabetic characters to uppercase
    strsplit("\n") %>% # split by newline
    unlist() %>% # unlist
    as.data.frame() %>% # convert to data.frame; 1 pair per row
    
    tidyr::separate(., ., c("Codon", "Translation"), ",") %>% # split into multiple columns
    tidyr::drop_na() %>% # drop missing values
    dplyr::mutate(., Codon = gsub("T", "U", Codon)) %>% # replace T with U in Codon column
    dplyr::filter(., nchar(Codon) == 3) %>% # remove codons that are not 3 characters long
    dplyr::filter(., nchar(Translation) == 1) %>% # remove translations that are not 1 characters long
    
    dplyr::mutate(., Codon = as.character(Codon), Translation = as.character(Translation)) # confirm both columns are characters
  
  # select appropriate filepath for JSON file
  filepath <- switch(
    selection,
    "Standard" = "trans/translations_standard.json",
    "Vertebrate mitochondrial" = "trans/translations_vertMito.json",
    "Yeast mitochondrial" = "trans/translations_yeastMito.json",
    "Mold, protozoan, and coelenterate mitochondrial + Mycoplasma / Spiroplasma" = "trans/translations_mold.json",
    "Invertebrate mitochondrial" = "trans/translations_invertMito.json",
    "Ciliate, dasycladacean and Hexamita nuclear" = "trans/translations_cilNucl.json",
    "Echinoderm and flatworm mitochondrial" = "trans/translations_echinoMito.json",
    "Euplotid nuclear" = "trans/translations_euplotidNucl.json",
    "Alternative yeast nuclear" = "trans/translations_altYeastNucl.json",
    "Ascidian mitochondrial" = "trans/translations_ascMito.json",
    "Alternative flatworm mitochondrial" = "trans/translations_altWormMito.json",
    "Blepharisma nuclear" = "trans/translations_blephNucl.json",
    "Chlorophycean mitochondrial" = "trans/translations_chloroMito.json",
    "Trematode mitochondrial" = "trans/translations_tremMito.json",
    "Scenedesmus obliquus mitochondrial" = "trans/translations_sceneMito.json",
    "Pterobranchia mitochondrial" = "trans/translations_pteroMito.json"
  )
  
  # read JSON file with standard genetic code
  translationDF <- jsonlite::fromJSON(filepath) %>% # read JSON file with translation mappings
    dplyr::mutate(., Codon = as.character(Codon), Translation = as.character(Translation)) # confirm both columns are characters
  
  # if input codons and translations are each of length 0, return unmodified data.frame
  if(nrow(translationDF_changed) == 0){
    message("No changes to translation mapping!")
    return(translationDF)
  }
  
  # iterate through dataframe of changes
  for(i in 1:nrow(translationDF_changed)){
    
    # if codon of interest already exists in dataframe, replace the translation
    if(translationDF_changed$Codon[i] %in% translationDF$Codon){
      translationDF[which(translationDF$Codon == translationDF_changed$Codon[i]),2] <- translationDF_changed$Translation[i]
    } else{
      
      # else, append to the data.frame
      translationDF <- rbind(
        data.frame(Codon = translationDF_changed$Codon[i], Translation = translationDF_changed$Translation[i], stringsAsFactors = FALSE),
        translationDF
      )
    }
  }
  
  # return modified mappings
  return(translationDF)
}

#' This function formats the motif query list as a Perl-like regular expression
fa_motif_format <- function(motifList = NULL, motifType = "Nucleotide"){
  # format motif(s) as a Perl query
  # if "Nucleotide" is specified, consider the degenerate codes
  if(motifType == "Nucleotide"){
    motifList <- motifList %>% toupper(.) %>%
      gsub(" ", "", .) %>%
      gsub('U', 'T', .) %>%
      gsub('R', '[AG]', .) %>%
      gsub('Y', '[CT]', .) %>%
      gsub('W', '[AT]', .) %>%
      gsub('S', '[GC]', .) %>%
      gsub('M', '[AC]', .) %>%
      gsub('K', '[GT]', .) %>%
      gsub('B', '[^A]', .) %>%
      gsub('D', '[^C]', .) %>%
      gsub('H', '[^G]', .) %>%
      gsub('V', '[^T]', .) %>%
      gsub('N', '[ACGT]', .) %>%
      gsub(",", "|", .)
  } else{
    motifList <- motifList %>% toupper(.) %>%
      gsub(" ", "", .) %>%
      gsub(",", "|", .)
  }
  
  # return formatted motif(s)
  return(motifList)
}