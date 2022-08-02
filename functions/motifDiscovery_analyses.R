# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

# All code in "./functions/functions_FSBC.R" comes from the following:
# "FSBC: FAST STRING-BASED CLUSTERING FOR HT-SELEX DATA" by Kato et al. (DOI: 10.1186/s12859-020-03607-1). 
source("./functions/functions_FSBC.R")

# This function is a wrapper around the FSBC code, which identifies over-enriched substrings from an input sequence set
fa_fsbc_motifDiscovery <- function(fastaInput = NULL, minReads = 10, lengthRange = c(5, 10)){
  
  # read FASTA file
  inputDF <- fa_formatInput(fastaInput = fastaInput) %>% dplyr::filter(Reads >= minReads)
  
  # make DNAStringSet object from sequences
  sequenceSet <- inputDF$seqs %>% Biostrings::DNAStringSet()
  
  # calculate sequence frequencies
  sequenceFreqs <- fsbc_calc_freq(sequenceSet)
  
  # calculate base frequencies
  baseFreqs <- fsbc_get_base_ratio(sequenceSet)
  
  # find over-enriched subsequences (motifs)
  motifSet <- fsbc_search_subseq(
    x = sequenceFreqs,
    freq = sequenceFreqs@metadata$freq,
    symbols = baseFreqs, lmin = min(lengthRange), lmax = max(lengthRange)
  ) %>%
    tibble::rownames_to_column("Motif") %>%
    dplyr::mutate(
      R = round(R, digits = 3),
      P = round(P, digits = 5),
      Z = round(Z, digits = 3),
      ZZ = round(ZZ, digits = 3)
    ) %>%
    dplyr::rename(Motif_Length = L, Rank = rank) %>%
    dplyr::select(Motif, P, ZZ, Motif_Length, Rank)
  
  # return set of motifs
  return(motifSet)
}
