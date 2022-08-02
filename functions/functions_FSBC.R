#' ALL FUNCTIONS IN THIS FILE COME FROM the folowing:
#' "FSBC: FAST STRING-BASED CLUSTERING FOR HT-SELEX DATA" by Kato et al. (DOI: 10.1186/s12859-020-03607-1). 

library(magrittr)
library(Biostrings)

#' Check self overlapping regions of the string.
#'
#' If the string is "ATATA", the set of self overlapping is "A" and "ATA".
#' The following is an example of string "ATATA".
#' Self overlapping with "A".
#' ATATA
#'     ATATA
#'
#' Self overlapping with "ATA".
#' ATATA
#'   ATATA
#'
#' @param x string.
#' @param ... do not use.
#'
#' @return the set of self overlapping region of string.
#'
#' @import magrittr
#'
fsbc_check_self_overlap <- function(x, ...) {
	if(nchar(x) == 1) return(NULL)

	res <- sapply(1:(nchar(x) - 1), function(i) {
			s1 <- substr(x, nchar(x) - i + 1, nchar(x))
			s2 <- substr(x, 1, i)
			c(s1, s2)
	}) %>% t
	res <- res[res[,1] == res[,2],1]

	res
}

#' Calculate the probability of string.
#'
#' If all nucleobases are equal to 1/4, the probability of "ATGC" is (1/4)^4 = 0.00390625.
#'
#' @param x string.
#' @param symbols probability of symbols.
#' @param ... do not use.
#'
#' @return probability of string.
#'
fsbc_calc_prob_subseq <- function(x, symbols, ...) {
	res <- 1
	for(i in 1:nchar(x)) res <- res * symbols[str_sub(x, i, i)]
	res <- as.numeric(res)
	res
}

#' Calculation of probability of sequences including specific string.
#'
#' This function calculates the probability of $L$-mer sequence which includes $l$-mer string.
#' For example, the probabilit of 5-mer sequence which includes "ATGC" (The probabilities of A, T, G, and C are 1/4.) is 0.0781250.
#'
#' @param x string.
#' @param symbols probability of symbols.
#' @param L length of sequence.
#' @param ... do not use.
#'
#' @examples
#'   fsbc_calc_prob_seq("ATGC", c(A = 0.25, T = 0.25, G = 0.25, C = 0.25), L = 5)
#'
#' @return probability of a sequence including specific string
fsbc_calc_prob_seq <- function(x, symbols, L = 40, ...) {
	op <- fsbc_check_self_overlap(x) # Get overlap regino

	P <- c() # Prob of l-mer sequence.
	l <- nchar(x) # Length of string
	p <- fsbc_calc_prob_subseq(x, symbols) # Prob of string

	# Set initial prob.
	P[1:(l-1)] <- 0
	P[l] <- p

	for(i in (l+1):L) {
		P[i] = p * (1 - P[i - l]) + P[i - 1]

		for(j in op) {
			ll <- i - nchar(x) + nchar(j)
			pp <- fsbc_calc_prob_subseq(j, symbols)
			P[i] = P[i] - (p / pp) * (P[ll] - P[ll - 1])
		}
	}

	res <- data.frame(P)
	return(res)
}

#' Extend string with each symbol.
#'
#' This function returns the extended strings.
#' If the string is "ATGC", this function returns {"ATGCA", "ATGCT", "ATGCC", "ATGCG"}.
#'
#' @param x data.frame of string information.
#' @param symbols probabity of symbols.
#' @param ... do not use.
#'
#' @return extended strings
fsbc_extend_subseq <- function(x, symbols, ...) {
	res <- lapply(symbols, function(i) {
			rownames(x) <- paste0(rownames(x), i)
			return(x)
	})
	res <- do.call(rbind, res)
	return(res)
}

#' Get ratio of nucleobases.
#'
#' This function calculates the ratio of nucleobases with all sequence data.
#'
#' @param x DNAStringSet of sequence data.
#' @param symbols DNA_BASES.
#' @param ... do not use.
#'
#' @return ratio of nucleobases.
#'
#' @import Biostrings
#'
#' @export
fsbc_get_base_ratio <- function(x, symbols = DNA_BASES, ...) {
	res <- apply(sapply(symbols, function(i) letterFrequency(x, i)), 2, sum)
	attr(res, "freq") <- res # Keep frequency.
	res <- res / sum(res[DNA_BASES]) # Calcualte ratio.
	return(res)
}

#' Get permutations.
#'
#' @param x vector of symbols.
#' @param n the number of permutations.
#' @param res characters
#' @param ... do not use.
#'
#' @return
fsbc_perm <- function(x, n, res = "", ...) {
	if ( n == 0 ) {
		return(res)
	} else {
		sapply(x, function(i) {
			res <- paste0(res, i)
			fsbc_perm(x, n - 1, res)
		})
	}
}

#' Get the set of initial strings.
#'
#' This function generates all $l$-mer strings.
#' If $l$ is 2, this function generates
#' { AA, AT, AG, AC, TA, TT, TG, TG, GA, GT, GG, GC, CA, CT, CG, CC }.
#'
#' @param symbols symbols (symbols are nucleobases).
#' @param r length of string.
#' @param l length of string.
#' @param ... do not use.
#'
#' @return set of initial strings.
#'
#' @examples
#' fsbc_get_initial_subseq(c("A", "T", "G", "C"), r = 2)
#'
#' @import gtools
#' @import Biostrings
#' @import magrittr
#'
fsbc_get_initial_subseq <- function(symbols, l, r = l, ...) {
	res <- fsbc_perm(symbols, r) %>% as.vector

	# Set $Z$-score as 0.
	df <- data.frame(Z = rep(0, length(res)))
	rownames(df) <- res

	return(df)
}

#' This is a wrapper function of "scale" to avoid errors of NULL.
#'
#' Scale function shows error with NULL.
#' This function return NULL, if the value is NULL for avoiding error of scale function.
#'
#' @param x numerical vector.
#' @param ... parameters for scale function.
#' @return results of scale.
#'
#' @import magrittr
fsbc_scale <- function(x, ...) {
	if(is.null(x)) return(NULL)
	res <- c(scale(x, ...)) %>% set_names(., x)
	return(res)
}

#' Calculate Z-score, of string.
#' This function calculates Z-score with sequence probability and observed frequency.
#' Z-score is provaded with frequency, ratio, and probability of sequences which include string.
#'
#' @param x string.
#' @param y DNAStringSet, sequence dataset.
#' @param symbols probability of symbols.
#' @param freq frequency of \code{y}.
#' @param L length of sequence.
#' @param fixed parameter for vcountPattern.
#' @param ... do not use.
#'
#' @import Biostrings
#' @import magrittr
#'
fsbc_calc_subseq_score <- function(x, y, symbols, freq,
	L = width(y), fixed = T, ...) {

	PL <- lapply(x, fsbc_calc_prob_seq, symbols, L = max(L)) %>%
		do.call(cbind, .) %>% set_names (., x)

	N <- ifelse(length(freq) == 1, length(y), sum(freq))

	P <- sapply(x, function(i) sum(PL[L, i] * freq) / N) # expected probability
	F <- sapply(x, function(i) sum((vcountPattern(i, y, fixed = fixed) > 0) * freq)) # observed frequency.
	R <- F / N # observed ratio
	Z <- (R - P) / sqrt( P * ( 1 - P ) / N ) # Z-score.
	Z[F == 1] <- 0 # If frequency is equal to 0, Z-score will be changed to 0.

	res <- data.frame(F, R, P, Z)
	rownames(res) <- x

	return(res)
}

#' Select strings with eliminating the search space.
#'
#' This function select strings with lengths from $l$-min to $l$-max with eliminating the search space.
#' strings of $l$-min length are all enumerated for initial subset of string searching.
#' The elimination rule of function is $x > y$, where $x$ and $y$ are pre and post extended string of Z-scores.
#' The elimination rule is flexible to change with the parameter \code{fun}.
#'
#' @param x DNAStringSet, non-redundant sequence file.
#' @param freq frequency of sequences of \code{x}. This have to be same order with DNAStringSet \code{x}.
#' @param lmin minimum length of string.
#' @param lmax maximum length of string.
#' @param symbols probability of symbols.
#' @param s0 initial string.
#' @param fun carry over criteria for branch and bound.
#' @param ... functions for fsbc_calc_subseq_score function.
#'
#' @return data.frame of selected strings.
#'
#' @import Biostrings
#' @import magrittr
#'
#' @export
fsbc_search_subseq <- function(x, freq,
	lmin = 5, lmax = 10,
	symbols = fsbc_get_base_ratio(x, symbols = DNA_BASES),
	s0 = fsbc_get_initial_subseq(names(symbols), lmin),
	fun = function(x, y) x > y, ...) {

	res <- list()

	# Select string from lmin to lmax with branch and bound method.
	for(i in lmin:lmax) {
		s1 <- fsbc_calc_subseq_score(rownames(s0), x, symbols, freq, ...)
		rows <- fun(s1$Z, s0$Z)
		s1 <- s1[rows,]
		s1$ZZ <- fsbc_scale(s1$Z)
		s1$L  <- nchar(rownames(s1))
		if(NROW(s1) == 0) return(break)
		s0 <- fsbc_extend_subseq(s1, names(symbols))

		res[[i]] <- s1
	}

	# Reformat to the data.frame object.
	res <- do.call(rbind, res)
	res <- res[sort.list(res$ZZ, decreasing = T),]
	res$rank <- 1:NROW(res)
	attr(res, "symbols") <- symbols

	return(res)
}

#' Generate clusters based on string information.
#'
#' This function generates clusters with ordered strings, recursively.
#'
#' @param x vector of strings.
#' @param y DNAStringSet, sequence dataset.
#' @param res result object for this recursive function. Default is list().
#' @param ... do not use
#'
#' @return list of sequence clusters
#'
#' @import Biostrings
#' @import magrittr
#'
#' @export
fsbc_seq_cluster <- function(x, y, res = list(), ...) {
	if(length(y) < 1 | length(x) < 1 ) return(res)

	rows <- vcountPattern(x[1], y) > 0

	r <- y[rows]
	res <- c(res, list(r)) %>% set_names(., c(names(res), x[1]))
	res <- fsbc_seq_cluster(x[-1], y[!rows], res = res)

	# Exclude small size of clusters less than 2.
	res <- Filter(Negate(function(i) length(i) < 2), res)

	# Clean up the deep layer of list structure to single layer.
	res <- unlist(res, recursive = F)

	return(res)
}

#' Add cluster IDs to names of DNAStringSet dataset.
#'
#' This function adds cluster ID to the DNAStringSet dataset, and also add cluster ID and cluster seeds of subsequece to \code{y}.
#'
#' @param x List of the clusters.
#' @param y DNSStringSet, sequence datset.
#' @param ... do not use
#'
#' @import Biostrings
#' @import magrittr
#'
#' @export
fsbc_label_cluster <- function(x, y, ...) {
	cs <- c()
	id <- c()

	for(i in 1:length(names(x))) {
		rows <- y %in% x[[names(x)[i]]]
		cs[rows] <- names(x)[i]
		id[rows] <- i
	}

	y@metadata <- c(y@metadata, list(cluster.id = id, cluster.seeds = cs))
	names(y) <- paste0(names(y), ".C:", id)

	return(y)
}

#' Calculate frequency of sequence.
#'
#' The frequency, ratio, and ranking are stored in metadata in the DNAStringSet object.
#'
#' @param x DNAStringset, sequence dataset.
#' @param ... do not use.
#'
#' @return DNAStringSet, unique sequence dataset with frequency.
#'
#' @import Biostrings
#'
#' @export
fsbc_calc_freq <- function(x, ...) {
  tab   <- table(as.data.frame(x))
  tab   <- sort(tab, decreasing = T)

  seq   <- names(tab)
  res   <- DNAStringSet(seq)
  freq  <- as.numeric(tab)
  ratio <- as.numeric(tab / sum(tab))
  rank  <- rank(- freq, ties.method = "min")
  names(res) <- paste0(names(res), "R:", rank, ".F:", freq)

  res@metadata <- list(freq = freq, ratio = ratio, rank = rank)
  return(res)
}
