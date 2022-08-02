# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates a PCA scatter plot (first two components) of a kmer matrix
fa_clusterDiversity_kmerPCA <- function(clusterFile = NULL, kmerSize = 3, clustersToPlot = NULL){
  
  # format clustered file
  clusterDF <- fa_formatInput(fastaInput = clusterFile, population = NULL)
  
  # convert ambiguous bases to X
  clusterDF$seqs <- gsub("[^ACGTU]", "X", clusterDF$seqs)
  
  # convert NC values to NA
  clusterDF[clusterDF == "NC"] <- NA
  clusterDF$Cluster <- as.numeric(as.character(clusterDF$Cluster))
  
  # only keep the specified clusters
  clusterDF <- clusterDF[which(clusterDF$Cluster %in% clustersToPlot),]
  
  # compute kmer matrix
  kmerMatrix <- clusterDF$seqs %>%
    as.character(.) %>%
    strsplit(., split = "") %>%
    kmer::kcount(., k = kmerSize, named = FALSE)
  rownames(kmerMatrix) <- NULL
  
  # pca with kmer matrix
  kmerPCA <- prcomp(kmerMatrix, scale. = FALSE)
  
  # make plot
  kmerPCA_plot <- factoextra::fviz_pca_ind(
    kmerPCA,
    addEllipses = F, geom = c("point"),
    col.ind = as.factor(clusterDF$Cluster), legend.title = "Cluster"
  ) +
    ggplot2::labs(title = "k-mer PCA") +
    ggplot2::theme_classic()
  
  # make bold text
  kmerPCA_plot <- boldPlots(p = kmerPCA_plot)
  
  # make interactive figure
  fig <- plotly::ggplotly(kmerPCA_plot)
  
  # return interactive PCA plot
  return(fig)
}

#' This function creates line plots of metadata per cluster
fa_clusterDiversity_metaplot <- function(diversityDF = NULL){
  
  # remove non-clustered sequences and convert cluster to character
  diversityDF <- diversityDF %>%
    dplyr::filter(Cluster != "NC") %>%
    dplyr::mutate(Cluster = as.character(Cluster))
  
  # line plot of unique sequences per cluster
  seqsPlot <- plotly::ggplotly(
    boldPlots(
      ggplot2::ggplot(diversityDF, ggplot2::aes(x = as.numeric(Cluster), y = TotalSequences)) +
        ggplot2::geom_line(colour = "blue", size = 2) +
        ggplot2::labs(x = "Cluster", y = "Seq. Count") +
        ggplot2::theme_bw()
    )
  )
  
  # line plot of total sequences per cluster
  readsPlot <- plotly::ggplotly(
    boldPlots(
      ggplot2::ggplot(diversityDF, ggplot2::aes(x = as.numeric(Cluster), y = TotalReads)) +
        ggplot2::geom_line(colour = "orange", size = 2) +
        ggplot2::labs(x = "Cluster", y = "Read Count") +
        ggplot2::theme_bw()
    )
  )
  
  # line plot of average LED to seed per cluster
  ledPlot <- plotly::ggplotly(
    boldPlots(
      ggplot2::ggplot(diversityDF, ggplot2::aes(x = as.numeric(Cluster), y = AverageLED)) +
        ggplot2::geom_line(colour = "forestgreen", size = 2) +
        ggplot2::labs(x = "Cluster", y = "Avg. LED") +
        ggplot2::theme_bw()
    )
  )
  
  # make interactive figure from subplots
  fig <- plotly::subplot(seqsPlot, readsPlot, ledPlot, titleY = TRUE, shareX = TRUE, nrows = 3) %>%
    plotly::layout(title = "<b>Cluster Metaplots</b>", margin = list(b = 50, l = 50, r = 50, t = 50))
  
  # return interactive figure
  return(fig)
}
