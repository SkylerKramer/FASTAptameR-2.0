# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

#' This function creates a histogram of fold changes
fa_enrich_histogram <- function(df = NULL){
  
  # make plot
  p <- ggplot2::ggplot(df, ggplot2::aes(log2E_ba)) +
    ggplot2::geom_histogram(color = "black", fill = "skyblue", bins = 30) +
    ggplot2::labs(x = "log2(Enrichment)", y = "Number of Unique Sequences", title = "log2E Histogram - b:a") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45))
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates a scatter plot of RPMs in two populations
fa_enrich_scatter <- function(df = NULL){
  
  # add epsilon factor if RPM is equal to 0; respectively changed RPM1 and RPM2 to RPM.a and RPM.b
  df <- df %>%
    dplyr::mutate(
      RPM.a = ifelse(RPM.a == 0, RPM.a + 0.01, RPM.a),
      RPM.b = ifelse(RPM.b == 0, RPM.b + 0.01, RPM.b)
    )
  
  # make plot
  p <- ggplot2::ggplot(df, ggplot2::aes(RPM.a, RPM.b, text = seqs)) +
    ggplot2::geom_point(colour = "skyblue", alpha = 0.5) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "RPM.a", y = "RPM.b", title = "RPM.a vs RPM.b") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates an RA plot from RPMs in two populations; epsilon = 0.1
fa_enrich_ra <- function(df = NULL){
  # add epsilon factor if RPM is equal to 0
  df <- df %>%
    dplyr::mutate(
      RPM.a = ifelse(RPM.a == 0, RPM.a + 0.01, RPM.a),
      RPM.b = ifelse(RPM.b == 0, RPM.b + 0.01, RPM.b)
    )
  
  # add fold change and average log RPM
  df <- df %>%
    dplyr::mutate(
      R = log2(RPM.b / RPM.a),
      A = 0.5 * log2(RPM.b * RPM.a)
    )
  
  # make plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = A, y = R, text = seqs)) +
    ggplot2::geom_point(colour = "skyblue", alpha = 0.5) +
    ggplot2::labs(x = "Average log2(RPM)", y = "Fold Change", title = "Enrichment RA Plot") +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p)
  
  # return interactive figure
  return(fig)
}

#' This function creates box plots of enrichment per cluster
fa_enrich_clusterBoxplots <- function(df = NULL){
  # save population information
  populations <- gsub(".*_", "", colnames(df)[4]) %>% toupper() %>% strsplit("") %>% unlist()
  
  # rename columns for plotting functionality
  colnames(df) <- c("Seqs", "Cluster", "SeedDistance", "Enrichment")
  
  # get longest cluster number string
  padLength <- df$Cluster %>% as.character() %>% nchar() %>% max()
  
  # left-pad cluster number and sort
  df_format <- df %>%
    dplyr::mutate(Cluster = stringr::str_pad(Cluster, padLength, pad = "0", side = "left")) %>%
    dplyr::arrange(Cluster) %>%
    dplyr::mutate(Cluster = factor(Cluster, levels = unique(Cluster)))
  
  # filter for seed sequences
  df_seeds <- df_format %>% filter(SeedDistance == 0)
  
  # make box plots
  p <- ggplot2::ggplot() +
    ggplot2::geom_boxplot(
      data = df_format,
      ggplot2::aes(Cluster, Enrichment),
      fill = "skyblue", outlier.shape = NA
    ) +
    ggplot2::geom_point(
      data = df_seeds,
      ggplot2::aes(Cluster, Enrichment, text = glue::glue(
        "
      Cluster no.: {Cluster}
      Enrichment: {Enrichment}
      Seed: {Seqs}
      "
      )),
      colour = "red", pch = 18
    ) +
    ggplot2::labs(
      x = "Cluster No.",
      y = "Enrichment",
      title = paste0("Sequence Enrichment by Cluster - ", populations[1], ":", populations[2])
    ) +
    ggplot2::theme_classic()
  
  # make bold text
  p <- boldPlots(p = p)
  
  # make interactive figure
  fig <- plotly::ggplotly(p, tooltip = "text", dynamicTicks = TRUE)
  
  # return interactive figure
  return(fig)
}
