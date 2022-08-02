# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

clusterDiversityTab <- tabPanel(
  "Diversity",
  
  sidebarLayout(
    sidebarPanel(
      
      # ask for input file
      fileInput(
        "clusterDiversityInput",
        label = strong("Input data:"),
        multiple = FALSE,
        placeholder = "Clustered FASTA file",
        accept = c('.fasta')
      ),
      
      # start button
      actionButton("clusterDiversityStart", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("clusterDiversityDownload", label = h5("Download"), style='padding:2px; font-size:80%'),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # start button
      actionButton("clusterDiversity_metaplotStart", label = h5("Cluster metaplots"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("clusterDiversity_metaplotStart", "Plots of unique-/total-sequence count and average LED per cluster"),
      
      tags$head(tags$style("#clusterDiversity_metaplotWindow .modal-body{ min-height:650px}")),
      shinyBS::bsModal(
        id = "clusterDiversity_metaplotWindow",
        title = "Cluster Diversity Metaplots",
        trigger = "clusterDiversity_metaplotStart",
        size = "large",
        shinycssloaders::withSpinner(plotly::plotlyOutput("clusterDiverstiyMetaplotOutput", height = "600px"))
      ),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # select value of k for kmer PCA
      radioButtons("kmerPCAButton_k", label = strong("kmer size for PCA plot:"), choices = 3:5, selected = 3, inline = TRUE),
      shinyBS::bsTooltip("kmerPCAButton_k", "Should the PCA be based on 3-, 4-, or 5-mers?"), 
      
      # slider for top number of clusters to plot
      selectizeInput("kmerPCA_clusters", label = strong("Plot which clusters?"), choices = "*", multiple = TRUE), 
      shinyBS::bsTooltip("kmerPCA_clusters", "Which clusters should be included in the plot?"),
      
      # start button
      actionButton("kmerPCAStart", label = h5("k-mer PCA"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("kmerPCAStart", "Plot PCA of k-mer matrix, colored by cluster"),
      h5("*Characters outside of [A,C,G,T,U] converted to 'X'."),
      
      shinyBS::bsModal(
        id = "kmerPCAWindow",
        title = "k-mer PCA",
        trigger = "kmerPCAStart",
        size = "large",
        shinycssloaders::withSpinner(plotly::plotlyOutput("kmerPCAOutput"))
      )
    ),
    
    # display count output as datatable
    mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("clusterDiversityOutput")))
  )
)