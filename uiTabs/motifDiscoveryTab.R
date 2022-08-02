# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

motifDiscoveryTab <- tabPanel(
  "Discovery",
  
  sidebarLayout(
    sidebarPanel(
      
      # ask for input file
      fileInput(
        "motifDiscovery_input",
        label = strong("Input data*:"),
        multiple = FALSE,
        placeholder = "FASTA file",
        accept = c('.fasta')
      ),
      
      # warning about supported alphabet
      h5("*This currently only supports DNA characters, but alternative alphabets will be supported soon!"),
      
      # slider for minimum number of reads to consider
      sliderInput("motifDiscovery_minReads", label = strong("Min. number of reads to consider:"), min = 0, max = 1000, value = 10, step = 5),
      shinyBS::bsTooltip("motifDiscovery_minReads", "Min. number of reads for sequences to be retained in analyis."),
      
      # set min and max lengths for motifs
      sliderInput(
        "motifDiscovery_lengthRange",
        label = strong("Motif length range:"),
        min = 3, max = 15,
        value = c(5, 10), step = 1
      ),
      shinyBS::bsTooltip("motifDiscovery_lengthRange", "The minimum and maximum string lengths to consider in the analysis."),
      
      # start button
      actionButton("motifDiscovery_start", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("motifDiscovery_download", label = h5("Download"), style='padding:2px; font-size:80%')
    ),
    
    # display output as plot and data.table
    mainPanel(
      
      # show plot of over-enriched strings
      shinycssloaders::withSpinner(plotly::plotlyOutput("motifDiscovery_plot", height = "500px")),
      
      # show over-enriched strings
      shinycssloaders::withSpinner(DT::dataTableOutput("motifDiscovery_output"))
    )
  )
)