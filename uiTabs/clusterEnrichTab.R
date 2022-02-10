# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

clusterEnrichTab <- tabPanel(
  "Enrichment",
  sidebarLayout(
    sidebarPanel(
      # ask for input file
      fileInput("clusterEnrichInput", label = strong("Input data:"),
                multiple = TRUE, placeholder = "Cluster analysis files",
                accept = c('.csv')),
      
      # note on selecting multiple files
      em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
      
      # reorder files because fileInput keeps them in alphabetical order
      selectizeInput("clusterEnrich_selectInput", label = strong("Select file order."),
                     choices = "*", multiple = TRUE),
      
      # start button
      actionButton("clusterEnrichStart", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button for summary
      downloadButton("clusterEnrichDownload_summary", label = h5("Download Summary"), style='padding:2px; font-size:80%'),
      
      # download button for enrich
      downloadButton("clusterEnrichDownload_enrich", label = h5("Download Enrichments"), style='padding:2px; font-size:80%')
    ),
    
    # display count output as datatable
    mainPanel(
      shinycssloaders::withSpinner(DT::dataTableOutput("clusterEnrichOutput")),
      tags$hr(style="border-color: black;"),
      shinycssloaders::withSpinner(DT::dataTableOutput("clusterEnrichOutput_enrichTable")),
      shinycssloaders::withSpinner(plotly::plotlyOutput("clusterEnrichTrackingOutput", height = "800px"))
    )
  )
)