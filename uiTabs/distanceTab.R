# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

distanceTab <- tabPanel(
  "Distance",
  sidebarLayout(
    sidebarPanel(
      # ask for input file
      fileInput("distanceInput", label = strong("Input data:"),
                multiple = FALSE, placeholder = "FASTA file",
                accept = c('.fasta')),
      
      # ask for query sequence
      textInput("querySequence", label = strong("Query sequence:")),
      shinyBS::bsTooltip("querySequence", "Sequence of interest"),
      
      # sequence trunation range; allows user to truncate sequences
      sliderInput("distanceTruncRange", label = strong("Sequence range to be analyzed:"),
                  min = 1, max = 300, value = c(1, 300), step = 1),
      shinyBS::bsTooltip("distanceTruncRange", "Adjust this range if you want to consider specific regions in your sequence."),
      
      # select file type for download
      radioButtons("distanceDownloadType", label = strong("FASTA or CSV download?"),
                   choices = c("FASTA", "CSV"), selected = "FASTA", inline = TRUE),
      shinyBS::bsTooltip("distanceDownloadType",
                         "FASTA is required for subsequent modules; CSV retains all features from data table"),
      
      # start button
      actionButton("distanceStart", label = h5("Start"),
                   style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("distanceDownload", label = h5("Download"),
                     style='padding:2px; font-size:80%'),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # start button
      actionButton("distance_histStart", label = h5("Distance Histogram"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("distance_histStart", "See a histogram of distances"),
      
      tags$head(tags$style("#distance_histWindow .modal-body{ min-height:650px}")),
      shinyBS::bsModal("distance_histWindow", "Distance Histogram", "distance_histStart", size = "large",
                       shinycssloaders::withSpinner(plotly::plotlyOutput("distance_histOutput", height = "600px")))
    ),
    
    # display distance output as datatable
    mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("distanceOutput")))
  )
)