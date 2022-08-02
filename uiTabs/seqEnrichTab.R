# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

seqEnrichTab <- tabPanel(
  "Sequence Enrichment",
  
  sidebarLayout(
    sidebarPanel(
      
      # ask for input file
      fileInput(
        "enrichInput",
        label = strong("Input data:"),
        multiple = TRUE,
        placeholder = "FASTA files",
        accept = c('.fasta')
      ),
      
      # note on selecting multiple files
      em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
      
      # reorder files because fileInput keeps them in alphabetical order
      selectizeInput("enrich_selectInput", label = strong("Select file order."), choices = "*", multiple = TRUE),
      
      # radiobutton for whether to keep missing sequences
      radioButtons("enrichKeepNA", label = strong("Keep missing sequences?"), choices = c("Yes", "No"), selected = "No", inline = TRUE),
      shinyBS::bsTooltip("enrichKeepNA", "If No, only keep sequences found in all files."),
      
      # start button
      actionButton("enrichStart", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("enrichDownload", label = h5("Download"), style='padding:2px; font-size:80%'),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # fold-change histogram
      actionButton("fcHistStart", label = h5("log2(Enrichment) Histogram"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("fcHistStart", "Histogram of fold-changes for every population comparison"),
      shinyBS::bsModal(
        id = "fcHistWindow",
        title = "log2(Enrichment) Histogram",
        trigger = "fcHistStart",
        size = "large",
        shinycssloaders::withSpinner(plotlyOutput("fcHistOutput"))
      ),
      
      # rpm scatter plot
      actionButton("rpmScatterStart", label = h5("RPM Scatter Plot"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("rpmScatterStart", "Scatter plot of RPMs from each supplied population"),
      shinyBS::bsModal(
        id = "rpmScatterWindow",
        title = "RPM Scatter Plot",
        trigger = "rpmScatterStart",
        size = "large",
        shinycssloaders::withSpinner(plotlyOutput("rpmScatterOutput"))
      ),
      
      # RA plot
      actionButton("raStart", label = h5("RA Plot"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("raStart", "RA plot to compare across input files."),
      shinyBS::bsModal(
        id = "raWindow",
        title = "RA Plot",
        trigger = "raStart",
        size = "large",
        shinycssloaders::withSpinner(plotlyOutput("raOutput"))
      ),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # box plots by cluster
      actionButton("boxplotStart", label = h5("Cluster Boxplot"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("boxplotStart", "Sequence enrichment per cluster."),
      shinyBS::bsModal(
        id = "boxplotWindow",
        title = "Sequence Enrichment per Cluster",
        trigger = "boxplotStart",
        size = "large",
        shinycssloaders::withSpinner(uiOutput("boxplotOutput"))
      )
    ),
    
    # display count output as datatable
    mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("enrichOutput")))
  )
)