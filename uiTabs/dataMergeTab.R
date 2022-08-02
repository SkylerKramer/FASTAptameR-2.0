# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

dataMergeTab <- tabPanel(
  "Data Merge",
  
  sidebarLayout(
    
    sidebarPanel(
      
      # ask for input files
      fileInput(
        "dataMerge_input",
        label = strong("Input data:"),
        multiple = TRUE,
        placeholder = "FASTA files",
        accept = c('.fasta')
      ),
      
      # note on selecting multiple files
      em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
      
      # reorder files because fileInput keeps them in alphabetical order
      selectizeInput("dataMerge_selectInput", label = strong("Select file order."), choices = "*", multiple = TRUE),
      
      # radio buttons to define desired merge type
      radioButtons(
        "dataMerge_mergeType",
        label = strong("How should the data be joined?"),
        choices = c("Union", "Intersection", "Left"),
        selected = "Union",
        inline = TRUE
      ),
      shinyBS::bsTooltip(
        "dataMerge_mergeType",
        "Do you want all sequences (Union), shared sequences (intersection), or only sequences from the first population (Left)?"
      ),
      
      # start button
      actionButton("dataMerge_start", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button for summary
      downloadButton("dataMerge_download", label = h5("Download"), style='padding:2px; font-size:80%')
    ),
    
    mainPanel(
      
      # row for plots
      # fluidRow(
      #   column(5, shinycssloaders::withSpinner(plotly::plotlyOutput("dataMerge_sequencePersistence"))),
      #   column(7, shinycssloaders::withSpinner(plotOutput("dataMerge_UpSetR")))
      # ),
      
      shinycssloaders::withSpinner(plotly::plotlyOutput("dataMerge_sequencePersistence")),
      shinycssloaders::withSpinner(plotOutput("dataMerge_UpSetR")),
      
      # merge results
      shinycssloaders::withSpinner(DT::dataTableOutput("dataMerge_output"))
    )
  )
)