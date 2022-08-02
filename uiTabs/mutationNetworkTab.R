# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

mutationNetworkTab <- tabPanel(
  "Mutation Network",
  
  sidebarLayout(
    sidebarPanel(
      
      # ask for input file
      fileInput(
        "mutationNetwork_input",
        label = strong("Input data:"),
        multiple = FALSE,
        placeholder = "FASTA file",
        accept = c('.fasta')
      ),
      
      # ask for start node
      textInput("mutationNetwork_startNode", label = strong("Start sequence:")),
      shinyBS::bsTooltip("mutationNetwork_startNode", "Starting point of mutation network."),
      
      # ask for end node
      textInput("mutationNetwork_endNode", label = strong("End sequence:")),
      shinyBS::bsTooltip("mutationNetwork_endNode", "Ending point of mutation network."),
      
      # maximum acceptable cost value
      sliderInput(
        "mutationNetwork_maxCost",
        label = strong("Maximum allowable distance:"),
        min = 1, max = 5,
        value = 1, step = 1
      ),
      shinyBS::bsTooltip("mutationNetwork_maxCost", "What is the maximum acceptable LED to consider between any two consecutive sequences in the network?"),
      
      # start button
      actionButton("mutationNetwork_start", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("mutationNetwork_download", label = h5("Download"), style='padding:2px; font-size:80%')
    ),
    
    mainPanel(
      
      # display output in text box if no path is found
      shinycssloaders::withSpinner(uiOutput("mutationNetwork_text_output")),
      
      # display output in data.table if a path is found
      shinycssloaders::withSpinner(DT::dataTableOutput("mutationNetwork_DT_output"))
    )
  )
)