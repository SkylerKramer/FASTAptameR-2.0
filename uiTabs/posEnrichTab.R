# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

posEnrichTab <- tabPanel(
  "Position Enrichment",
  
  sidebarLayout(
    sidebarPanel(
      # ask for input file
      fileInput("posEnrichInput", label = strong("Choose data to analyze:"),
                multiple = FALSE, placeholder = "Sequence Enrichment CSV file",
                accept = c(".csv")),
      
      # get reference sequence
      textInput("posEnrich_refSequence", label = strong("Reference sequence:"), value = ""),
      shinyBS::bsTooltip("posEnrich_refSequence", "What is your main sequence of interest?"),
      
      # radio buttons for alphabet adjustments
      radioButtons("alphabetMods", label = strong("Do you want to adjust the standard alphabet?"),
                   choices = c("Yes","No"), selected = "No", inline = TRUE),
      
      # only show this panel if the user wants to use non-standard translations
      conditionalPanel(
        condition = "input.alphabetMods == 'Yes'",
        
        # ask for motif
        textAreaInput("posEnrich_mods",
                      label = strong("Comma-separated pairs for replacements or single letters for additions."),
                      height = '100px'),
        shinyBS::bsTooltip("posEnrich_mods", "'C,Z' replaces C with Z, whereas 'Z' adds Z to the alphabet.")
      ),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # enrichment range
      sliderInput("slider_enrichmentRange", label = strong("Range of enrichment values:"),
                  min = 0, max = 200, value = c(0,10), step = 1),
      shinyBS::bsTooltip("slider_enrichmentRange", "What is the range of enrichment values to consider?"),
      
      # get sequence type
      radioButtons("posEnrich_seqType", label = strong("Type of sequences?"),
                   choices = c("Nucleotide", "AminoAcid"), selected = "Nucleotide", inline = TRUE),
      shinyBS::bsTooltip("posEnrich_seqType", "Nucleotide or amino acid sequences?"),
      
      # colour pickers
      colourpicker::colourInput("lowCol", "Select low colour", "red2", returnName = TRUE),
      colourpicker::colourInput("midCol", "Select middle colour", "gold", returnName = TRUE),
      colourpicker::colourInput("highCol", "Select high colour", "yellow", returnName = TRUE),
      
      # text input for users to change breaks
      textAreaInput("posEnrich_breakpoints", label = strong("Comma-separated breakpoints."), placeholder = "7,47"),
      shinyBS::bsTooltip("posEnrich_breakpoints", "For example: 7,47."),
      
      # plot start
      actionButton("posEnrich_start", label = h5("Start"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("posEnrich_start", "Show plots of enrichment per position.")
    ),
    
    # show the datatable that you make in the server
    mainPanel(
      shinycssloaders::withSpinner(plotly::plotlyOutput("avPlot_output", height = "600px")),
      shinycssloaders::withSpinner(plotly::plotlyOutput("heatmap_output", height = "500px"))
    )
  )
)