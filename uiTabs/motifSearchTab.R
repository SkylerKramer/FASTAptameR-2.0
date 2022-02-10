# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

motifSearchTab <- tabPanel(
  "Search",
  sidebarLayout(
    sidebarPanel(
      # ask for input file
      fileInput("motifSearchInput", label = strong("Input data:"),
                multiple = FALSE, placeholder = "FASTA file",
                accept = c('.fasta')),
      
      # ask for motif
      textInput("motifInput_search", label = strong("Comma-separated patterns:")),
      shinyBS::bsTooltip("motifInput_search", "E.g., aaa,gtg"),
      
      # radio buttons for parentheses
      radioButtons("motifSearchButton_highlight", label = strong("Place patterns in parentheses in output?"),
                   choices = c("Yes","No"), selected = "No", inline = TRUE),
      shinyBS::bsTooltip("motifSearchButton_highlight", "Does not capture overlapping patterns"),
      
      # radio buttons for type of filtering
      radioButtons("motifSearchButton_partial",
                   label = strong("If multiple patterns, return partial matches?"),
                   choices = c("Yes","No"), selected = "No", inline = TRUE),
      shinyBS::bsTooltip("motifSearchButton_partial",
                         "If No, returned sequences must have ALL patterns."),
      
      # radio buttons for nt or aa searching
      radioButtons("motifSearchButton_motifType", label = strong("Type of pattern?"),
                   choices = c("Nucleotide","AminoAcid", "String"), selected = "Nucleotide",
                   inline = T),
      shinyBS::bsTooltip("motifSearchButton_motifType",
                         "Nucleotide option considers degenerate codes; other options do not."),
      
      # select file type for download
      radioButtons("motifSearchDownloadType", label = strong("FASTA or CSV download?"),
                   choices = c("FASTA", "CSV"), selected = "FASTA", inline = TRUE),
      shinyBS::bsTooltip("motifSearchDownloadType",
                         "FASTA is required for subsequent modules; CSV retains all features from data table"),
      
      # start button
      actionButton("motifSearchStart", label = h5("Start"),
                   style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("motifSearchDownload", label = h5("Download"),
                     style='padding:2px; font-size:80%')
    ),
    
    # display distance output as datatable
    mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("motifSearchOutput")))
  )
)