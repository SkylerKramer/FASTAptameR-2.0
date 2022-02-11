# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

translateTab <- tabPanel(
  "Translate",
  sidebarLayout(
    sidebarPanel(
      # ask for input file
      fileInput("translateInput", label = strong("Choose data to translate:"),
                multiple = FALSE, placeholder = "FASTA file",
                accept = c('.fasta')),
      
      # radio buttons for ORF
      radioButtons("orfButton", label = strong("Open reading frame:"), choices = 1:3, selected = 1, inline = TRUE),
      shinyBS::bsTooltip("orfButton", "What is the open reading frame?"),
      
      # radio buttons for converge
      radioButtons("convergeButton", label = strong("Should non-unique sequences be merged?"),
                   choices = c("Yes","No"), selected = "Yes", inline = TRUE),
      shinyBS::bsTooltip("convergeButton", "Merging may show sequence convergence"),
      
      # dropdown to select genetic code
      selectInput(
        "translateSelection", label = strong("Genetic code selection:"),
        choices = c(
          "Standard",
          "Vertebrate mitochondrial",
          "Yeast mitochondrial",
          "Mold, protozoan, and coelenterate mitochondrial + Mycoplasma / Spiroplasma",
          "Invertebrate mitochondrial",
          "Ciliate, dasycladacean and Hexamita nuclear",
          "Echinoderm and flatworm mitochondrial",
          "Euplotid nuclear",
          "Alternative yeast nuclear",
          "Ascidian mitochondrial",
          "Alternative flatworm mitochondrial",
          "Blepharisma nuclear",
          "Chlorophycean mitochondrial",
          "Trematode mitochondrial",
          "Scenedesmus obliquus mitochondrial",
          "Pterobranchia mitochondrial"
        )
      ),
      shinyBS::bsTooltip("translateSelection", "Which code should be used for translating?"),
      
      # radio buttons for non-standard translations
      radioButtons("nonstandardTranslations", label = strong("Do you want to customize the translation table?"),
                   choices = c("Yes","No"), selected = "No", inline = TRUE),
      
      # only show this panel if the user wants to use non-standard translations
      conditionalPanel(
        condition = "input.nonstandardTranslations == 'Yes'",
        
        # ask for motif
        textAreaInput("translateInput_changes", label = strong("Comma-separated codon / translation pairs."), height = '100px'),
        shinyBS::bsTooltip("translateInput_changes", "Please supply one pair per line: e.g., UAG,B.")
      ),
      
      # select file type for download
      radioButtons("translateDownloadType", label = strong("FASTA or CSV download?"),
                   choices = c("FASTA", "CSV"), selected = "FASTA", inline = TRUE),
      shinyBS::bsTooltip("translateDownloadType",
                         "FASTA is required for subsequent modules; CSV retains all features from data table"),
      
      # start button
      actionButton("translateStart", label = h5("Start"),
                   style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("translateDownload", label = h5("Download"),
                     style='padding:2px; font-size:80%'),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # slider for minimum number of reads in "Reads per Rank" plot
      sliderInput("translateSlider_minReads", label = strong("Min. number of reads to plot:"),
                  min = 0, max = 1000, value = 10, step = 10),
      shinyBS::bsTooltip("translateSlider_minReads", "What is the min. number of reads to plot?"),
      
      # slider for maximum rank in "Reads per Rank" plot
      sliderInput("translateSlider_maxRanks", label = strong("Max. rank to plot:"),
                  min = 10, max = 1000, value = 100, step = 10),
      shinyBS::bsTooltip("translateSlider_maxRanks", "How many of the top ranks should be plotted?"),
      
      # Reads per Rank plot
      actionButton("translate_rprPlotStart", label = h5("Reads per Rank"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("translate_rprPlotStart",
                         "Shows a line plot comparing reads and ranks of unique sequences"),
      shinyBS::bsModal("translate_rprPlotWindow", "Reads per Rank", "translate_rprPlotStart", size = "large",
                       shinycssloaders::withSpinner(plotly::plotlyOutput("translate_rprPlotOutput"))),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # start button
      actionButton("translate_seqHistStart", label = h5("Sequence-Length Histogram"),
                   style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("translate_seqHistStart", "See a histogram of sequence lengths in the counted data"),
      shinyBS::bsModal("translate_seqHistWindow", "Sequence-Length Histogram", "translate_seqHistStart",
                       size = "large",
                       shinycssloaders::withSpinner(plotly::plotlyOutput("translate_seqHistOutput", height = "650px")))
    ),
    
    # display translate output as datatable
    mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("translateOutput")))
  )
)