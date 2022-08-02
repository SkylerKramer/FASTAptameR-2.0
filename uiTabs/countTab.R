# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

countTab <- tabPanel(
  "Count",
  
  # settings for error notifications
  tags$head(tags$style(HTML(".shiny-notification {position:fixed; top: calc(50%); left: calc(50%);}"))),
  
  sidebarLayout(
    sidebarPanel(
      # ask for input file
      fileInput(
        "countInput",
        label = strong("Choose data to count*:"),
        multiple = FALSE,
        placeholder = "FASTQ or FASTA file",
        accept = c(
          '.FASTQ', '.FQ', '.fastq', '.fq',
          '.FASTA', '.FA', '.fasta', '.fa'
        )
      ),
      
      # optional file upload via github
      p(HTML(paste0(
        "Sample data is available from ",
        a(href = "https://github.com/SkylerKramer/AptamerLibrary", "here", .noWS = "outside", target = "_blank"),
        "."
      ))),
      
      # optionally generate reverse complement of sequences
      radioButtons("reverseComplement", label = strong("Return reverse complement of sequences?"), choices = c("Yes", "No"), selected = "No", inline = TRUE),
      shinyBS::bsTooltip("reverseComplement", "Optionally return the reverse complement of the input sequences?"),
      
      # select file type for download
      radioButtons("countDownloadType", label = strong("FASTA or CSV download?"), choices = c("FASTA", "CSV"), selected = "FASTA", inline = TRUE),
      shinyBS::bsTooltip("countDownloadType", "FASTA is required for subsequent modules; CSV retains all features from data table"),
      
      # start button
      actionButton("countStart", label = h5("Start"), style='padding:11px; font-size:80%'),
      
      # download button
      downloadButton("countDownload", label = h5("Download"), style='padding:2px; font-size:80%'),
      
      # add space
      tags$br(),
      tags$br(),
      
      # file upload warning
      strong("*Do not start until loading bar shows 'Upload complete'."),
      
      # Count messages
      strong(uiOutput("countUI_seqCounts")),
      
      # show console output
      shinyjs::useShinyjs(),
      strong(textOutput("countTextOutput")),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # slider for minimum number of reads in "Reads per Rank" plot
      sliderInput(
        "countSlider_minReads",
        label = strong("Min. number of reads to plot:"),
        min = 0, max = 1000,
        value = 10, step = 10
      ),
      shinyBS::bsTooltip("countSlider_minReads", "What is the min. number of reads to plot?"),
      
      # slider for maximum rank in "Reads per Rank" plot
      sliderInput(
        "countSlider_maxRanks",
        label = strong("Max. rank to plot:"),
        min = 10, max = 1000,
        value = 100, step = 10
      ),
      shinyBS::bsTooltip("countSlider_maxRanks", "How many of the top ranks should be plotted?"),
      
      # Reads per Rank plot
      actionButton("count_rprPlotStart", label = h5("Reads per Rank"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("count_rprPlotStart", "Shows a line plot comparing reads and ranks of unique sequences"),
      shinyBS::bsModal(
        id = "count_rprPlotWindow",
        title = "Reads per Rank",
        trigger = "count_rprPlotStart",
        size = "large",
        shinycssloaders::withSpinner(plotly::plotlyOutput("count_rprPlotOutput"))
      ),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # start button for seq. length histogram
      actionButton("count_seqHistStart", label = h5("Sequence-Length Histogram"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("count_seqHistStart", "Shows a histogram of sequence lengths."),
      
      shinyBS::bsModal(
        id = "count_seqHistWindow",
        title = "Sequence-Length Histogram",
        trigger = "count_seqHistStart",
        size = "large",
        shinycssloaders::withSpinner(plotly::plotlyOutput("count_seqHistOutput", height = "650px"))
      ),
      
      # horizontal line
      tags$hr(style="border-color: black;"),
      
      # radio button to ask user if they want to adjust the bins in the abundance plot
      radioButtons("abundButton", label = strong("Adjust default bins?"), choices = c("Yes","No"), selected = "No", inline = TRUE),
      shinyBS::bsTooltip("abundButton", "Adjust bins in the binned abundance plot."),
      
      # only show this panel if the user wants to use non-standard translations
      conditionalPanel(
        condition = "input.abundButton == 'Yes'",
        
        # radio button to ask user about using singletons
        radioButtons(
          "singletonButton",
          label = strong("Should singletons be a separate category?"),
          choices = c("Yes","No"),
          selected = "Yes",
          inline = TRUE
        ),
        shinyBS::bsTooltip("singletonButton", "If 'Yes' (DEFAULT), then singletons are a separate category."),
        
        # text input for users to change breaks
        textAreaInput("count_newBreaks", label = strong("Comma-separated breakpoints."), placeholder = "10,100,1000"),
        shinyBS::bsTooltip("count_newBreaks", "For example: 10,100,1000.")
      ),
      
      # start button for abundance plot
      actionButton("count_abPlotStart", label = h5("Abundance Plot"), style='padding:11px; font-size:80%'),
      shinyBS::bsTooltip("count_abPlotStart", "Shows the binned abundance of sequences."),
      shinyBS::bsModal(
        id = "count_abPlotWindow",
        title = "Binned Abundance Plot",
        trigger = "count_abPlotStart",
        size = "large",
        shinycssloaders::withSpinner(plotly::plotlyOutput("count_abPlotOutput", height = "650px"))
      )
    ),
    
    # display count output as datatable
    mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("countOutput")))
  )
)