### Skyler Kramer - 2021/01/14
### FASTAptameR 2.0 - user interface

## imports
library(shiny)
library(shinyBS)
library(shinyFiles)
library(colourpicker)

library(dplyr)
library(purrr)
library(ggplot2)
library(plotly)

## source files
source("./functions_support.R")
source("./functions_plots.R")
source("./functions_backend.R")
source("./functions_ui.R")

## change limit for file sizes
options(shiny.maxRequestSize=2000*1024^2)

## sanitize error messages
options(shiny.sanitize.errors = TRUE)

## define ui
ui <- navbarPage("FASTAptameR 2.0",
                 
                 ## application theme
                 theme = shinythemes::shinytheme("cosmo"),
                 
                 ## count function
                 tabPanel("Count",
                          
                          # settings for error notifications
                          tags$head(
                              tags$style(
                                  HTML(".shiny-notification {
                               position:fixed;
                               top: calc(50%);
                               left: calc(50%);
                               }"
                                  )
                              )
                          ),
                          
                          ## color of progress bar
                          tags$head(tags$style(".progress-bar{max-height:100%;}")),
                          
                          sidebarLayout(
                              sidebarPanel(
                                  # ask for input file
                                  fileInput("countInput", label = h5("Choose data to count*:"),
                                            multiple = FALSE, placeholder = "FASTQ or FASTA file",
                                            accept = c(
                                                '.FASTQ', '.FQ', '.fastq', '.fq',
                                                '.FASTA', '.FA', '.fasta', '.fa'
                                            )),
                                  
                                  # optional file upload via github
                                  textInput("countInput_URL", label = h5("Online source:"),
                                            value = "https://raw.githubusercontent.com/SkylerKramer/AptamerLibrary/main/sample100.fastq"),
                                  shinyBS::bsTooltip("countInput_URL", "Link to pre-uploaded data set"),
                                  
                                  # select file type for download
                                  radioButtons("countDownloadType", label = h5("FASTA or CSV download?"),
                                               choices = c("FASTA", "CSV"), selected = "FASTA", inline = TRUE),
                                  shinyBS::bsTooltip("countDownloadType",
                                            "FASTA is required for subsequent modules; CSV retains all features from data table"),
                                  
                                  # start button
                                  actionButton("countStart", label = h5("Start"), style='padding:11px; font-size:80%'),
                                  
                                  # download button
                                  downloadButton("countDownload", label = h5("Download"), style='padding:2px; font-size:80%'),
                                  
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
                                  sliderInput("countSlider_minReads", label = h5("Min. number of reads to plot:"),
                                              min = 0, max = 1000, value = 10, step = 10),
                                  shinyBS::bsTooltip("countSlider_minReads", "What is the min. number of reads to plot?"),
                                  
                                  # slider for maximum rank in "Reads per Rank" plot
                                  sliderInput("countSlider_maxRanks", label = h5("Max. rank to plot:"),
                                              min = 10, max = 1000, value = 100, step = 10),
                                  shinyBS::bsTooltip("countSlider_maxRanks", "How many of the top ranks should be plotted?"),
                                  
                                  # Reads per Rank plot
                                  actionButton("count_rprPlotStart", label = h5("Reads per Rank"), style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("count_rprPlotStart",
                                                     "Shows a line plot comparing reads and ranks of unique sequences"),
                                  shinyBS::bsModal("count_rprPlotWindow", "Reads per Rank", "count_rprPlotStart", size = "large",
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("count_rprPlotOutput"))),
                                  
                                  # horizontal line
                                  tags$hr(style="border-color: black;"),
                                  
                                  # start button for seq. length histogram
                                  actionButton("count_seqHistStart", label = h5("Sequence-Length Histogram"), 
                                               style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("count_seqHistStart", "Shows a histogram of sequence lengths."),
                                  
                                  shinyBS::bsModal("count_seqHistWindow", "Sequence-Length Histogram", "count_seqHistStart", size = "large",
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("count_seqHistOutput", height = "650px"))),
                                  
                                  # horizontal line
                                  tags$hr(style="border-color: black;"),
                                  
                                  # radio button to ask user if they want to adjust the bins in the abundance plot
                                  radioButtons("abundButton", label = h5("Adjust default bins?"),
                                               choices = c("Yes","No"), selected = "No", inline = TRUE),
                                  shinyBS::bsTooltip("abundButton", "Adjust bins in the binned abundance plot."),
                                  
                                  # only show this panel if the user wants to use non-standard translations
                                  conditionalPanel(
                                      condition = "input.abundButton == 'Yes'",
                                      
                                      # radio button to ask user about using singletons
                                      radioButtons("singletonButton", label = h5("Should singletons be a separate category?"),
                                                   choices = c("Yes","No"), selected = "Yes", inline = TRUE),
                                      shinyBS::bsTooltip("singletonButton", "If 'Yes' (DEFAULT), then singletons are a separate category."),
                                      
                                      # text input for users to change breaks
                                      textAreaInput("count_newBreaks", label = h5("Comma-separated breakpoints."), placeholder = "10,100,1000"),
                                      shinyBS::bsTooltip("count_newBreaks", "For example: 10,100,1000.")
                                  ),
                                  
                                  # start button for abundance plot
                                  actionButton("count_abPlotStart", label = h5("Abundance Plot"), 
                                               style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("count_abPlotStart", "Shows the binned abundance of sequences."),
                                  shinyBS::bsModal("count_abPlotWindow", "Binned Abundance Plot", "count_abPlotStart", size = "large",
                                                   shinycssloaders::withSpinner(plotly::plotlyOutput("count_abPlotOutput", height = "650px")))
                              ),
                              
                              # display count output as datatable
                              mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("countOutput")))
                          ),
                 ),
                 tabPanel("Translate",
                          sidebarLayout(
                              sidebarPanel(
                                  # ask for input file
                                  fileInput("translateInput", label = h5("Choose data to translate:"),
                                            multiple = FALSE, placeholder = "FASTA file",
                                            accept = c('.fasta')),
                                  
                                  # radio buttons for ORF
                                  radioButtons("orfButton", label = h5("Open reading frame:"), choices = 1:3, selected = 1, inline = TRUE),
                                  shinyBS::bsTooltip("orfButton", "What is the open reading frame?"),
                                  
                                  # radio buttons for converge
                                  radioButtons("convergeButton", label = h5("Should non-unique sequences be merged?"),
                                               choices = c("Yes","No"), selected = "Yes", inline = TRUE),
                                  shinyBS::bsTooltip("convergeButton", "Merging may show sequence convergence"),
                                  
                                  # dropdown to select genetic code
                                  selectInput(
                                      "translateSelection", label = h5("Genetic code selection:"),
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
                                  radioButtons("nonstandardTranslations", label = h5("Do you want to customize the translation table?"),
                                               choices = c("Yes","No"), selected = "No", inline = TRUE),
                                  
                                  # only show this panel if the user wants to use non-standard translations
                                  conditionalPanel(
                                      condition = "input.nonstandardTranslations == 'Yes'",
                                      
                                      # ask for motif
                                      textAreaInput("translateInput_changes", label = h5("Comma-separated codon / translation pairs."), height = '100px'),
                                      shinyBS::bsTooltip("translateInput_changes", "Please supply one pair per line: e.g., UAG,B.")
                                  ),
                                  
                                  # select file type for download
                                  radioButtons("translateDownloadType", label = h5("FASTA or CSV download?"),
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
                                  sliderInput("translateSlider_minReads", label = h5("Min. number of reads to plot:"),
                                              min = 0, max = 1000, value = 10, step = 10),
                                  shinyBS::bsTooltip("translateSlider_minReads", "What is the min. number of reads to plot?"),
                                  
                                  # slider for maximum rank in "Reads per Rank" plot
                                  sliderInput("translateSlider_maxRanks", label = h5("Max. rank to plot:"),
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
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("translate_seqHistOutput")))
                              ),
                              
                              # display translate output as datatable
                              mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("translateOutput")))
                          ),
                 ),
                 tabPanel("Motif",
                          tabsetPanel(
                              tabPanel("Search",
                                       sidebarLayout(
                                           sidebarPanel(
                                               # ask for input file
                                               fileInput("motifSearchInput", label = h5("Input data:"),
                                                         multiple = FALSE, placeholder = "FASTA file",
                                                         accept = c('.fasta')),
                                               
                                               # ask for motif
                                               textInput("motifInput_search", label = h5("Comma-separated patterns:")),
                                               shinyBS::bsTooltip("motifInput_search", "E.g., aaa,gtg"),
                                               
                                               # radio buttons for parentheses
                                               radioButtons("motifSearchButton_highlight", label = h5("Place patterns in parentheses in output?"),
                                                            choices = c("Yes","No"), selected = "No", inline = TRUE),
                                               shinyBS::bsTooltip("motifSearchButton_highlight", "Does not capture overlapping patterns"),
                                               
                                               # radio buttons for type of filtering
                                               radioButtons("motifSearchButton_partial",
                                                            label = h5("If multiple patterns, return partial matches?"),
                                                            choices = c("Yes","No"), selected = "No", inline = TRUE),
                                               shinyBS::bsTooltip("motifSearchButton_partial",
                                                                  "If No, returned sequences must have ALL patterns."),
                                               
                                               # radio buttons for nt or aa searching
                                               radioButtons("motifSearchButton_motifType", label = h5("Type of pattern?"),
                                                            choices = c("Nucleotide","AminoAcid", "String"), selected = "Nucleotide",
                                                            inline = T),
                                               shinyBS::bsTooltip("motifSearchButton_motifType",
                                                         "Nucleotide option considers degenerate codes; other options do not."),
                                               
                                               # select file type for download
                                               radioButtons("motifSearchDownloadType", label = h5("FASTA or CSV download?"),
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
                                       ),
                              ),
                              tabPanel("Tracker",
                                       sidebarLayout(
                                           sidebarPanel(
                                               # ask for input files
                                               fileInput("motifTrackerInput", label = h5("Input data:"),
                                                         multiple = TRUE, placeholder = "FASTA files",
                                                         accept = c('.fasta')),
                                               
                                               # note on selecting multiple files
                                               em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
                                               
                                               # reorder files because fileInput keeps them in alphabetical order
                                               selectizeInput("motifTracker_selectInput", label = h5("Select file order."),
                                                              choices = "*", multiple = TRUE),
                                               
                                               # ask for motif
                                               textAreaInput("motifInput_query", label = h5("Motif or sequence list:"), height = '100px'),
                                               shinyBS::bsTooltip("motifInput_query", "Please supply one motif or sequence of interest per line."),
                                               
                                               # ask for aliases
                                               textAreaInput("motifInput_alias", label = h5("Alias list (1 per line):"), height = '100px'),
                                               shinyBS::bsTooltip("motifInput_alias", "User-defined IDs for queries."),
                                               
                                               # radio buttons for nt or aa searching
                                               radioButtons("motifTrackerButton_queryType", label = h5("Search for motifs or whole sequences?"),
                                                            choices = c("Motif", "Sequence"), selected = "Motif", inline = TRUE),
                                               shinyBS::bsTooltip("motifTrackerButton_queryType",
                                                                  "Are you searching for motifs or sequences?"),
                                               
                                               conditionalPanel(
                                                   condition = "input.motifTrackerButton_queryType == 'Motif'",
                                                   
                                                   # radio buttons for nt or aa searching
                                                   radioButtons("motifTrackerButton_motifType", label = h5("Type of pattern?"),
                                                                choices = c("Nucleotide", "AminoAcid", "String"), selected = "Nucleotide",
                                                                inline = TRUE),
                                                   shinyBS::bsTooltip("motifTrackerButton_motifType",
                                                                      "Nucleotide option considers degenerate codes; other options do not.")
                                               ),
                                               
                                               # start button
                                               actionButton("motifTrackerStart", label = h5("Start"), style='padding:11px; font-size:80%'),
                                               
                                               # download button for summary
                                               downloadButton("motifTrackerDownload_summary", label = h5("Download Summary"), style='padding:2px; font-size:80%'),
                                               
                                               # download button for enrich
                                               downloadButton("motifTrackerDownload_enrich", label = h5("Download Enrichments"), style='padding:2px; font-size:80%')
                                           ),
                                           
                                           # display distance output as datatable
                                           mainPanel(
                                               shinycssloaders::withSpinner(DT::dataTableOutput("motifTrackerOutput")),
                                               shinycssloaders::withSpinner(DT::dataTableOutput("motifTrackerOutput_enrichmentTable")),
                                               shinycssloaders::withSpinner(plotly::plotlyOutput("motifTrackerPlot"))
                                           )
                                       )
                              )
                          )
                 ),
                 
                 ## distance function
                 tabPanel("Distance",
                          sidebarLayout(
                              sidebarPanel(
                                  # ask for input file
                                  fileInput("distanceInput", label = h5("Input data:"),
                                            multiple = FALSE, placeholder = "FASTA or CSV file",
                                            accept = c('.fasta', '.csv')),
                                  
                                  # ask for query sequence
                                  textInput("querySequence", label = h5("Query sequence:")),
                                  shinyBS::bsTooltip("querySequence", "Sequence of interest"),
                                  
                                  # sequence trunation range; allows user to truncate sequences
                                  sliderInput("distanceTruncRange", label = h5("Sequence range to be analyzed:"),
                                              min = 1, max = 300, value = c(1, 300), step = 1),
                                  shinyBS::bsTooltip("distanceTruncRange", "Adjust this range if you want to consider specific regions in your sequence."),
                                  
                                  # select file type for download
                                  radioButtons("distanceDownloadType", label = h5("FASTA or CSV download?"),
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
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("distance_histOutput")))
                              ),
                              
                              # display distance output as datatable
                              mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("distanceOutput")))
                          ),
                 ),
                 tabPanel("Enrichment",
                          tabsetPanel(
                              tabPanel(
                                  "Sequence Enrichment",
                                  
                                  sidebarLayout(
                                      sidebarPanel(
                                          # ask for input file
                                          fileInput("enrichInput", label = h5("Input data:"),
                                                    multiple = TRUE, placeholder = "FASTA files",
                                                    accept = c('.fasta')),
                                          
                                          # note on selecting multiple files
                                          em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
                                          
                                          # reorder files because fileInput keeps them in alphabetical order
                                          selectizeInput("enrich_selectInput", label = h5("Select file order."),
                                                         choices = "*", multiple = TRUE),
                                          
                                          # radiobutton for whether to keep missing sequences
                                          radioButtons("enrichKeepNA", label = h5("Keep missing sequences?"),
                                                       choices = c("Yes", "No"), selected = "No", inline = TRUE),
                                          shinyBS::bsTooltip("enrichKeepNA", "If No, only keep sequences found in all files."),
                                          
                                          # start button
                                          actionButton("enrichStart", label = h5("Start"),
                                                       style='padding:11px; font-size:80%'),
                                          
                                          # download button
                                          downloadButton("enrichDownload", label = h5("Download"),
                                                         style='padding:2px; font-size:80%'),
                                          
                                          # horizontal line
                                          tags$hr(style="border-color: black;"),
                                          
                                          # slider for minimum number of reads in "Reads per Rank" plot
                                          sliderInput("enrich_minReads", label = h5("Minimum number of reads to consider for persistence plot?"),
                                                      min = 0, max = 1000, value = 0, step = 10),
                                          shinyBS::bsTooltip("enrich_minReads", "What is the min. number of reads to plot?"),
                                          
                                          # sequence persistence barplot
                                          actionButton("seqPersStart", label = h5("Seq. Persistence"), style='padding:11px; font-size:80%'),
                                          shinyBS::bsTooltip("seqPersStart", "Barplot of sequence persistence."),
                                          shinyBS::bsModal("seqPersWindow", "Sequence Persistence Analysis", "seqPersStart", size = "large",
                                                           shinycssloaders::withSpinner(plotly::plotlyOutput("seqPersOutput"))),
                                          
                                          # horizontal line
                                          tags$hr(style="border-color: black;"),
                                          
                                          # fold-change histogram
                                          actionButton("fcHistStart", label = h5("log2(Enrichment) Histogram"), style='padding:11px; font-size:80%'),
                                          shinyBS::bsTooltip("fcHistStart", "Histogram of fold-changes for every population comparison"),
                                          shinyBS::bsModal("fcHistWindow", "log2(Enrichment) Histogram", "fcHistStart", size = "large",
                                                           shinycssloaders::withSpinner(uiOutput("fcHistOutput"))),
                                          
                                          # rpm scatter plot
                                          actionButton("rpmScatterStart", label = h5("RPM Scatter Plot"), style='padding:11px; font-size:80%'),
                                          shinyBS::bsTooltip("rpmScatterStart", "Scatter plot of RPMs from each supplied population"),
                                          shinyBS::bsModal("rpmScatterWindow", "RPM Scatter Plot", "rpmScatterStart", size = "large",
                                                           shinycssloaders::withSpinner(uiOutput("rpmScatterOutput"))),
                                          
                                          # RA plot
                                          actionButton("raStart", label = h5("RA Plot"), style='padding:11px; font-size:80%'),
                                          shinyBS::bsTooltip("raStart",
                                                             "RA plot to compare across input files."),
                                          shinyBS::bsModal("raWindow", "RA Plot", "raStart", size = "large",
                                                           shinycssloaders::withSpinner(uiOutput("raOutput"))),
                                          
                                          # horizontal line
                                          tags$hr(style="border-color: black;"),
                                          
                                          # box plots by cluster
                                          actionButton("boxplotStart", label = h5("Cluster Boxplot"), style='padding:11px; font-size:80%'),
                                          shinyBS::bsTooltip("boxplotStart",
                                                             "Sequence enrichment per cluster."),
                                          shinyBS::bsModal("boxplotWindow", "Sequence Enrichment per Cluster", "boxplotStart", size = "large",
                                                           shinycssloaders::withSpinner(uiOutput("boxplotOutput")))
                                      ),
                                      
                                      # display count output as datatable
                                      mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("enrichOutput")))
                                  )
                              ),
                              
                              tabPanel(
                                  "Position Enrichment",
                                  
                                  sidebarLayout(
                                      sidebarPanel(
                                          # ask for input file
                                          fileInput("posEnrichInput", label = h5("Choose data to analyze:"),
                                                    multiple = FALSE, placeholder = "Sequence Enrichment CSV file",
                                                    accept = c(".csv")),
                                          
                                          # get reference sequence
                                          textInput("posEnrich_refSequence", label = h5("Reference sequence:"), value = ""),
                                          shinyBS::bsTooltip("posEnrich_refSequence", "What is your main sequence of interest?"),
                                          
                                          # radio buttons for alphabet adjustments
                                          radioButtons("alphabetMods", label = h5("Do you want to adjust the standard alphabet?"),
                                                       choices = c("Yes","No"), selected = "No", inline = TRUE),
                                          
                                          # only show this panel if the user wants to use non-standard translations
                                          conditionalPanel(
                                              condition = "input.alphabetMods == 'Yes'",
                                              
                                              # ask for motif
                                              textAreaInput("posEnrich_mods",
                                                            label = h5("Comma-separated pairs for replacements or single letters for additions."),
                                                            height = '100px'),
                                              shinyBS::bsTooltip("posEnrich_mods", "'C,Z' replaces C with Z, whereas 'Z' adds Z to the alphabet.")
                                          ),
                                          
                                          # horizontal line
                                          tags$hr(style="border-color: black;"),
                                          
                                          # enrichment range
                                          sliderInput("slider_enrichmentRange", label = h5("Range of enrichment values:"),
                                                      min = 0, max = 200, value = c(0,10), step = 1),
                                          shinyBS::bsTooltip("slider_enrichmentRange", "What is the range of enrichment values to consider?"),
                                          
                                          # get sequence type
                                          radioButtons("posEnrich_seqType", label = h5("Type of sequences?"),
                                                       choices = c("Nucleotide", "AminoAcid"), selected = "Nucleotide", inline = TRUE),
                                          shinyBS::bsTooltip("posEnrich_seqType", "Nucleotide or amino acid sequences?"),
                                          
                                          # colour pickers
                                          colourpicker::colourInput("lowCol", "Select low colour", "red2", returnName = TRUE),
                                          colourpicker::colourInput("midCol", "Select middle colour", "gold", returnName = TRUE),
                                          colourpicker::colourInput("highCol", "Select high colour", "yellow", returnName = TRUE),
                                          
                                          # text input for users to change breaks
                                          textAreaInput("posEnrich_breakpoints", label = h5("Comma-separated breakpoints."), placeholder = "10,100,1000"),
                                          shinyBS::bsTooltip("posEnrich_breakpoints", "For example: 7,47."),
                                          
                                          # plot start
                                          actionButton("posEnrich_start", label = h5("Start"), style='padding:11px; font-size:80%'),
                                          shinyBS::bsTooltip("posEnrich_start", "Show plots of enrichment per position.")
                                      ),
                                      
                                      # show the datatable that you make in the server
                                      mainPanel(
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("avPlot_output")),
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("heatmap_output"))
                                      )
                                  )
                              )
                          )
                    
                 ),
                 tabPanel("Cluster",
                          tabsetPanel(
                              tabPanel("Cluster",
                                       sidebarLayout(
                                           sidebarPanel(
                                               # ask for input file
                                               fileInput("clusterInput", label = h5("Input data:"),
                                                         multiple = FALSE, placeholder = "FASTA file",
                                                         accept = c('.fasta')),
                                               
                                               # slider for minimum number of reads to cluster
                                               sliderInput("clusterSlider_minReads", label = h5("Min. number of reads to cluster:"),
                                                           min = 0, max = 1000, value = 10, step = 5),
                                               shinyBS::bsTooltip("clusterSlider_minReads",
                                                                  "Min. number of reads for sequences to be clustered"),
                                               
                                               # slider for max LED
                                               sliderInput("clusterSlider_maxLED", label = h5("Max. LED:"),
                                                           min = 0, max = 20, value = 7, step = 1),
                                               shinyBS::bsTooltip("clusterSlider_maxLED", "Max. edit distance from seed sequence"),
                                               
                                               # slider for total number of desired clusters
                                               sliderInput("clusterSlider_totalClusters", label = h5("Max. number of clusters to generate:"),
                                                           min = 5, max = 1000, value = 20, step = 5),
                                               shinyBS::bsTooltip("clusterSlider_totalClusters", "Total number of desired clusters"),
                                               
                                               # button to optionally remove non-clustered sequences
                                               radioButtons("clusterButton_keepNC", label = h5("Keep non-clustered sequences?"),
                                                            choices = c("Yes", "No"), selected = "No", inline = TRUE),
                                               shinyBS::bsTooltip("clusterButton_keepNC",
                                                                  "Non-clustered sequences will have NC in their IDs."),
                                               
                                               # button to determine number of desired outputs
                                               radioButtons("clusterButton_outputs", label = h5("Write output to one file per cluster?"),
                                                            choices = c("Yes","No"), selected = "No", inline = TRUE),
                                               shinyBS::bsTooltip("clusterButton_outputs",
                                                                  "One file per cluster (Yes) or one file for all clusters (No)"),
                                               
                                               # note to user
                                               tags$hr(style="border-color: black;"),
                                               h5("If Yes, please provide an absolute path to a directory below. No output will be displayed if Yes."),
                                               tags$hr(style="border-color: black;"),
                                               
                                               # directory input
                                               textInput("clusterInput_directory", "Directory path:", value = NULL),
                                               shinyBS::bsTooltip("clusterInput_directory",
                                                                  "Only required if you want one cluster per file. Must use backward slashes (`/`)!"),
                                               
                                               # select file type for download
                                               radioButtons("clusterDownloadType", label = h5("FASTA or CSV download?"),
                                                            choices = c("FASTA", "CSV"), selected = "FASTA", inline = TRUE),
                                               shinyBS::bsTooltip("clusterDownloadType",
                                                         "FASTA is required for subsequent modules; CSV retains all features from data table"),
                                               
                                               # start button
                                               actionButton("clusterStart", label = h5("Start"), style='padding:11px; font-size:80%'),
                                               
                                               # download button for summary
                                               downloadButton("clusterDownload", label = h5("Download"), style='padding:2px; font-size:80%'),
                                               
                                               # show console output
                                               shinyjs::useShinyjs(),
                                               strong(textOutput("clusterTextOutput"))
                                           ),
                                           
                                           # display count output as datatable
                                           mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("clusterOutput")))
                                       )
                              ),
                              tabPanel("Diversity",
                                       sidebarLayout(
                                           sidebarPanel(
                                               # ask for input file
                                               fileInput("clusterDiversityInput", label = h5("Input data:"),
                                                         multiple = FALSE, placeholder = "Clustered FASTA file",
                                                         accept = c('.fasta')),
                                               
                                               # start button
                                               actionButton("clusterDiversityStart", label = h5("Start"),
                                                            style='padding:11px; font-size:80%'),
                                               
                                               # download button
                                               downloadButton("clusterDiversityDownload", label = h5("Download"),
                                                              style='padding:2px; font-size:80%'),
                                               
                                               # horizontal line
                                               tags$hr(style="border-color: black;"),
                                               
                                               # start button
                                               actionButton("clusterDiversity_metaplotStart", label = h5("Cluster metaplots"),
                                                            style='padding:11px; font-size:80%'),
                                               shinyBS::bsTooltip("clusterDiversity_metaplotStart",
                                                                  "Plots of unique-/total-sequence count and average LED per cluster"),
                                               
                                               tags$head(tags$style("#clusterDiversity_metaplotWindow .modal-body{ min-height:650px}")),
                                               shinyBS::bsModal("clusterDiversity_metaplotWindow", "Cluster Diversity Metaplots",
                                                                "clusterDiversity_metaplotStart", size = "large",
                                                                shinycssloaders::withSpinner(plotly::plotlyOutput("clusterDiverstiyMetaplotOutput"))),
                                               
                                               # horizontal line
                                               tags$hr(style="border-color: black;"),
                                               
                                               # select value of k for kmer PCA
                                               radioButtons("kmerPCAButton_k", label = h5("kmer size for PCA plot:"),
                                                            choices = 3:5, selected = 3, inline = TRUE),
                                               shinyBS::bsTooltip("kmerPCAButton_k",
                                                                  "Should the PCA be based on 3-, 4-, or 5-mers?"), 
                                               
                                               # slider for top number of clusters to plot
                                               sliderInput("kmerPCASlider_topClusters", label = h5("Number of top clusters to plot:"),
                                                           min = 1, max = 20, value = 10, step = 1),
                                               shinyBS::bsTooltip("kmerPCASlider_topClusters",
                                                                  "How many of the top clusters do you want to plot?"),
                                               
                                               # select value of k for kmer PCA
                                               radioButtons("kmerPCAButton_keepNC", label = h5("Keep non-clustered sequences?"),
                                                            choices = c("Yes", "No"), selected = "No", inline = TRUE),
                                               shinyBS::bsTooltip("kmerPCAButton_keepNC",
                                                                  "If Yes, keep top N clusters AND the non-clustered sequences."), 
                                               
                                               # start button
                                               actionButton("kmerPCAStart", label = h5("k-mer PCA"), style='padding:11px; font-size:80%'),
                                               shinyBS::bsTooltip("kmerPCAStart", "Plot PCA of k-mer matrix, colored by cluster"),
                                               h5("*Characters outside of [A,C,G,T,U] converted to 'X'."),
                                               
                                               shinyBS::bsModal("kmerPCAWindow", "k-mer PCA", "kmerPCAStart", size = "large",
                                                       shinycssloaders::withSpinner(plotly::plotlyOutput("kmerPCAOutput")))
                                           ),
                                           
                                           # display count output as datatable
                                           mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("clusterDiversityOutput")))
                                       )
                              ),
                              tabPanel("Enrichment",
                                       sidebarLayout(
                                           sidebarPanel(
                                               # ask for input file
                                               fileInput("clusterEnrichInput", label = h5("Input data:"),
                                                         multiple = TRUE, placeholder = "Cluster analysis files",
                                                         accept = c('.csv')),
                                               
                                               # note on selecting multiple files
                                               em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
                                               
                                               # reorder files because fileInput keeps them in alphabetical order
                                               selectizeInput("clusterEnrich_selectInput", label = h5("Select file order."),
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
                                               shinycssloaders::withSpinner(DT::dataTableOutput("clusterEnrichOutput_enrichTable")),
                                               shinycssloaders::withSpinner(plotly::plotlyOutput("clusterEnrichTrackingOutput"))
                                           )
                                       ))
                          )
                 ),
                 
                 ## help panel; displays the README file
                 tabPanel("About",
                          # uiOutput("pdfView")
                          sidebarLayout(
                              sidebarPanel(
                                  # module connections
                                  div(img(src = "circosPlot.png", height = 500, width = 500), style = "text-align: center;")
                              ),
                              
                              mainPanel(
                                  # title
                                  h1("About FASTAptameR 2.0"),
                                  
                                  # general description
                                  p(
                                      "FASTAptameR 2.0 (FA2) is an R-based update of FASTAptamer. Like its predecessor, FA2 is an open-source
                                      toolkit designed to analyze populations of sequences resulting from combinatorial selections. This
                                      updated version features a user interface (UI), interactive graphics, more modules, and a faster
                                      implementation of the original clustering algorithm."
                                  ),
                                  
                                  # inputs and outputs
                                  p(HTML(paste0(
                                      "To get started, supply a FASTA/Q file to the Count module. Sample FASTQ files can be downloaded ",
                                      a(href = "https://burkelab.missouri.edu/fastaptamer.html", "here", .noWS = "outside", target = "_blank"),
                                      "."
                                  ))),
                                  
                                  # availability section
                                  h3("Availability"),
                                  
                                  # list of accession points
                                  tags$ul(
                                      tags$li(HTML(paste0(
                                          "Website (direct access): ",
                                          a(href = "https://fastaptamer2.missouri.edu/", "fastaptamer2.missouri.edu",
                                            target = "blank")
                                      ))),
                                      
                                      tags$li(HTML(paste0(
                                          "Docker (local download of platform): ",
                                          a(href = "https://hub.docker.com/repository/docker/skylerkramer/fastaptamer2", "skylerkramer/fastaptamer2",
                                            target = "blank")
                                      ))),
                                      
                                      tags$li(HTML(paste0(
                                          "Github (code download): ",
                                          a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/", "SkylerKramer/FASTAptameR-2.0",
                                            target = "blank")
                                      ))),
                                      
                                      tags$li(HTML(paste0(
                                          "User guide: ",
                                          a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/blob/main/UserGuide.pdf", "UserGuide.pdf",
                                            target = "_blank")
                                      )))
                                  ),
                                  
                                  # citation
                                  h3("Citation"),
                                  p("If you use or modify FA2, please cite: CITATION."),
                                  
                                  # contact section
                                  h3("Contact"),
                                  p(
                                      "FA2 is maintained by Skyler Kramer and Donald Burke. To report bugs or request features, please use one of
                                      the following: "
                                  ),
                                  
                                  # list of contacts
                                  tags$ul(
                                      tags$li(HTML(paste0(
                                          "Github: ",
                                          a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/issues", "SkylerKramer/FASTAptameR-2.0",
                                            target = "blank")
                                      ))),
                                      
                                      tags$li(HTML(paste0(
                                          "Twitter: ",
                                          a(href = "https://twitter.com/BurkeLabRNA", "@BurkeLabRNA", target = "_blank")
                                      ))),
                                      
                                      tags$li("Email: burkelab@missouri.edu, stk7c9@umsystem.edu")
                                  )
                              )
                          )
                 ),
                 
                 ## favicon
                 tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
                 
                 ## source HTML code for GA if file is found (only found if running through web app)
                 if(file.exists("google-analytics.html")){
                     tags$head(includeHTML(("google-analytics.html")))
                 }
)
server <- function(input, output, session) {
    
    ## COUNT - DATA GENERATION
    countDF <- eventReactive(input$countStart, {
        if(is.null(isolate(input$countInput)) & isolate(input$countInput_URL) == ""){
            showNotification("No file or link provided!", type = "error", duration = NULL)
            return(NULL)
        } else{
            # capture output
            withCallingHandlers({
                shinyjs::html("countTextOutput", "")
                fa_count(dataInput = ifelse(!is.null(isolate(input$countInput$datapath)),
                                            isolate(input$countInput$datapath),
                                            isolate(input$countInput_URL)))
            },
            # redirect output to text in UI
            message = function(m){
                shinyjs::html(id = "countTextOutput", html = m$message, add = FALSE)
            })
        }
    })
    
    ## COUNT - DATA DISPLAY
    output$countOutput <- DT::renderDataTable(DT::datatable({
        countDF()
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## COUNT - METADATA
    output$countUI_seqCounts <- renderUI({
        if(is.null(countDF())){
            HTML(paste(""))
        } else{
            HTML(paste(fa_count_metadata(countData = countDF()), collapse = "<br/>"))
        }
    })
    
    ## COUNT - DOWNLOAD
    output$countDownload <- downloadHandler(
        filename = function(){
            fa_count_outputName(inputFile = isolate(input$countInput$name), outputType = isolate(input$countDownloadType))
        },
        
        content = function(file){
            if(isolate(input$countDownloadType) == "FASTA"){
                write.table(fa_formatOutput(outputData = countDF()[input[["countOutput_rows_all"]],]), file,
                            quote = FALSE, row.names = FALSE, col.names = FALSE)
            }else{
                write.csv(countDF()[input[["countOutput_rows_all"]],], file, quote = FALSE, row.names = FALSE)
            }
        }
    )
    
    ## COUNT - READS PER RANK - RENDER
    count_rpr <- eventReactive(input$count_rprPlotStart, {
        if(is.null(countDF())){
            return(NULL)
        } else{
            fa_count_rpr(countData = countDF(),
                         minReads = isolate(input$countSlider_minReads),
                         maxRanks = isolate(input$countSlider_maxRanks))
        }
    })
    
    ## COUNT - READS PER RANK - PLOT
    output$count_rprPlotOutput <- plotly::renderPlotly({
        if(is.null(count_rpr())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            count_rpr()
        }
    })
    
    ## COUNT - SEQUENCE LENGTH HISTOGRAM - RENDER
    count_histogram <- eventReactive(input$count_seqHistStart, {
        if(is.null(countDF())){
            return(NULL)
        } else{
            fa_count_histogram(countData = countDF())
        }
    })
    
    ## COUNT - SEQUENCE LENGTH HISTOGRAM - PLOT
    output$count_seqHistOutput <- plotly::renderPlotly({
        if(is.null(count_histogram())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            count_histogram()
        }
    })
    
    ## COUNT - BINNED ABUNDANCE PLOT - RENDER
    count_abundance <- eventReactive(input$count_abPlotStart, {
        if(is.null(countDF())){
            return(NULL)
        } else if(input$abundButton == "No"){
            fa_count_binnedAbundance(countData = countDF())
        } else if(input$abundButton == "Yes" & grepl("[^0-9,]", gsub("\\s", "", input$count_newBreaks))){
            showNotification("New break points must be numeric and in comma-separated list!")
            return(NULL)
        } else{
            fa_count_binnedAbundance(
                countData = countDF(),
                useSingleton = ifelse(isolate(input$singletonButton) == "Yes", TRUE, FALSE),
                breaks = isolate(input$count_newBreaks) %>% gsub("\\s", "", .) %>% strsplit(., split = ",") %>% unlist() %>% as.numeric()
            )
        }
    })
    
    ## COUNT - BINNED ABUNDANCE PLOT - PLOT
    output$count_abPlotOutput <- plotly::renderPlotly({
        if(is.null(count_abundance())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            count_abundance()
        }
    })
    
    ## TRANSLATE - DATA GENERATION
    translateDF <- eventReactive(input$translateStart, {
        if(is.null(isolate(input$translateInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9,]", gsub("\\s", "", isolate(input$translateInput_changes)))){
            showNotification("Modifications must be alphanumeric!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_translate(fastaInput = isolate(input$translateInput$datapath),
                         orf = isolate(input$orfButton),
                         converge = ifelse(isolate(input$convergeButton) == "Yes", TRUE, FALSE),
                         inputChanges = isolate(input$translateInput_changes),
                         translateSelection = "Standard")
        }
    })
    
    ## TRANSLATE - DATA DISPLAY
    output$translateOutput <- DT::renderDataTable(DT::datatable({
        #translateDF()
        if(sum(is.na(translateDF())) != 0){
            showNotification("Sequences must only have [A, C, G, T, U]!", type = "error", duration = NULL)
            return(NULL)
        } else{
            translateDF()
        }
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## TRANSLATE - DATA DOWNLOAD
    output$translateDownload <- downloadHandler(
        # set filename
        filename = function(){
            inputFile <- input$translateInput$name
            
            # add "-translate" to filename
            fileExt <- rightSubstr(inputFile, 6)
            fname <- sub(fileExt, paste0("-translate", fileExt), inputFile)
            
            if(isolate(input$translateDownloadType) == "FASTA"){
                fname
            }else{
                sub("fasta([^fasta]*)$", "csv", fname)
            }
        },
        
        # set file content
        content = function(file){
            if(isolate(input$translateDownloadType) == "FASTA"){
                write.table(fa_formatOutput(outputData = translateDF()[input[["translateOutput_rows_all"]],]), file,
                            quote = FALSE, row.names = FALSE, col.names = FALSE)
            }else{
                write.csv(translateDF()[input[["translateOutput_rows_all"]],], file, quote = FALSE, row.names = FALSE)
            }
        }
    )
    
    ## TRANSLATE - READS PER RANK - RENDER
    translate_rpr <- eventReactive(input$translate_rprPlotStart, {
        if(is.null(translateDF())){
            return(NULL)
        } else{
            fa_count_rpr(countData = translateDF(),
                         minReads = isolate(input$translateSlider_minReads),
                         maxRanks = isolate(input$translateSlider_maxRanks))
        }
    })
    
    ## TRANSLATE - READS PER RANK - PLOT
    output$translate_rprPlotOutput <- plotly::renderPlotly({
        if(is.null(translate_rpr())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            translate_rpr()
        }
    })
    
    ## TRANSLATE - SEQUENCE LENGTH HISTOGRAM - RENDER
    translate_hist <- eventReactive(input$translate_seqHistStart, {
        if(is.null(translateDF())){
            return(NULL)
        } else{
            fa_count_histogram(countData = translateDF())
        }
    })
    
    ## TRANSLATE - SEQUENCE LENGTH HISTOGRAM - PLOT
    output$translate_seqHistOutput <- plotly::renderPlotly({
        if(is.null(translate_hist())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            translate_hist()
        }
    })
    
    ## MOTIF SEARCH - DATA GENERATION
    motifSearchDF <- eventReactive(input$motifSearchStart, {
        if(is.null(isolate(input$motifSearchInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$motifInput_search) == ""){
            showNotification("Must supply valid pattern(s)!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9]", gsub(",| ", "", isolate(input$motifInput_search)))){
            showNotification("Pattern list must be alphanumeric!", type = "error",
                             duration = NULL)
            return(NULL)
        }else{
            ## function call
            fa_motifSearch(fastaInput = isolate(input$motifSearchInput$datapath),
                           motif = isolate(input$motifInput_search),
                           highlight = ifelse(isolate(input$motifSearchButton_highlight) == "Yes", TRUE, FALSE),
                           partial = ifelse(isolate(input$motifSearchButton_partial) == "Yes", TRUE, FALSE),
                           motifType = isolate(input$motifSearchButton_motifType))
        }
    })
    
    ## MOTIF SEARCH - DATA DISPLAY
    output$motifSearchOutput <- DT::renderDataTable(DT::datatable({
        motifSearchDF()
    },
    filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE,
    options = list(searchHighlight = TRUE, search = list(
        regex = TRUE,
        search = fa_motif_format(motifList = isolate(input$motifInput_search), motifType = isolate(input$motifSearchButton_motifType))
    ))
    ))
    
    ## MOTIF SEARCH - DATA DOWNLOAD
    output$motifSearchDownload <- downloadHandler(
        # set filename
        filename = function(){
            # get input file name
            inputFile <- input$motifSearchInput$name
            
            # add "-motif" to filename
            fileExt <- rightSubstr(inputFile, 6)
            fname <- sub(fileExt, paste0("-motif", fileExt), inputFile)
            
            if(isolate(input$motifSearchDownloadType) == "FASTA"){
                fname
            }else{
                sub("fasta([^fasta]*)$", "csv", fname)
            }
        },
        
        # set file content
        content = function(file){
            # format data for output as FASTA file
            if(isolate(input$motifSearchDownloadType) == "FASTA"){
                write.table(fa_formatOutput(outputData = motifSearchDF()[input[["motifSearchOutput_rows_all"]],]), file,
                            quote = FALSE, row.names = FALSE, col.names = FALSE)
            }else{
                write.csv(motifSearchDF()[input[["motifSearchOutput_rows_all"]],], file,
                          quote = FALSE, row.names = FALSE)
            }
        }
    )
    
    ## MOTIF TRACKER - UPDATE FILE SELECTIONS
    observe({
        updateSelectizeInput(session = session, inputId = "motifTracker_selectInput", choices = input$motifTrackerInput$name)
    })

    ## MOTIF TRACKER - DATA GENERATION
    motifTrackerDF <- eventReactive(input$motifTrackerStart, {
        if(is.null(isolate(input$motifTrackerInput))){
            showNotification("No files provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$motifTrackerInput)) < 2){
            showNotification("Supply at least 2 files!", type = "error", duration = NULL)
            return(NULL)
        } else if(length(isolate(input$motifTracker_selectInput)) < 2){
            showNotification("Please order at least 2 files!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$motifInput_query) == ""){
            showNotification("Must supply valid motif(s)!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9]", gsub(",| |\n", "", isolate(input$motifInput_query)))){
            showNotification("Query list may not contain special characters other than commas!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9]", gsub(" |\n", "", isolate(input$motifInput_alias)))){
            showNotification("Alias list may not contain any special characters!", type = "error", duration = NULL)
            return(NULL)
        } else if(length(unlist(strsplit(isolate(input$motifInput_query), split = "\n"))) != length(unlist(strsplit(isolate(input$motifInput_alias), split = "\n"))) & isolate(input$motifInput_alias) != ""){
            showNotification("When aliases are provided, there must be one for each query!", type = "error", duration = NULL)
            return(NULL)
        } else{
            if(isolate(input$motifTrackerButton_queryType) == "Motif"){
                fa_motif_motifTracker(
                    fastaInputs = isolate(input$motifTrackerInput[match(input$motifTracker_selectInput, input$motifTrackerInput$name),]$datapath),
                    fileNames = isolate(input$motifTrackerInput[match(input$motifTracker_selectInput, input$motifTrackerInput$name),]$name),
                    queryList = isolate(input$motifInput_query),
                    queryAliases = if(isolate(input$motifInput_alias) == ""){NULL} else{isolate(input$motifInput_alias)},
                    motifType = isolate(input$motifTrackerButton_motifType)
                )
            } else{
                fa_motif_sequenceTracker(
                    fastaInputs = isolate(input$motifTrackerInput[match(input$motifTracker_selectInput, input$motifTrackerInput$name),]$datapath),
                    fileNames = isolate(input$motifTrackerInput[match(input$motifTracker_selectInput, input$motifTrackerInput$name),]$name),
                    queryList = isolate(input$motifInput_query),
                    queryAliases = if(isolate(input$motifInput_alias) == ""){NULL} else{isolate(input$motifInput_alias)}
                )
            }
        }
    })
    
    ## MOTIF TRACKER - DATA DISPLAY
    output$motifTrackerOutput <- DT::renderDataTable(DT::datatable({
        motifTrackerDF()
    }, rownames = FALSE
    ))
    
    ## MOTIF TRACKER - DATA DOWNLOAD
    output$motifTrackerDownload_summary <- downloadHandler(
        # set filename
        filename = function(){
            "motifTracker_summary.csv"
        },
        
        # set file content
        content = function(file){
            # format data for output as CSV file
            write.csv(motifTrackerDF(), file, row.names = FALSE, quote = FALSE)
        }
    )
    
    ## MOTIF TRACKER - ENRICHMENT GENERATION
    motifTrackerDF_enrich <- eventReactive(input$motifTrackerStart, {
      if(is.null(motifTrackerDF())){
        return(NULL)
      } else{
        fa_motif_trackerEnrichment(trackerDF = motifTrackerDF())
      }
    })
    
    ## MOTIF TRACKER - ENRICHMENT DISPLAY
    output$motifTrackerOutput_enrichmentTable <- DT::renderDataTable(DT::datatable({
      motifTrackerDF_enrich()
    }, rownames = FALSE
    ))
    
    ## MOTIF TRACKER - ENRICHMENT DOWNLOAD
    output$motifTrackerDownload_enrich <- downloadHandler(
      # set filename
      filename = function(){
        "motifTracker_enrich.csv"
      },
      
      # set file content
      content = function(file){
        # format data for output as CSV file
        write.csv(motifTrackerDF_enrich(), file, row.names = FALSE, quote = FALSE)
      }
    )
    
    ## MOTIF TRACKER - PLOT GENERATION
    trackerPlot <- eventReactive(input$motifTrackerStart, {
        if(is.null(motifTrackerDF())){
            return(NULL)
        } else{
            if(isolate(input$motifTrackerButton_queryType) == "Motif"){
                fa_motif_motifTrackerPlot(targetDF = motifTrackerDF())
            } else{
                fa_motif_sequenceTrackerPlot(targetDF = motifTrackerDF())
            }
        }
    })
    
    ## MOTIF TRACKER - RENDER PLOT
    output$motifTrackerPlot <- plotly::renderPlotly({
        if(is.null(trackerPlot())){
            return(NULL)
        } else{
            trackerPlot()
        }
    })
    
    ## DISTANCE - UPDATE SLIDER BAR
    observe({
      updateSliderInput(
        session = session,
        inputId = "distanceTruncRange",
        min = 1, max = nchar(input$querySequence), 
        value = c(1, nchar(input$querySequence)),
        step = 1
      )
    })
    
    ## DISTANCE - DATA GENERATION
    distanceDF <- eventReactive(input$distanceStart, {
        if(is.null(isolate(input$distanceInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$querySequence) == ""){
            showNotification("Must supply valid query sequence!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9]", gsub(" ", "", isolate(input$querySequence)))){
            showNotification("Query sequence must be alphanumeric!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_distance(dataInput = isolate(input$distanceInput$datapath),
                        querySequence = isolate(input$querySequence),
                        seqRange = isolate(input$distanceTruncRange))
        }
    })
    
    ## DISTANCE - DATA OUTPUT
    output$distanceOutput <- DT::renderDataTable(DT::datatable({
        distanceDF()
    },
    filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## DISTANCE - DATA DOWNLOAD
    output$distanceDownload <- downloadHandler(
        # set filename
        filename = function(){
            # get input file name
            inputFile <- input$distanceInput$name
            
            # add "-distance" to filename
            fileExt <- ifelse(grepl("\\.csv", inputFile), rightSubstr(inputFile, 4), rightSubstr(inputFile, 6))
            fname <- sub(fileExt, paste0("-distance", fileExt), inputFile)
            
            if(isolate(input$distanceDownloadType) == "FASTA"){
                fname
            }else{
                sub("fasta([^fasta]*)$", "csv", fname)
            }
        },
        
        # set file content
        content = function(file){
            if(isolate(input$distanceDownloadType) == "FASTA"){
                write.table(fa_formatOutput(outputData = distanceDF()[input[["distanceOutput_rows_all"]],]), file,
                            quote = FALSE, row.names = FALSE, col.names = FALSE)
            }else{
                write.csv(distanceDF()[input[["distanceOutput_rows_all"]],], file,
                          quote = FALSE, row.names = FALSE)
            }
        }
    )
    
    ## DISTANCE - DISTANCE HISTOGRAM - RENDER
    distance_histogram <- eventReactive(input$distance_histStart, {
        if(is.null(distanceDF())){
            return(NULL)
        } else{
            fa_distance_histogram(distanceData = distanceDF(),
                                  querySequence = isolate(input$querySequence))
        }
    })
    
    ## DISTANCE - DISTANCE HISTOGRAM - PLOT
    output$distance_histOutput <- plotly::renderPlotly({
        if(is.null(distance_histogram())){
            showNotification("Please generate a data table from a single population before plotting!",
                             type = "error", duration = NULL)
            return(NULL)
        } else{
            distance_histogram()
        }
    })
    
    ## SEQUENCE ENRICHMENT - UPDATE FILE SELECTIONS
    observe({
        updateSelectizeInput(session = session, inputId = "enrich_selectInput", choices = input$enrichInput$name)
    })
    
    ## SEQUENCE ENRICHMENT - DATA GENERATION
    enrichDF <- eventReactive(input$enrichStart, {
        if(is.null(isolate(input$enrichInput))){
            showNotification("No files provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$enrichInput)) < 2){
            showNotification("Supply at least 2 files!", type = "error", duration = NULL)
            return(NULL)
        } else if(length(isolate(input$enrich_selectInput)) < 2){
            showNotification("Please order at least 2 files!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_enrich(fastaInputs = isolate(input$enrichInput[match(input$enrich_selectInput, input$enrichInput$name),]$datapath),
                      keepNA = ifelse(isolate(input$enrichKeepNA) == "Yes", TRUE, FALSE))
        }
    })
    
    ## SEQUENCE ENRICHMENT - DATA OUTPUT
    output$enrichOutput <- DT::renderDataTable(DT::datatable({
        enrichDF()
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## SEQUENCE ENRICHMENT - DATA DOWNLOAD
    output$enrichDownload <- downloadHandler(
        # set filename
        filename = function(){
            "enrich.csv"
        },
        
        # set file content
        content = function(file){
            # format data for output as CSV file
            write.csv(enrichDF()[input[["enrichOutput_rows_all"]],], file, row.names = FALSE, quote = FALSE)
        }
    )
    
    ## SEQUENCE ENRICHMENT - GATHER SEQUENCE PERSISTENCE DATA
    enrich_seqPers <- eventReactive(input$seqPersStart, {
        if(is.null(isolate(input$enrichInput))){
            showNotification("No files provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$enrichInput)) < 2 | nrow(isolate(input$enrichInput)) > 3){
            showNotification("Supply 2-3 files!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_enrich_seqPersistence(fastaInputs = isolate(input$enrichInput$datapath), minReads = isolate(input$enrich_minReads))
        }
    })
    
    ## SEQUENCE ENRICHMENT - RENDER SEQUENCE PERSISENCE ANALYSIS
    output$seqPersOutput <- plotly::renderPlotly({
        if(is.null(enrich_seqPers())){
            return(NULL)
        } else{
            enrich_seqPers()
        }
    })
    
    ## SEQUENCE ENRICHMENT - PREP DATA FOR LOG2(ENRICHMENT) HISTOGRAMS
    enrich_histDataPrep <- eventReactive(input$fcHistStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else{
            # select columns with log2E in their names
            log2E <- enrichDF() %>% dplyr::select(., dplyr::starts_with("log2E"))
            
            # make list of selected columns
            dataList <- list()
            for(i in 1:ncol(log2E)){
                dataList[[i]] <- log2E %>% dplyr::select(., i)
                
                # add sequence column
                dataList[[i]]$seqs <- enrichDF()$seqs
            }
            
            # name the list according to the column names
            names(dataList) <- colnames(log2E)
            
            # return dataList
            dataList
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY DATA FOR LOG2(ENRICHMENT) HISTOGRAMS
    observeEvent(input$fcHistStart, {
        output$fcHistOutput <- renderUI({
            # plot requires prepped data
            pd <- req(enrich_histDataPrep())
            
            # create taglist
            tagList(map(
                pd,
                ~ plotly::renderPlotly(fa_enrich_histogram(.))
            ))
        })
    })
    
    ## SEQUENCE ENRICHMENT - PREP DATA FOR RPM SCATTER PLOTS
    enrich_scatterDataPrep <- eventReactive(input$rpmScatterStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else{
            # select columns with RPM in their names
            rpm <- enrichDF() %>% dplyr::select(., dplyr::starts_with("RPM"))
            
            # make list of consecutive column pairs; rename the list elements with the appended letters
            dataList <- list()
            for(i in 1:(ncol(rpm) - 1)){
                dataList[[i]] <- rpm[,c(i,i+1)]
                names(dataList)[i] <- paste0(letters[i], letters[i+1])
                
                # add sequence column
                dataList[[i]]$seqs <- enrichDF()$seqs
            }
            
            # return dataList
            dataList
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY DATA FOR RPM SCATTER PLOTS
    observeEvent(input$rpmScatterStart, {
        output$rpmScatterOutput <- renderUI({
            # plot requires prepped data
            pd <- req(enrich_scatterDataPrep())
            
            # create taglist
            tagList(map(
                pd,
                ~ plotly::renderPlotly(fa_enrich_scatter(.))
            ))
        })
    })
    
    ## SEQUENCE ENRICHMENT - PREP DATA FOR RA PLOTS
    enrich_raDataPrep <- eventReactive(input$raStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else{
            # select columns with RPM in their names
            rpm <- enrichDF() %>% dplyr::select(., dplyr::starts_with("RPM"))
            
            # make list of consecutive column pairs; rename the list elements with the appended letters
            dataList <- list()
            for(i in 1:(ncol(rpm) - 1)){
                dataList[[i]] <- rpm[,c(i,i+1)]
                names(dataList)[i] <- paste0(letters[i], letters[i+1])
                
                # add sequence column
                dataList[[i]]$seqs <- enrichDF()$seqs
            }
            
            # return dataList
            dataList
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY DATA FOR RA PLOTS
    observeEvent(input$raStart, {
        output$raOutput <- renderUI({
            # plot requires prepped data
            pd <- req(enrich_raDataPrep())
            
            # create taglist
            tagList(map(
                pd,
                ~ plotly::renderPlotly(fa_enrich_ra(.))
            ))
        })
    })
    
    ## SEQUENCE ENRICHMENT - PREP DATA FOR CLUSTER BOX PLOTS
    enrich_clusterBoxDataPrep <- eventReactive(input$boxplotStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else if(enrichDF() %>% dplyr::select(dplyr::contains("Cluster")) %>% ncol() == 0){
            showNotification("No cluster information!", type = "error")
            return(NULL)
        } else{
            # get number of populations
            numPops <- enrichDF() %>% dplyr::select(dplyr::contains("Cluster.")) %>% ncol()
            
            # initialize empty list to hold data.frames
            dataList <- list()
            
            # iterate through all populations, starting with the 2nd one
            for(i in 2:numPops){
                # select columns with sequence, cluster, distance to seed, and sequence enrichment
                dataList[[i-1]] <- enrichDF() %>%
                    dplyr::select(
                        "seqs",
                        paste0("Cluster.", letters[i]),
                        paste0("SeedDistance.", letters[i]),
                        paste0("enrichment_", letters[i], letters[i-1])
                    )
            }
            
            # return dataList
            dataList
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY DATA FOR CLUSTER BOX PLOTS
    observeEvent(input$boxplotStart, {
        output$boxplotOutput <- renderUI({
            # plot requires prepped data
            pd <- req(enrich_clusterBoxDataPrep())
            
            # create taglist
            tagList(map(
                pd,
                ~ plotly::renderPlotly(fa_enrich_clusterBoxplots(.))
            ))
        })
    })
    
    ## POSITIONAL ENRICHMENT - AVERAGE ENRICHMENT BARPLOT
    posEnrich_avPlot <- eventReactive(input$posEnrich_start, {
        if(is.null(isolate(input$posEnrichInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9]", gsub(" ", "", isolate(input$posEnrich_refSequence)))){
            showNotification("Reference sequence must be alphanumeric!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9,]", gsub("\\s", "", isolate(input$posEnrich_mods)))){
            showNotification("Modifications must be alphanumeric!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$posEnrich_refSequence) == ""){
            showNotification("No reference sequence provided!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_enrich_avgSequenceBar(
              dataPath = isolate(input$posEnrichInput$datapath),
              refSeq = isolate(input$posEnrich_refSequence),
              modList = isolate(input$posEnrich_mods),
              seqType = isolate(input$posEnrich_seqType),
              enrichRange = isolate(input$slider_enrichmentRange),
              lowCol = isolate(input$lowCol), midCol = isolate(input$midCol), highCol = isolate(input$highCol),
              breakpoints = isolate(input$posEnrich_breakpoints) %>% strsplit(",") %>% unlist() %>% as.numeric()
            )
        }
    })
    
    ## POSITIONAL ENRICHMENT - DISPLAY AVERAGE ENRICHMENT BARPLOT
    output$avPlot_output <- plotly::renderPlotly({
        if(is.null(posEnrich_avPlot())){
            return(NULL)
        } else{
            posEnrich_avPlot()
        }
    })
    
    ## POSITIONAL ENRICHMENT - HEATMAP
    posEnrich_heatMap <- eventReactive(input$posEnrich_start, {
        if(is.null(isolate(input$posEnrichInput))){
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9]", gsub(" ", "", isolate(input$posEnrich_refSequence)))){
            return(NULL)
        } else if(isolate(input$posEnrich_refSequence) == ""){
            return(NULL)
        } else if(grepl("[^a-zA-Z0-9,]", gsub("\\s", "", isolate(input$posEnrich_mods)))){
            return(NULL)
        } else{
            fa_enrich_heatMap(dataPath = isolate(input$posEnrichInput$datapath),
                              refSeq = isolate(input$posEnrich_refSequence),
                              modList = isolate(input$posEnrich_mods),
                              enrichRange = isolate(input$slider_enrichmentRange),
                              seqType = isolate(input$posEnrich_seqType),
                              lowCol = isolate(input$lowCol), midCol = isolate(input$midCol), highCol = isolate(input$highCol))
        }
    })
    
    ## POSITIONAL ENRICHMENT - DISPLAY HEATMAP
    output$heatmap_output <- plotly::renderPlotly({
        if(is.null(posEnrich_heatMap())){
            return(NULL)
        } else{
            posEnrich_heatMap()
        }
    })
    
    ## CLUSTER - DATA GENERATION
    clusterDF <- eventReactive(input$clusterStart, {
        if(is.null(isolate(input$clusterInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$clusterSlider_totalClusters == 0)){
            showNotification("Must generate at least one cluster to proceed!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$clusterButton_outputs) == "Yes" & isolate(input$clusterInput_directory) == ""){
            showNotification("Must supply a directory if multiple outputs are desired!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$clusterButton_outputs) == "Yes" & !dir.exists(isolate(input$clusterInput_directory))){
            showNotification("Directory must already exist!", type = "error", duration = NULL)
            return(NULL)
        } else{
            
            # capture output
            withCallingHandlers({
                shinyjs::html("clusterTextOutput", "")
                fa_clusterLED(fastaInput = isolate(input$clusterInput$datapath),
                              minReads = isolate(input$clusterSlider_minReads),
                              maxLED = isolate(input$clusterSlider_maxLED),
                              totalClusters = isolate(input$clusterSlider_totalClusters),
                              multipleOutputs = ifelse(isolate(input$clusterButton_outputs) == "Yes", TRUE, FALSE),
                              outputDirectory = isolate(input$clusterInput_directory),
                              keepNC = ifelse(isolate(input$clusterButton_keepNC) == "Yes", TRUE, FALSE))
            },
            # redirect output to text in UI
            message = function(m){
                shinyjs::html(id = "clusterTextOutput", html = m$message, add = FALSE)
            })
        }
    })
    
    ## CLUSTER - DATA OUTPUT
    output$clusterOutput <- DT::renderDataTable(DT::datatable({
        clusterDF()
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## CLUSTER - DATA DOWNLOAD
    output$clusterDownload <- downloadHandler(
        # set filename
        filename = function(){
            # get input file name
            inputFile <- input$clusterInput$name

            # add "-cluster" to filename
            fileExt <- rightSubstr(inputFile, 6)
            fname <- sub(fileExt, "-cluster.fasta", inputFile)
            
            if(isolate(input$clusterDownloadType) == "FASTA"){
                fname
            }else{
                sub("fasta([^fasta]*)$", "csv", fname)
            }
        },
        
        # set file content
        content = function(file){
            # format data for output as FASTA file
            if(isolate(input$clusterDownloadType) == "FASTA"){
                write.table(fa_formatOutput(outputData = clusterDF()[input[["clusterOutput_rows_all"]],]), file,
                            quote = FALSE, row.names = FALSE, col.names = FALSE)
            }else{
                write.csv(clusterDF()[input[["clusterOutput_rows_all"]],], file, 
                          quote = FALSE, row.names = FALSE)
            }
        }
    )
    
    ## CLUSTER DIVERSITY - DATA GENERATION
    clusterDiversityDF <- eventReactive(input$clusterDiversityStart, {
        if(is.null(isolate(input$clusterDiversityInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(lengths(strsplit(readLines(isolate(input$clusterDiversityInput$datapath), n = 1), split = "-")) != 6){
            showNotification("Please provided a file from FASTAptameR-Cluster!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_clusterDiversity(clusterFASTA = isolate(input$clusterDiversityInput$datapath))
        }
    })
    
    ## CLUSTER DIVERSITY - DATA OUTPUT
    output$clusterDiversityOutput <- DT::renderDataTable(DT::datatable({
        clusterDiversityDF()
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## CLUSTER DIVERSITY - DATA DOWNLOAD
    output$clusterDiversityDownload <- downloadHandler(
        # set filename
        filename = function(){
            # get input file name
            inputFile <- input$clusterDiversityInput$name
            
            # add "-cluster" to filename
            fileExt <- rightSubstr(inputFile, 6)
            sub(fileExt, "-clusterDiversity.csv", inputFile)
        },
        
        # set file content
        content = function(file){
            # format data for output as CSV file
            write.csv(clusterDiversityDF()[input[["clusterDiversityOutput_rows_all"]],], file, row.names = FALSE, quote = FALSE)
        }
    )
    
    ## CLUSTER DIVERSITY - METAPLOTS - RENDER
    clusterDiversityMetaplot <- eventReactive(input$clusterDiversity_metaplotStart, {
        if(is.null(isolate(input$clusterDiversityInput))){
            return(NULL)
        } else if(lengths(strsplit(readLines(isolate(input$clusterDiversityInput$datapath), n = 1), split = "-")) != 6){
            showNotification("Please provide a file from FASTAptameR-Cluster!", type = "error", duration = NULL)
            return(NULL)
        } else if(is.null(clusterDiversityDF())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_clusterDiversity_metaplot(diversityDF = clusterDiversityDF())
        }
    })
    
    ## CLUSTER DIVERSITY - METADATA - PLOT
    output$clusterDiverstiyMetaplotOutput <- plotly::renderPlotly({
        if(is.null(clusterDiversityMetaplot())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            clusterDiversityMetaplot()
        }
    })
    
    ## CLUSTER DIVERSITY - K-MER PCA - RENDER
    kmerPCA <- eventReactive(input$kmerPCAStart, {
        if(is.null(isolate(input$clusterDiversityInput))){
            return(NULL)
        } else if(lengths(strsplit(readLines(isolate(input$clusterDiversityInput$datapath), n = 1), split = "-")) != 6){
            showNotification("Please provide a file from FASTAptameR-Cluster!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_clusterDiversity_kmerPCA(clusterFile = isolate(input$clusterDiversityInput$datapath),
                                        kmerSize = as.numeric(isolate(input$kmerPCAButton_k)),
                                        topClusters = as.numeric(isolate(input$kmerPCASlider_topClusters)),
                                        keepNC = ifelse(isolate(input$kmerPCAButton_keepNC) == "Yes", TRUE, FALSE))
        }
    })
    
    ## CLUSTER DIVERSITY - K-MER PCA - PLOT
    output$kmerPCAOutput <- plotly::renderPlotly({
        kmerPCA()
    })
    
    ## CLUSTER DIVERSITY - UPDATE FILE SELECTIONS
    observe({
        updateSelectizeInput(session = session, inputId = "clusterEnrich_selectInput", choices = input$clusterEnrichInput$name)
    })
    
    ## CLUSTER ENRICH - DATA GENERATION
    clusterEnrichDF <- eventReactive(input$clusterEnrichStart, {
        if(is.null(isolate(input$clusterEnrichInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$clusterEnrichInput)) < 2){
            showNotification("Supply at least 2 files!", type = "error", duration = NULL)
            return(NULL)
        } else if(length(isolate(input$clusterEnrich_selectInput)) < 2){
            showNotification("Please order at least 2 files!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_clusterEnrich(
                # clusterCSVs = isolate(input$clusterEnrichInput$datapath),
                clusterCSVs = isolate(input$clusterEnrichInput[match(input$clusterEnrich_selectInput, input$clusterEnrichInput$name),]$datapath),
                # fileNames = isolate(input$clusterEnrichInput$name)
                fileNames = isolate(input$clusterEnrichInput[match(input$clusterEnrich_selectInput, input$clusterEnrichInput$name),]$name)
            )
        }
    })
    
    ## CLUSTER ENRICH - DATA OUTPUT
    output$clusterEnrichOutput <- DT::renderDataTable(DT::datatable({
        clusterEnrichDF()
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## CLUSTER ENRICH - SUMMARY DOWNLOAD
    output$clusterEnrichDownload_summary <- downloadHandler(
        # set filename
        filename = "clusterEnrich_summary.csv",
        
        # set file content
        content = function(file){
            # format data for output as CSV file
            write.csv(clusterEnrichDF()[input[["clusterEnrichOutput_rows_all"]],], file, row.names = FALSE, quote = FALSE)
        }
    )
    
    ## CLUSTER ENRICH - ENRICHMENT COMPUTATION
    clusterEnrichOutput_enrichDF <- eventReactive(input$clusterEnrichStart, {
      if(is.null(clusterEnrichDF())){
        return(NULL)
      } else{
        fa_clusterEnrich_tracker(trackerDF = clusterEnrichDF())
      }
    })
    
    ## CLUSTER ENRICH - DATA OUTPUT FOR ENRICHMENT TABLE
    output$clusterEnrichOutput_enrichTable <- DT::renderDataTable(DT::datatable({
      clusterEnrichOutput_enrichDF()
    }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
    ))
    
    ## CLUSTER ENRICH - ENRICH DOWNLOAD
    output$clusterEnrichDownload_enrich <- downloadHandler(
      # set filename
      filename = "clusterEnrich_enrich.csv",
      
      # set file content
      content = function(file){
        # format data for output as CSV file
        write.csv(clusterEnrichOutput_enrichDF()[input[["clusterEnrichOutput_enrichDF_rows_all"]],], file, row.names = FALSE, quote = FALSE)
      }
    )
    
    ## CLUSTER ENRICH - RENDER LINE PLOT
    clusterEnrichTrackingPlot <- eventReactive(input$clusterEnrichStart, {
        if(is.null(clusterEnrichDF())){
            return(NULL)
        } else{
            fa_clusterEnrichTracker(clusterEnrichDF = clusterEnrichDF())
        }
    })
    
    ## CLUSTER ENRICH - DISPLAY LINE PLOT
    output$clusterEnrichTrackingOutput <- plotly::renderPlotly(
        clusterEnrichTrackingPlot()
    )
    
    ## DISPLAY USER GUIDE
    # output$pdfView <- renderUI({
    #     tags$iframe(style="height:800px; width:100%; scrolling=yes",
    #                 src="https://docs.google.com/gview?url=https://raw.githubusercontent.com/SkylerKramer/FASTAptameR-2.0/main/UserGuide.pdf&embedded=true")
    # })
    
    session$onSessionEnded(stopApp)
}

shinyApp(ui = ui, server = server)