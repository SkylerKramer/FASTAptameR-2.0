### Skyler Kramer - 2021/01/14
### FASTAptameR 2.0 - user interface

## imports
library(shiny)
library(shinyBS)
source("./Functions.R")

## change limit for file sizes
options(shiny.maxRequestSize=2000*1024^2)

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
                          
                          sidebarLayout(
                              sidebarPanel(
                                  # ask for input file
                                  fileInput("countInput", label = h5("Choose data to count*:"),
                                            multiple = F, placeholder = "FASTQ or FASTA file",
                                            accept = c('.FASTQ', '.FQ', '.fastq', '.fq',
                                                       '.FASTA', '.FA', '.fasta', '.fa')),
                                  
                                  # change color of progress bar
                                  tags$style(".progress-bar {background-color: red;}"),
                                  
                                  # file upload warning
                                  strong("*Do not start until loading bar shows 'Upload complete'."),
                                  
                                  # optional file upload via github
                                  textInput("countInput_URL", label = h5("Online source:"),
                                            value = "https://raw.githubusercontent.com/SkylerKramer/AptamerLibrary/main/sample100.fastq"),
                                  shinyBS::bsTooltip("countInput_URL", "Link to pre-uploaded data set"),
                                  
                                  # select file type for download
                                  radioButtons("countDownloadType", label = h5("FASTA or CSV download?"),
                                               choices = c("FASTA", "CSV"), selected = "FASTA", inline = T),
                                  shinyBS::bsTooltip("countDownloadType",
                                            "FASTA is required for subsequent modules; CSV retains all features from data table"),
                                  
                                  # start button
                                  actionButton("countStart", label = h5("Start"), style='padding:11px; font-size:80%'),
                                  
                                  # download button
                                  downloadButton("countDownload", label = h5("Download"), style='padding:2px; font-size:80%'),
                                  
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
                                  
                                  # start button
                                  actionButton("count_seqHistStart", label = h5("Sequence-Length Histogram"), 
                                               style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("count_seqHistStart", "See a histogram of sequence lengths in the counted data"),
                                  
                                  tags$head(tags$style("#count_seqHistWindow .modal-body{ min-height:650px}")),
                                  shinyBS::bsModal("count_seqHistWindow", "Sequence-Length Histogram", "count_seqHistStart", size = "large",
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("count_seqHistOutput")))
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
                                            multiple = F, placeholder = "FASTA file",
                                            accept = c('.fasta')),
                                  
                                  # radio buttons for ORF
                                  radioButtons("orfButton", label = h5("Open reading frame:"), choices = 1:3, selected = 1, inline = T),
                                  shinyBS::bsTooltip("orfButton", "What is the open reading frame?"),
                                  
                                  # radio buttons for converge
                                  radioButtons("convergeButton", label = h5("Should non-unique sequences be merged?"),
                                               choices = c("Yes","No"), selected = "Yes", inline = T),
                                  shinyBS::bsTooltip("convergeButton", "Merging may show sequence convergence"),
                                  
                                  # select file type for download
                                  radioButtons("translateDownloadType", label = h5("FASTA or CSV download?"),
                                               choices = c("FASTA", "CSV"), selected = "FASTA", inline = T),
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
                                                         multiple = F, placeholder = "FASTA file",
                                                         accept = c('.fasta')),
                                               
                                               # ask for motif
                                               textInput("motifInput_search", label = h5("Comma-separated patterns:")),
                                               shinyBS::bsTooltip("motifInput_search", "E.g., aaa,gtg"),
                                               
                                               # radio buttons for parentheses
                                               radioButtons("motifSearchButton_highlight", label = h5("Place patterns in parentheses in output?"),
                                                            choices = c("Yes","No"), selected = "No", inline = T),
                                               shinyBS::bsTooltip("motifSearchButton_highlight", "Does not capture overlapping patterns"),
                                               
                                               # radio buttons for type of filtering
                                               radioButtons("motifSearchButton_partial",
                                                            label = h5("If multiple patterns, return partial matches?"),
                                                            choices = c("Yes","No"), selected = "No", inline = T),
                                               shinyBS::bsTooltip("motifSearchButton_partial",
                                                                  "If No, returned sequences must have ALL patterns."),
                                               
                                               # radio buttons for nt or aa searching
                                               radioButtons("motifSearchButton_motifType", label = h5("Type of pattern?"),
                                                            choices = c("Nucleotide","AminoAcid", "String"), selected = "Nucleotide",
                                                            inline = T),
                                               shinyBS::bsTooltip("motifSearchButton_motifType",
                                                         "Nucleotide option considers degenerate codes; other options do not"),
                                               
                                               # select file type for download
                                               radioButtons("motifSearchDownloadType", label = h5("FASTA or CSV download?"),
                                                            choices = c("FASTA", "CSV"), selected = "FASTA", inline = T),
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
                                                         multiple = T, placeholder = "FASTA files",
                                                         accept = c('.fasta')),
                                               
                                               # note on selecting multiple files
                                               em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
                                               
                                               # ask for motif
                                               textInput("motifInput_enrich", label = h5("Comma-separated patterns:")),
                                               shinyBS::bsTooltip("motifInput_enrich", "E.g., aaa,gtg"),
                                               
                                               # radio buttons for nt or aa searching
                                               radioButtons("motifTrackerButton_motifType", label = h5("Type of pattern?"),
                                                            choices = c("Nucleotide", "AminoAcid", "String"), selected = "Nucleotide",
                                                            inline = T),
                                               shinyBS::bsTooltip("motifTrackerButton_motifType",
                                                                  "Yes for nucleotides; No for amino acids"),
                                               
                                               # start button
                                               actionButton("motifTrackerStart", label = h5("Start"),
                                                            style='padding:11px; font-size:80%'),
                                               
                                               # download button
                                               downloadButton("motifTrackerDownload", label = h5("Download"),
                                                              style='padding:2px; font-size:80%'),
                                               
                                               # Enrichment messages
                                               strong(uiOutput("motifTrackerUI"))
                                           ),
                                           
                                           # display distance output as datatable
                                           mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("motifTrackerOutput")))
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
                                            multiple = F, placeholder = "FASTA file",
                                            accept = c('.fasta')),
                                  
                                  # ask for query sequence
                                  textInput("querySequence", label = h5("Query sequence:")),
                                  shinyBS::bsTooltip("querySequence", "Sequence of interest"),
                                  
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
                 tabPanel("Sequence enrichment",
                          sidebarLayout(
                              sidebarPanel(
                                  # ask for input file
                                  fileInput("enrichInput", label = h5("Input data:"),
                                            multiple = T, placeholder = "FASTA files",
                                            accept = c('.fasta')),
                                  
                                  # note on selecting multiple files
                                  em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
                                  
                                  # radio buttons for parentheses
                                  radioButtons("enrichButton", label = h5("Remove sequences not found in all populations?"),
                                               choices = c("Yes","No"), selected = "Yes", inline = T),
                                  shinyBS::bsTooltip("enrichButton", "These sequences will not be included in plots or calculations."),

                                  # start button
                                  actionButton("enrichStart", label = h5("Start"),
                                               style='padding:11px; font-size:80%'),
                                  
                                  # download button
                                  downloadButton("enrichDownload", label = h5("Download"),
                                                 style='padding:2px; font-size:80%'),
                                  
                                  # note to user
                                  h5("*Numeric columns filtered by slider or typing (e.g., 1 ... 10)"),
                                  
                                  # horizontal line
                                  tags$hr(style="border-color: black;"),
                                  
                                  # fold-change histogram
                                  actionButton("fcHistStart", label = h5("log2(Enrichment) Histogram"), style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("fcHistStart", "Histogram of fold-changes for every population comparison"),
                                  shinyBS::bsModal("fcHistWindow", "log2(Enrichment) Histogram", "fcHistStart", size = "large",
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("fcHistOutput"))),

                                  # horizontal line
                                  tags$hr(style="border-color: black;"),
                                  
                                  # rpm scatter plot
                                  actionButton("rpmScatterStart", label = h5("RPM Scatter Plot"), style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("rpmScatterStart", "Scatter plot of RPMs from each supplied population"),
                                  shinyBS::bsModal("rpmScatterWindow", "RPM Scatter Plot", "rpmScatterStart", size = "large",
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("rpmScatterOutput"))),
                                  
                                  # horizontal line
                                  tags$hr(style="border-color: black;"),
                                  
                                  # volcano plot
                                  actionButton("volcanoStart", label = h5("Volcano Plot"), style='padding:11px; font-size:80%'),
                                  shinyBS::bsTooltip("volcanoStart",
                                                     "Volcano plot of fold-change significance for every population comparison"),
                                  shinyBS::bsModal("volcanoWindow", "Volcano Plot", "volcanoStart", size = "large",
                                          shinycssloaders::withSpinner(plotly::plotlyOutput("volcanoOutput")))
                              ),
                              
                              # display count output as datatable
                              mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("enrichOutput")))
                          )
                 ),
                 tabPanel("Cluster",
                          tabsetPanel(
                              tabPanel("Cluster",
                                       sidebarLayout(
                                           sidebarPanel(
                                               # ask for input file
                                               fileInput("clusterInput", label = h5("Input data:"),
                                                         multiple = F, placeholder = "FASTA file",
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
                                               sliderInput("clusterSlider_totalClusters", label = h5("Total clusters:"),
                                                           min = 0, max = 1000, value = 20, step = 5),
                                               shinyBS::bsTooltip("clusterSlider_totalClusters", "Total number of desired clusters"),
                                               
                                               # button to optionally remove non-clustered sequences
                                               radioButtons("clusterButton_keepNC", label = h5("Keep non-clustered sequences?"),
                                                            choices = c("Yes", "No"), selected = "No", inline = T),
                                               shinyBS::bsTooltip("clusterButton_keepNC",
                                                                  "Non-clustered sequences will have NC in their IDs."),
                                               
                                               # button to determine number of desired outputs
                                               radioButtons("clusterButton_outputs", label = h5("One file per cluster?"),
                                                            choices = c("Yes","No"), selected = "No", inline = T),
                                               shinyBS::bsTooltip("clusterButton_outputs",
                                                                  "One file per cluster (Yes) or one file for all clusters (No)"),
                                               
                                               # note to user
                                               tags$hr(style="border-color: black;"),
                                               h5("If Yes, please provide an absolute path to a directory below. No output will be displayed if Yes."),
                                               tags$hr(style="border-color: black;"),
                                               
                                               # directory input
                                               textInput("clusterInput_directory", "Directory path:", value = NULL),
                                               shinyBS::bsTooltip("clusterInput_directory",
                                                                  "Only required if you want one cluster per file"),
                                               
                                               # select file type for download
                                               radioButtons("clusterDownloadType", label = h5("FASTA or CSV download?"),
                                                            choices = c("FASTA", "CSV"), selected = "FASTA", inline = T),
                                               shinyBS::bsTooltip("clusterDownloadType",
                                                         "FASTA is required for subsequent modules; CSV retains all features from data table"),
                                               
                                               # start button
                                               actionButton("clusterStart", label = h5("Start"),
                                                            style='padding:11px; font-size:80%'),
                                               
                                               # download button
                                               downloadButton("clusterDownload", label = h5("Download"),
                                                              style='padding:2px; font-size:80%'),
                                               
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
                                                         multiple = F, placeholder = "Clustered FASTA file",
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
                                                            choices = 3:5, selected = 3, inline = T),
                                               shinyBS::bsTooltip("kmerPCAButton_k",
                                                                  "Should the PCA be based on 3-, 4-, or 5-mers?"), 
                                               
                                               # slider for top number of clusters to plot
                                               sliderInput("kmerPCASlider_topClusters", label = h5("Number of top clusters to plot:"),
                                                           min = 1, max = 20, value = 10, step = 1),
                                               shinyBS::bsTooltip("kmerPCASlider_topClusters",
                                                                  "How many of the top clusters do you want to plot?"),
                                               
                                               # select value of k for kmer PCA
                                               radioButtons("kmerPCAButton_keepNC", label = h5("Keep non-clustered sequences?"),
                                                            choices = c("Yes", "No"), selected = "Yes", inline = T),
                                               shinyBS::bsTooltip("kmerPCAButton_keepNC",
                                                                  "If Yes, keep top N clusters AND the non-clustered sequences."), 
                                               
                                               # start button
                                               actionButton("kmerPCAStart", label = h5("k-mer PCA"), style='padding:11px; font-size:80%'),
                                               shinyBS::bsTooltip("kmerPCAStart", "Plot PCA of k-mer matrix, colored by cluster"),
                                               h5("*Alphabets outside of [A,C,G,T] are not yet supported."),
                                               
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
                                                         multiple = T, placeholder = "Cluster analysis files",
                                                         accept = c('.csv')),
                                               
                                               # note on selecting multiple files
                                               em("Holding ctrl (Windows) or command (Mac) will allow you to click multiple files."),
                                               tags$br(),
                                               
                                               # start button
                                               actionButton("clusterEnrichStart", label = h5("Start"),
                                                            style='padding:11px; font-size:80%'),
                                               
                                               # download button
                                               downloadButton("clusterEnrichDownload", label = h5("Download"),
                                                              style='padding:2px; font-size:80%')
                                           ),
                                           
                                           # display count output as datatable
                                           mainPanel(shinycssloaders::withSpinner(DT::dataTableOutput("clusterEnrichOutput")))
                                       ))
                          )
                 ),
                 
                 ## help panel; displays the README file
                 tabPanel("Help",
                          tags$iframe(style="height:800px; width:100%; scrolling=yes", src="UserGuide.pdf")
                 )
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
                shinyjs::html(id = "countTextOutput", html = m$message, add = F)
            })
        }
    })
    
    ## COUNT - DATA DISPLAY
    output$countOutput <- DT::renderDataTable(DT::datatable({
        countDF()
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
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
                write.table(fa_formatOutput(outputData = countDF()), file, quote = F, row.names = F, col.names = F)
            }else{
                write.csv(countDF(), file, quote = F, row.names = F)
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
    
    ## TRANSLATE - DATA GENERATION
    translateDF <- eventReactive(input$translateStart, {
        if(is.null(isolate(input$translateInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_translate(fastaInput = isolate(input$translateInput$datapath),
                         orf = isolate(input$orfButton),
                         converge = ifelse(isolate(input$convergeButton) == "Yes", T, F))
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
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
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
                write.table(fa_formatOutput(outputData = translateDF()), file, quote = F, row.names = F, col.names = F)
            }else{
                write.csv(translateDF(), file, quote = F, row.names = F)
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
        } else if(grepl("[^a-zA-Z]", gsub(",| ", "", isolate(input$motifInput_search)))){
            showNotification("Pattern list may not contains numerics or special characters other than commas!", type = "error",
                             duration = NULL)
            return(NULL)
        }else{
            ## function call
            fa_motifSearch(fastaInput = isolate(input$motifSearchInput$datapath),
                           motif = isolate(input$motifInput_search),
                           highlight = ifelse(isolate(input$motifSearchButton_highlight) == "Yes", T, F),
                           partial = ifelse(isolate(input$motifSearchButton_partial) == "Yes", T, F),
                           motifType = isolate(input$motifSearchButton_motifType))
        }
    })
    
    ## MOTIF SEARCH - DATA DISPLAY
    output$motifSearchOutput <- DT::renderDataTable(DT::datatable({
        motifSearchDF()
    },
    filter = list(position = "top", plain = T, clear = F), rownames = F,
    options = list(searchHighlight = T, search = list(regex = T,
                                                      search = fa_motif_format(motifList = isolate(input$motifInput_search),
                                                                               motifType = isolate(input$motifSearchButton_motifType))))
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
                write.table(fa_formatOutput(outputData = motifSearchDF()), file, quote = F, row.names = F, col.names = F)
            }else{
                write.csv(motifSearchDF(), file, quote = F, row.names = F)
            }
        }
    )

    ## MOTIF TRACKER - DATA GENERATION
    motifTrackerDF <- eventReactive(input$motifTrackerStart, {
        if(is.null(isolate(input$motifTrackerInput))){
            showNotification("No files provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$motifTrackerInput)) < 2 | nrow(isolate(input$motifTrackerInput)) > 3){
            showNotification("Supply 2-3 files!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$motifInput_enrich) == ""){
            showNotification("Must supply valid motif(s)!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z]", gsub(",| ", "", isolate(input$motifInput_enrich)))){
            showNotification("Motif list may not contains numerics or special characters other than commas!", type = "error",
                             duration = NULL)
            return(NULL)
        } else{
            fa_motifTracker(fastaInputs = isolate(input$motifTrackerInput$datapath),
                            motif = isolate(input$motifInput_enrich),
                            fileNames = isolate(input$motifTrackerInput$name),
                            motifType = isolate(input$motifTrackerButton_motifType))
        }
    })
    
    ## MOTIF TRACKER - DATA DISPLAY
    output$motifTrackerOutput <- DT::renderDataTable(DT::datatable({
        motifTrackerDF()
    }, rownames = F
    ))
    
    ## MOTIF TRACKER - METADATA
    output$motifTrackerUI <- renderUI({
        HTML(paste(fa_motifTracker_scores(RPM = motifTrackerDF()$Motif.RPM), collapse = "<br/>"))
    })
    
    ## MOTIF TRACKER - DATA DOWNLOAD
    output$motifTrackerDownload <- downloadHandler(
        # set filename
        filename = function(){
            "motifTracker.csv"
        },
        
        # set file content
        content = function(file){
            # format data for output as CSV file
            write.csv(motifTrackerDF(), file, row.names = F, quote = F)
            write.table("", file, append = T, row.names = F, col.names = F)
            write.table(HTML(paste(fa_motifTracker_scores(RPM = motifTrackerDF()$Motif.RPM), collapse = "<br/>")),
                        file, append = T, row.names = F, col.names = F)
        }
    )
    
    ## DISTANCE - DATA GENERATION
    distanceDF <- eventReactive(input$distanceStart, {
        if(is.null(isolate(input$distanceInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(isolate(input$querySequence) == ""){
            showNotification("Must supply valid query sequence!", type = "error", duration = NULL)
            return(NULL)
        } else if(grepl("[^a-zA-Z]", gsub(" ", "", isolate(input$querySequence)))){
            showNotification("Query sequence may not contains numerics or special characters!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_distance(dataInput = isolate(input$distanceInput$datapath),
                        querySequence = isolate(input$querySequence))
        }
    })
    
    ## DISTANCE - DATA OUTPUT
    output$distanceOutput <- DT::renderDataTable(DT::datatable({
        distanceDF()
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
    ))
    
    ## DISTANCE - DATA DOWNLOAD
    output$distanceDownload <- downloadHandler(
        # set filename
        filename = function(){
            fa_distance_outputName(inputFile = input$distanceInput$name)
        },
        
        # set file content
        content = function(file){
            write.csv(distanceDF(), file, row.names = F, quote = F)
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
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            distance_histogram()
        }
    })
    
    ## SEQUENCE ENRICHMENT - DATA GENERATION
    enrichDF <- eventReactive(input$enrichStart, {
        if(is.null(isolate(input$enrichInput))){
            showNotification("No files provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$enrichInput)) < 2 | nrow(isolate(input$enrichInput)) > 3){
            showNotification("Supply 2-3 files!", type = "error", duration = NULL)
            return(NULL)
        } else{
            fa_enrich(fastaInputs = isolate(input$enrichInput$datapath),
                      removeNA = ifelse(isolate(input$enrichButton) == "Yes", T, F))
        }
    })
    
    ## SEQUENCE ENRICHMENT - DATA OUTPUT
    output$enrichOutput <- DT::renderDataTable(DT::datatable({
        enrichDF()
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
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
            write.csv(enrichDF()[input[["enrichOutput_rows_all"]],], file, row.names = F, quote = F)
        }
    )
    
    ## SEQUENCE ENRICHMENT - FOLD-CHANGE HISTOGRAMs
    enrich_fcHistData <- eventReactive(input$fcHistStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else{
            fa_enrich_histogram(df = enrichDF()[input[["enrichOutput_rows_all"]],])
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY FOLD-CHANGE HISTOGRAMS
    output$fcHistOutput <- plotly::renderPlotly({
        if(is.null(enrich_fcHistData())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            enrich_fcHistData()
        }
    })
    
    ## SEQUENCE ENRICHMENT - RPM SCATTER PLOTS
    enrich_rpmScatterData <- eventReactive(input$rpmScatterStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else{
            fa_enrich_scatter(df = df <- enrichDF()[input[["enrichOutput_rows_all"]],])
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY RPM SCATTER PLOTS
    output$rpmScatterOutput <- plotly::renderPlotly({
        if(is.null(enrich_rpmScatterData())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            enrich_rpmScatterData()
        }
    })
    
    ## SEQUENCE ENRICHMENT - VOLCANO PLOTS
    enrich_volcanoData <- eventReactive(input$volcanoStart, {
        if(is.null(enrichDF())){
            return(NULL)
        } else{
            fa_enrich_volcano(df = enrichDF()[input[["enrichOutput_rows_all"]],])
        }
    })
    
    ## SEQUENCE ENRICHMENT - DISPLAY VOLCANO PLOTS
    output$volcanoOutput <- plotly::renderPlotly({
        if(is.null(enrich_volcanoData())){
            showNotification("Please generate a data table before plotting!", type = "error", duration = NULL)
            return(NULL)
        } else{
            enrich_volcanoData()
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
        } else{
            # capture output
            withCallingHandlers({
                shinyjs::html("clusterTextOutput", "")
                fa_clusterLED(fastaInput = isolate(input$clusterInput$datapath),
                              minReads = isolate(input$clusterSlider_minReads),
                              maxLED = isolate(input$clusterSlider_maxLED),
                              totalClusters = isolate(input$clusterSlider_totalClusters),
                              multipleOutputs = ifelse(isolate(input$clusterButton_outputs) == "Yes", T, F),
                              outputDirectory = isolate(input$clusterInput_directory),
                              keepNC = ifelse(isolate(input$clusterButton_keepNC) == "Yes", T, F))
            },
            # redirect output to text in UI
            message = function(m){
                shinyjs::html(id = "clusterTextOutput", html = m$message, add = F)
            })
        }
    })
    
    ## CLUSTER - DATA OUTPUT
    output$clusterOutput <- DT::renderDataTable(DT::datatable({
        clusterDF()
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
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
                write.table(fa_formatOutput(outputData = clusterDF()), file, quote = F, row.names = F, col.names = F)
            }else{
                write.csv(clusterDF(), file, quote = F, row.names = F)
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
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
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
            write.csv(clusterDiversityDF(), file, row.names = F, quote = F)
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
                                        keepNC = ifelse(isolate(input$kmerPCAButton_keepNC) == "Yes", T, F))
        }
    })
    
    ## CLUSTER DIVERSITY - K-MER PCA - PLOT
    output$kmerPCAOutput <- plotly::renderPlotly({
        kmerPCA()
    })
    
    ## CLUSTER ENRICH - DATA GENERATION
    clusterEnrichDF <- eventReactive(input$clusterEnrichStart, {
        if(is.null(isolate(input$clusterEnrichInput))){
            showNotification("No file provided!", type = "error", duration = NULL)
            return(NULL)
        } else if(nrow(isolate(input$clusterEnrichInput)) < 2 | nrow(isolate(input$clusterEnrichInput)) > 3){
            showNotification("Supply 2-3 files!", type = "error", duration = NULL)
            return(NULL)
        }else{
            fa_clusterEnrich(clusterCSVs = isolate(input$clusterEnrichInput$datapath))
        }
    })
    
    ## CLUSTER ENRICH - DATA OUTPUT
    output$clusterEnrichOutput <- DT::renderDataTable(DT::datatable({
        clusterEnrichDF()
    }, filter = list(position = "top", plain = T, clear = F), rownames = F
    ))
    
    ## CLUSTER ENRICH - DATA DOWNLOAD
    output$clusterEnrichDownload <- downloadHandler(
        # set filename
        filename = "clusterEnrich.csv",
        
        # set file content
        content = function(file){
            # format data for output as CSV file
            write.csv(clusterEnrichDF(), file, row.names = F, quote = F)
        }
    )
}
shinyApp(ui = ui, server = server)