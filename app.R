# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

## imports
library(shiny)
library(shinyBS)
library(shinyFiles)
library(colourpicker)

library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(plotly)
library(UpSetR)

## source files for support functions
source("./functions/functions_support.R")
source("./functions/functions_ui.R")

## source files for analytical functions
source("./functions/count_analyses.R")
source("./functions/translate_analyses.R")
source("./functions/motifSearch_analyses.R")
source("./functions/motifTracker_analyses.R")
source("./functions/motifDiscovery_analyses.R")
source("./functions/mutationNetwork_analyses.R")
source("./functions/distance_analyses.R")
source("./functions/dataMerge_analyses.R")
source("./functions/seqEnrich_analyses.R")
source("./functions/cluster_analyses.R")
source("./functions/clusterDiversity_analyses.R")
source("./functions/clusterEnrich_analyses.R")

## source files for visualization functions
source("./functions/count_plots.R")
source("./functions/motifTracker_plots.R")
source("./functions/motifDiscovery_plots.R")
source("./functions/distance_plots.R")
source("./functions/dataMerge_plots.R")
source("./functions/seqEnrich_plots.R")
source("./functions/posEnrich_plots.R")
source("./functions/clusterDiversity_plots.R")
source("./functions/clusterEnrich_plots.R")

## source files for tabs
source("./uiTabs/countTab.R")
source("./uiTabs/translateTab.R")
source("./uiTabs/motifTab.R")
source("./uiTabs/distanceTab.R")
source("./uiTabs/mutationNetworkTab.R")
source("./uiTabs/dataMergeTab.R")
source("./uiTabs/enrichTab.R")
source("./uiTabs/clusterTab.R")
source("./uiTabs/aboutTab.R")

## change limit for file sizes
options(shiny.maxRequestSize=2000*1024^2)

## sanitize error messages
# options(shiny.sanitize.errors = TRUE)

## define ui
ui <- navbarPage(
  "FASTAptameR 2.0",
  
  ## application theme
  theme = shinythemes::shinytheme("cosmo"),
  
  countTab,
  translateTab,
  motifTab,
  mutationNetworkTab,
  distanceTab,
  dataMergeTab,
  enrichTab,
  clusterTab,
  aboutTab,
  
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
    
    if(is.null(isolate(input$countInput))){
      showNotification("No file or link provided!", type = "error", duration = NULL)
      return(NULL)
    } else if(!(toupper(tools::file_ext(isolate(input$countInput$datapath))) %in% c("FQ", "FA", "FASTQ", "FASTA"))){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      # capture output
      withCallingHandlers({
        shinyjs::html("countTextOutput", "")
        fa_count(
          dataInput = isolate(input$countInput$datapath),
          reverseComplement = ifelse(isolate(input$reverseComplement) == "Yes", TRUE, FALSE)
        )
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
        write.table(
          fa_formatOutput(outputData = countDF()[input[["countOutput_rows_all"]],]),
          file,
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )
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
      fa_count_rpr(
        countData = countDF(),
        minReads = isolate(input$countSlider_minReads),
        maxRanks = isolate(input$countSlider_maxRanks)
      )
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
    } else if(tools::file_ext(isolate(input$translateInput$datapath)) != "fasta"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_translate(
        fastaInput = isolate(input$translateInput$datapath),
        orf = isolate(input$orfButton),
        converge = ifelse(isolate(input$convergeButton) == "Yes", TRUE, FALSE),
        inputChanges = isolate(input$translateInput_changes),
        translateSelection = isolate(input$translateSelection)
      )
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
        write.table(
          fa_formatOutput(outputData = translateDF()[input[["translateOutput_rows_all"]],]),
          file,
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )
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
      fa_count_rpr(
        countData = translateDF(),
        minReads = isolate(input$translateSlider_minReads),
        maxRanks = isolate(input$translateSlider_maxRanks)
      )
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
      showNotification("Pattern list must be alphanumeric!", type = "error", duration = NULL)
      return(NULL)
    } else if(tools::file_ext(isolate(input$motifSearchInput$datapath)) != "fasta"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      
      ## function call
      fa_motifSearch(
        fastaInput = isolate(input$motifSearchInput$datapath),
        motif = isolate(input$motifInput_search),
        highlight = ifelse(isolate(input$motifSearchButton_highlight) == "Yes", TRUE, FALSE),
        partial = ifelse(isolate(input$motifSearchButton_partial) == "Yes", TRUE, FALSE),
        motifType = isolate(input$motifSearchButton_motifType)
      )
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
        write.table(
          fa_formatOutput(outputData = motifSearchDF()[input[["motifSearchOutput_rows_all"]],]),
          file,
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )
      }else{
        write.csv(
          motifSearchDF()[input[["motifSearchOutput_rows_all"]],],
          file,
          quote = FALSE, row.names = FALSE
        )
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
    } else if(isolate(input$motifTracker_selectInput)[1] == "*"){
      showNotification("Asterisk is a placeholder and not a valid file ordering!", duration = NULL)
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
    } else if(
      ifelse(
        sum(
          ifelse(tools::file_ext(isolate(input$motifTrackerInput$datapath)) == "fasta", TRUE, FALSE)) != nrow(isolate(input$motifTrackerInput)
          ),
        TRUE,
        FALSE
      )
    ){
      showNotification("Please confirm that all files have the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      
      # check if user is supplying motifs or sequences
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
      write.csv(motifTrackerDF_enrich(), file, row.names = FALSE, quote = TRUE)
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
  
  ## MOTIF DISCOVERY - DATA GENERATION
  motifDiscoveryDF <- eventReactive(input$motifDiscovery_start, {
    if(is.null(isolate(input$motifDiscovery_input))){
      showNotification("No file provided!", type = "error", duration = NULL)
      return(NULL)
    } else if(tools::file_ext(isolate(input$motifDiscovery_input$datapath)) != "fasta"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_fsbc_motifDiscovery(
        fastaInput = isolate(input$motifDiscovery_input$datapath),
        minReads = isolate(input$motifDiscovery_minReads),
        lengthRange = isolate(input$motifDiscovery_lengthRange)
      )
    }
  })
  
  ## MOTIF DISCOVERY - BUBBLE PLOT - RENDER
  motifDiscovery_plot <- eventReactive(input$motifDiscovery_start, {
    if(is.null(motifDiscoveryDF())){
      return(NULL)
    } else{
      fa_fsbc_motifDiscoveryPlot(motifDF = motifDiscoveryDF())
    }
  })
  
  ## MOTIF DISCOVERY - BUBBLE PLOT - PLOT
  output$motifDiscovery_plot <- plotly::renderPlotly({
    if(is.null(motifDiscovery_plot())){
      return(NULL)
    } else{
      motifDiscovery_plot()
    }
  })
  
  ## MOTIF DISCOVERY - DATA OUTPUT
  output$motifDiscovery_output <- DT::renderDataTable(DT::datatable({
    motifDiscoveryDF()
  },
  filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
  ))
  
  ## MOTIF DISCOVERY DOWNLOAD
  output$motifDiscovery_download <- downloadHandler(
    # set filename
    filename = function(){
      "motifDiscovery.csv"
    },
    
    # set file content
    content = function(file){
      
      # format data for output as CSV file
      write.csv(motifDiscoveryDF(), file, row.names = FALSE, quote = TRUE)
    }
  )
  
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
    } else if(tools::file_ext(isolate(input$distanceInput$datapath)) != "fasta"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_distance(
        dataInput = isolate(input$distanceInput$datapath),
        querySequence = isolate(input$querySequence),
        seqRange = isolate(input$distanceTruncRange)
      )
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
        write.table(
          fa_formatOutput(outputData = distanceDF()[input[["distanceOutput_rows_all"]],]),
          file,
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )
      }else{
        write.csv(
          distanceDF()[input[["distanceOutput_rows_all"]],],
          file,
          quote = FALSE, row.names = FALSE
        )
      }
    }
  )
  
  ## DISTANCE - DISTANCE HISTOGRAM - RENDER
  distance_histogram <- eventReactive(input$distance_histStart, {
    if(is.null(distanceDF())){
      return(NULL)
    } else{
      fa_distance_histogram(distanceData = distanceDF(), querySequence = isolate(input$querySequence))
    }
  })
  
  ## DISTANCE - DISTANCE HISTOGRAM - PLOT
  output$distance_histOutput <- plotly::renderPlotly({
    if(is.null(distance_histogram())){
      showNotification("Please generate a data table from a single population before plotting!", type = "error", duration = NULL)
      return(NULL)
    } else{
      distance_histogram()
    }
  })
  
  # ## MUTATION NETWORK - DATA GENERATION
  # mutationNetworkDF <- eventReactive(input$mutationNetwork_start, {
  #   if(is.null(isolate(input$mutationNetwork_input))){
  #     showNotification("No file provided!", type = "error", duration = NULL)
  #     return(NULL)
  #   } else if(isolate(input$mutationNetwork_startNode) == ""){
  #     showNotification("Must supply valid start sequence!", type = "error", duration = NULL)
  #     return(NULL)
  #   } else if(grepl("[^a-zA-Z0-9]", gsub(" ", "", isolate(input$mutationNetwork_startNode)))){
  #     showNotification("Start sequence must be alphanumeric!", type = "error", duration = NULL)
  #     return(NULL)
  #   } else if(isolate(input$mutationNetwork_endNode) == ""){
  #     showNotification("Must supply valid end sequence!", type = "error", duration = NULL)
  #     return(NULL)
  #   } else if(grepl("[^a-zA-Z0-9]", gsub(" ", "", isolate(input$mutationNetwork_endNode)))){
  #     showNotification("End sequence must be alphanumeric!", type = "error", duration = NULL)
  #     return(NULL)
  #   } else{
  #     fa_mutationalIntermediates(
  #       fastaInput = isolate(input$mutationNetwork_input$datapath),
  #       startNode = isolate(input$mutationNetwork_startNode),
  #       endNode = isolate(input$mutationNetwork_endNode),
  #       maxCost = isolate(input$mutationNetwork_maxCost)
  #     )
  #   }
  # })
  # 
  # ## MUTATION NETWORK - DATA OUTPUT; this is a data.table if a path was found or a character string if not
  # output$mutationNetwork_DT_output <- DT::renderDataTable(DT::datatable({
  #   if(is.data.frame(mutationNetworkDF())){
  #     mutationNetworkDF()
  #   }
  # },
  # filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
  # ))
  # 
  # output$mutationNetwork_text_output <- renderUI(
  #   if(!is.data.frame(mutationNetworkDF())){
  #     mutationNetworkDF()
  #   }
  # )
  # 
  # ## MUTATION NETWORK DOWNLOAD
  # output$mutationNetwork_download <- downloadHandler(
  #   # set filename
  #   filename = function(){
  #     "mutationNetwork.csv"
  #   },
  #   
  #   # set file content
  #   content = function(file){
  #     
  #     # format data for output as CSV file
  #     write.csv(mutationNetworkDF(), file, row.names = FALSE, quote = TRUE)
  #   }
  # )
  
  ## DATA MERGE - UPDATE FILE SELECTIONS
  observe({
    updateSelectizeInput(session = session, inputId = "dataMerge_selectInput", choices = input$dataMerge_input$name)
  })
  
  ## DATA MERGE - DATA GENERATION
  dataMergeDF <- eventReactive(input$dataMerge_start, {
    if(is.null(isolate(input$dataMerge_input))){
      showNotification("No files provided!", type = "error", duration = NULL)
      return(NULL)
    } else if(nrow(isolate(input$dataMerge_input)) < 2){
      showNotification("Please supply at least 2 files!", type = "error", duration = NULL)
      return(NULL)
    } else if(length(isolate(input$dataMerge_selectInput)) < 2){
      showNotification("Please order all of your files!", type = "error", duration = NULL)
      return(NULL)
    } else if(isolate(input$dataMerge_selectInput)[1] == "*"){
      showNotification("Asterisk is a placeholder and not a valid file ordering!", duration = NULL)
      return(NULL)
    } else if(
      ifelse(
        sum(
          ifelse(tools::file_ext(isolate(input$dataMerge_input$datapath)) == "fasta", TRUE, FALSE)) != nrow(isolate(input$dataMerge_input)
          ),
        TRUE,
        FALSE
      )
    ){
      showNotification("Please confirm that all files have the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      
      fa_dataMerge(
        fastaInputs = isolate(input$dataMerge_input[match(input$dataMerge_selectInput, input$dataMerge_input$name),]$datapath),
        mergeType = isolate(input$dataMerge_mergeType)
      )
    }
  })
  
  ## DATA MERGE - DATA OUTPUT
  output$dataMerge_output <- DT::renderDataTable(DT::datatable({
    dataMergeDF()
  }, filter = list(position = "top", plain = TRUE, clear = FALSE), rownames = FALSE
  ))
  
  ## DATA MERGE - DATA DOWNLOAD
  output$dataMerge_download <- downloadHandler(
    # set filename
    filename = function(){
      "dataMerge.csv"
    },
    
    # set file content
    content = function(file){
      # format data for output as CSV file
      write.csv(dataMergeDF()[input[["dataMerge_output_rows_all"]],], file, row.names = FALSE, quote = FALSE)
    }
  )
  
  ## DATA MERGE - GENERATE SEQUENCE PERSISTENCE PLOT
  dataMerge_sequencePersistencePlot <- eventReactive(input$dataMerge_start, {
    if(is.null(isolate(input$dataMerge_input))){
      return(NULL)
    } else if(nrow(isolate(input$dataMerge_input)) < 2){
      return(NULL)
    } else if(length(isolate(input$dataMerge_selectInput)) < 2){
      return(NULL)
    } else if(isolate(input$dataMerge_selectInput)[1] == "*"){
      return(NULL)
    } else{
      fa_dataMerge_seqPersistence(fastaInputs = isolate(input$dataMerge_input[match(input$dataMerge_selectInput, input$dataMerge_input$name),]$datapath))
    }
  })
  
  ## DATA MERGE - RENDER SEQUENCE PERSISENCE PLOT
  output$dataMerge_sequencePersistence <- plotly::renderPlotly({
    if(is.null(dataMerge_sequencePersistencePlot())){
      return(NULL)
    } else{
      dataMerge_sequencePersistencePlot()
    }
  })
  
  ## DATA MERGE - GENERATE UpSetR PLOT
  dataMerge_UpSetRPlot <- eventReactive(input$dataMerge_start, {
    if(is.null(isolate(input$dataMerge_input))){
      return(NULL)
    } else if(nrow(isolate(input$dataMerge_input)) < 2){
      return(NULL)
    } else if(length(isolate(input$dataMerge_selectInput)) < 2){
      return(NULL)
    } else if(isolate(input$dataMerge_selectInput)[1] == "*"){
      return(NULL)
    } else{
      fa_dataMerge_UpSetR(fastaInputs = isolate(input$dataMerge_input$datapath), fastaNames = isolate(input$dataMerge_input$name))
    }
  })
  
  ## DATA MERGE - RENDER UpSetR PLOT
  output$dataMerge_UpSetR <- renderPlot({
    if(is.null(dataMerge_UpSetRPlot())){
      return(NULL)
    } else{
      dataMerge_UpSetRPlot()
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
    } else if(nrow(isolate(input$enrichInput)) != 2){
      showNotification("Supply exactly 2 files!", type = "error", duration = NULL)
      return(NULL)
    } else if(length(isolate(input$enrich_selectInput)) != 2){
      showNotification("Please order your files!", type = "error", duration = NULL)
      return(NULL)
    } else if(isolate(input$enrich_selectInput)[1] == "*"){
      showNotification("Asterisk is a placeholder and not a valid file ordering!", duration = NULL)
      return(NULL)
    } else if(tools::file_ext(isolate(input$enrichInput$datapath[1])) != "fasta"){
      showNotification("File 1 does not have the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else if(tools::file_ext(isolate(input$enrichInput$datapath[2])) != "fasta"){
      showNotification("File 2 does not have the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_enrich(
        fastaInputs = isolate(input$enrichInput[match(input$enrich_selectInput, input$enrichInput$name),]$datapath),
        keepNA = ifelse(isolate(input$enrichKeepNA) == "Yes", TRUE, FALSE)
      )
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
  
  ## SEQUENCE ENRICHMENT - GENERATE LOG2(ENRICHMENT) HISTOGRAM
  enrich_histogram <- eventReactive(input$fcHistStart, {
    if(is.null(enrichDF())){
      showNotification("Generate an enrichment table!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_enrich_histogram(df = enrichDF())
    }
  })
  
  ## SEQUENCE ENRICHMENT - RENDER LOG2(ENRICHMENT) HISTOGRAM
  output$fcHistOutput <- plotly::renderPlotly({
    if(is.null(enrich_histogram())){
      return(NULL)
    } else{
      enrich_histogram()
    }
  })
  
  ## SEQUENCE ENRICHMENT - GENERATE RPM SCATTER PLOT
  enrich_rpmScatterPlot <- eventReactive(input$rpmScatterStart, {
    if(is.null(enrichDF())){
      showNotification("Generate an enrichment table!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_enrich_scatter(df = enrichDF())
    }
  })
  
  ## SEQUENCE ENRICHMENT - RENDER RPM SCATTER PLOT
  output$rpmScatterOutput <- plotly::renderPlotly({
    if(is.null(enrich_rpmScatterPlot())){
      return(NULL)
    } else{
      enrich_rpmScatterPlot()
    }
  })
  
  ## SEQUENCE ENRICHMENT - GENERATE RA PLOT
  enrich_raPlot <- eventReactive(input$raStart, {
    if(is.null(enrichDF())){
      showNotification("Generate an enrichment table!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_enrich_ra(df = enrichDF())
    }
  })
  
  ## SEQUENCE ENRICHMENT - RENDER RA PLOT
  output$raOutput <- plotly::renderPlotly({
    if(is.null(enrich_raPlot())){
      return(NULL)
    } else{
      enrich_raPlot()
    }
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
    } else if(tools::file_ext(isolate(input$posEnrichInput$datapath)) != "csv"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
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
    } else if(tools::file_ext(isolate(input$posEnrichInput$datapath)) != "csv"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_enrich_heatMap(
        dataPath = isolate(input$posEnrichInput$datapath),
        refSeq = isolate(input$posEnrich_refSequence),
        modList = isolate(input$posEnrich_mods),
        enrichRange = isolate(input$slider_enrichmentRange),
        seqType = isolate(input$posEnrich_seqType),
        lowCol = isolate(input$lowCol), midCol = isolate(input$midCol), highCol = isolate(input$highCol)
      )
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
    } else if(tools::file_ext(isolate(input$clusterInput$datapath)) != "fasta"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      
      # capture output
      withCallingHandlers({
        shinyjs::html("clusterTextOutput", "")
        fa_clusterLED(
          fastaInput = isolate(input$clusterInput$datapath),
          minReads = isolate(input$clusterSlider_minReads),
          maxLED = isolate(input$clusterSlider_maxLED),
          totalClusters = isolate(input$clusterSlider_totalClusters),
          multipleOutputs = ifelse(isolate(input$clusterButton_outputs) == "Yes", TRUE, FALSE),
          outputDirectory = isolate(input$clusterInput_directory),
          keepNC = ifelse(isolate(input$clusterButton_keepNC) == "Yes", TRUE, FALSE)
        )
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
        write.table(
          fa_formatOutput(outputData = clusterDF()[input[["clusterOutput_rows_all"]],]),
          file,
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )
      }else{
        write.csv(
          clusterDF()[input[["clusterOutput_rows_all"]],],
          file, 
          quote = FALSE, row.names = FALSE
        )
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
    } else if(tools::file_ext(isolate(input$clusterDiversityInput$datapath)) != "fasta"){
      showNotification("Please supply a file with the correct extension!", type = "error", duration = NULL)
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
  
  ## CLUSTER DIVERSITY - UPDATE CLUSTER SELECTIONS
  observe({
    updateSelectizeInput(session = session, inputId = "kmerPCA_clusters", choices = clusterDiversityDF()$Cluster)
  })
  
  ## CLUSTER DIVERSITY - K-MER PCA - RENDER
  kmerPCA <- eventReactive(input$kmerPCAStart, {
    if(is.null(isolate(input$clusterDiversityInput))){
      return(NULL)
    } else if(lengths(strsplit(readLines(isolate(input$clusterDiversityInput$datapath), n = 1), split = "-")) != 6){
      showNotification("Please provide a file from FASTAptameR-Cluster!", type = "error", duration = NULL)
      return(NULL)
    } else if(isolate(input$kmerPCA_clusters[1]) == "*" | length(isolate(input$kmerPCA_clusters)) < 1 | length(isolate(input$kmerPCA_clusters)) > 15){
      showNotification("Please select 1-15 clusters to visualize!", type = "error", duration = NULL)
      return(NULL)
    } else if(is.null(clusterDiversityDF())){
      showNotification("Please generate a table before generating this plot!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_clusterDiversity_kmerPCA(
        clusterFile = isolate(input$clusterDiversityInput$datapath),
        kmerSize = as.numeric(isolate(input$kmerPCAButton_k)),
        clustersToPlot = isolate(input$kmerPCA_clusters)
      )
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
      showNotification("Please order your files!", type = "error", duration = NULL)
      return(NULL)
    } else if(length(isolate(input$clusterEnrich_selectInput)) == 1 & isolate(input$clusterEnrich_selectInput[1]) == "*"){
      showNotification("Asterisk is a placeholder and not a valid file ordering!", duration = NULL)
      return(NULL)
    } else if(
      ifelse(
        sum(
          ifelse(tools::file_ext(isolate(input$clusterEnrichInput$datapath)) == "csv", TRUE, FALSE)) != nrow(isolate(input$clusterEnrichInput)
          ),
        TRUE,
        FALSE
      )
    ){
      showNotification("Please confirm that all files have the correct extension!", type = "error", duration = NULL)
      return(NULL)
    } else{
      fa_clusterEnrich(
        clusterCSVs = isolate(input$clusterEnrichInput[match(input$clusterEnrich_selectInput, input$clusterEnrichInput$name),]$datapath),
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
      write.csv(
        clusterEnrichOutput_enrichDF()[input[["clusterEnrichOutput_enrichDF_rows_all"]],],
        file,
        row.names = FALSE, quote = FALSE
      )
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
  
  session$onSessionEnded(stopApp)
}

shinyApp(ui = ui, server = server)