function(input, output, session) {
  options(shiny.maxRequestSize=1000*1024^2)
  
  consoleMessage <- reactiveVal("")
  output$consoleOut <- renderText({
    consoleMessage()
  })
  
  pdbInfo <- reactiveVal(FALSE)
  
  output$svmThreshSliders <- renderUI({
    if (input$separateFDRs) {
      tagList(
        sliderInput("svmThresholdIntra", "IntraProtein", min = -5, max = 10,
                    step = 0.1, value = -5),
        sliderInput("svmThresholdInter", "InterProtein", min = -5, max = 10,
                    step = 0.1, value = 0)
      )
    } else {
      sliderInput("svmThreshold", "SVM Score", min = -5, max = 10,
                  step = 0.1, value = 0)
    }
  })

  shinyFileChoose(input, "clmsData", roots=exDir, filetypes=c('', 'txt'))
  shinyFileChoose(input, "modules", roots=exDir, filetypes=c('', 'txt'))
  shinyFileChoose(input, "ms2pkls", roots=exDir, filetypes=c('', 'txt', 'mgf'))
  shinyFileChoose(input, "ms3pkls", roots=exDir, filetypes=c('', 'txt', 'mgf'))

  output$saveClassified <- downloadHandler(
    filename = function() {
      str_c(parseFilePaths(exDir, input$clmsData)$name)
    },
    content = function(file) {
      write_tsv(xlTable(), file)
    }
  )
 
  output$clmsDataFileName <- renderPrint({
    if (is.integer(input$clmsData)) {
      cat("No file selected")
    } else {
      cat(parseFilePaths(exDir, input$clmsData)$name)
    }
  })
  
  output$moduleFileName <- renderPrint({
    if (is.integer(input$modules)) {
      cat("No file selected")
    } else {
      cat(parseFilePaths(exDir, input$modules)$name)
      consoleMessage("*** Parsing Module File ***")
    }
  })
  
  output$ms2pklsName <- renderText({
    if (is.integer(input$ms2pkls)) {
      return("No file selected")
    } else {
      fnames <- parseFilePaths(exDir, input$ms2pkls)$name
      fnames <- str_replace(fnames, "\\.[[a-zA-Z0-9]]+?$", "")
      str_c(fnames, collapse="\n")
    }
  })
  
  output$ms3pklsName <- renderText({
    if (is.integer(input$ms3pkls)) {
      return("No file selected")
    } else {
      fnames <- parseFilePaths(exDir, input$ms3pkls)$name
      fnames <- str_replace(fnames, "\\.[[a-zA-Z0-9]]+?$", "")
      str_c(fnames, collapse="\n")
    }
  })
  
  scReader <- reactive({
    if (!is.integer(input$clmsData)) {
      if (input$experimentType == "ms3") {
        if (!is.integer(input$ms2pkls) & !is.integer(input$ms3pkls)) {
        inFile <- parseFilePaths(exDir, input$clmsData)
        ms3Files <- parseFilePaths(exDir, input$ms3pkls)$datapath
        ms2Files <- parseFilePaths(exDir, input$ms2pkls)$datapath
        masterScanFile <- createMasterScanFile(ms3Files, ms2Files)
        searchCompareMS3 <- readMS3results(inFile$datapath)
        searchCompareMS3 <- addMasterScanInfo(searchCompareMS3, masterScanFile)
        scTable <- processMS3xlinkResults(searchCompareMS3)
        scTable <- scTable %>% mutate(SVM.score = Score.Diff / 10)
        consoleMessage("*** Loading Search Compare File ***")
        return(scTable)
        }
      }
      else {
      inFile <- parseFilePaths(exDir, input$clmsData)
      consoleMessage("*** Loading Search Compare File ***")
      scTable <- readProspectorXLOutput(inFile$datapath)
      consoleMessage("*** Building SVM Classifier ***")
      scTable <- buildClassifier(scTable)
      return(scTable)
      }
    }
  })
  
  scResults <- reactiveVal()
  observe({scResults(scReader())})
  
  moduleFile <- reactive({
    if (!is.integer(input$modules)) {
      inFile <- parseFilePaths(exDir, input$modules)
      return(inFile$datapath)
    }
  })
  
  observeEvent(req(moduleFile()), {
    consoleMessage("*** Assigning Modules ***")
    msvFilePath <- parseFilePaths(exDir, input$clmsData)$datapath
    projFolder <- dirname(dirname(dirname(dirname(msvFilePath))))
    pdbFolder <- file.path(projFolder, "tstone", "pdb")
    scResults(processModuleFile(scResults(), moduleFile(), pdbFileDir = pdbFolder))
    if ("PDB" %in% names(scResults()) & sum(!is.na(scResults()$PDB)) > 0) {
      pdbInfo(TRUE)
    }
  })

  csmTab <- reactive({
    req(scResults())
    datTab <- scResults()
    # The links should be generated in touchStone upon reading the SC file:
    msvFilePath <- parseFilePaths(exDir, input$clmsData)$datapath
    msvFiles <- list.dirs(dirname(msvFilePath), recursive = F)
    btName <- dirname(dirname(msvFilePath))
    btNameDir <- dirname(btName)
    btParamFile <- dir(btNameDir, str_c(basename(btName), ".xml"))
    if (length(btParamFile) > 0) {
      btParams <- readParamsFile(file.path(btNameDir, btParamFile))
    } else {
      btParams <- NA
    }
    if (!is.na(btParams)) {
      instrumentType <- btParams %>%
        xml_find_all("instrument_name") %>%
        xml_text()
    } else {
      instrumentType <- NA
    }
    if (input$experimentType == "ms3") {
      ms3Files <- parseFilePaths(exDir, input$ms3pkls)$datapath
      ms2Files <- parseFilePaths(exDir, input$ms2pkls)$datapath
      ms3Files <- map_chr(ms3Files, function(x) str_replace(basename(x), ".mgf" ,""))
      names(ms2Files) <- map_chr(ms3Files, function(x) str_replace(basename(x), "\\.(txt|mgf)$" , ""))
      datTab <- datTab %>%
        mutate(Fraction.ms2 = map_chr(Fraction, function(x) ms2Files[x]),
               specMS3.1 = pmap_chr(list(msvFiles, Fraction, z.1, Peptide.1, 
                                         MSMS.Info.1), generateMSViewerLink.ms3),
               specMS3.2 = pmap_chr(list(msvFiles, Fraction, z.2, Peptide.2, 
                                         MSMS.Info.2), generateMSViewerLink.ms3),
               specMS2 = pmap_chr(list(dirname(Fraction.ms2), basename(Fraction.ms2), z, Peptide.1, Peptide.2,
                                       MSMS.Info, linkType="DSSOms3"), generateMSViewerLink))
    } else {
      datTab <- datTab %>%
        mutate(specMS2 = pmap_chr(list(msvFiles, Fraction, z, Peptide.1, Peptide.2,
                                    MSMS.Info, instrumentType), generateMSViewerLink))
    }
    datTab.len <- nrow(datTab)
    datTab <- datTab %>% mutate(id = str_c("xl", 1:nrow(datTab)),
                                keep = map_chr(id, function(x) {
                                  as.character(checkboxInput(x, label=NULL, value=T, width='20px'))
                                  }))
    minPPM = mmin(min(datTab$ppm, na.rm=T))
    maxPPM = mmax(max(datTab$ppm, na.rm=T))
    updateSliderInput(session, "ms1MassError", value = c(minPPM, maxPPM),
                      min = minPPM, max = maxPPM)
    minSVM = mmin(min(datTab$SVM.score, na.rm=T), 1)
    maxSVM = mmax(max(datTab$SVM.score, na.rm=T), 1)
    if (input$separateFDRs) {
      updateSliderInput(session, "svmThresholdInter", value = 0,
                       min = minSVM, max = maxSVM)
      updateSliderInput(session, "svmThresholdIntra", value = -5,
                       min = minSVM, max = maxSVM)
    } else {
      updateSliderInput(session, "svmThreshold", value = 0,
                        min = minSVM, max = maxSVM)
    }
    return(datTab)
  })
  
  clTab <- reactive({
    consoleMessage("*** Calculating Residue-Pairs ***")
    bestResPair(csmTab())
  })
  
  plTab <- reactive({
    consoleMessage("*** Calculating Protein-Pairs ***")
    bestProtPair(csmTab())
  })
  
  modTab <- reactive({
    req(moduleFile())
    consoleMessage("*** Calculating Module-Pairs ***")
    bestModPair(csmTab())
  })
  
  tabLevel <- reactive({
    req(csmTab())
    switch(input$summaryLevel,
           "CSMs" = csmTab(),
           "Unique Residue Pairs" = clTab(),
           "Protein Pairs" = plTab(),
           "Module Pairs" = modTab()
    )
  })
  
  tabLevelFiltered <- reactive({
    req(csmTab())
    tabLevel() %>% filter(Score.Diff >= input$scoreDiffThreshold,
                          between(Len.Pep.1, input$peptideLengthFilter[1], input$peptideLengthFilter[2]),
                          between(Len.Pep.2, input$peptideLengthFilter[1], input$peptideLengthFilter[2]),
                          between(ppm, input$ms1MassError[1], input$ms1MassError[2]))
  })
  
  xlTable <- reactive({
    req(csmTab())
    if (input$separateFDRs) {
      req(input$svmThresholdIntra, input$svmThresholdInter)
      tabLevelFiltered() %>% 
        filter(Decoy=="Target",
               (str_detect(xlinkClass, "intraProtein") & SVM.score >= input$svmThresholdIntra) |
                 (str_detect(xlinkClass, "interProtein") & SVM.score >= input$svmThresholdInter))
    } else {
      req(input$svmThreshold)
      tabLevelFiltered() %>% 
        filter(Decoy=="Target", SVM.score >= input$svmThreshold)
    }
  })
  
  xlTableDecoy <- reactive({
    req(csmTab())
    if (input$separateFDRs) {
      req(input$svmThresholdIntra, input$svmThresholdInter)
      tabLevelFiltered() %>% 
        filter(Decoy!="Target",
               (str_detect(xlinkClass, "intraProtein") & SVM.score >= input$svmThresholdIntra) |
                 (str_detect(xlinkClass, "interProtein") & SVM.score >= input$svmThresholdInter))
    } else {
      req(input$svmThreshold)
      tabLevelFiltered() %>% 
        filter(Decoy!="Target", SVM.score >= input$svmThreshold)
    }
  })
  
  classRatio <- reactiveVal()
  output$numberClassedLinks <- renderText({
    req(tabLevel())
    numIntra <- sum(str_count(xlTable()$xlinkClass, "intraProtein"))
    numInter <- sum(str_count(xlTable()$xlinkClass, "interProtein"))
    classRatio(round(numInter / (numIntra + numInter), 2))
    str_c("Number of ", input$summaryLevel, ": ", nrow(xlTable()),
    "; intraProtein: ", numIntra, "; interProtein: ", numInter)
  })
  
  prepDT <- function(xltab) {
    displayTable <- formatXLTable(xltab)
    wideCols <- which(names(displayTable) %in% c("xlinkeResPair",
                                                 "DB.Peptide.1",
                                                 "DB.Peptide.2",
                                                 "Protein.1",
                                                 "Protein.2",
                                                 "Peptide.1",
                                                 "Peptide.2",
                                                 "Module.1",
                                                 "Module.2",
                                                 "Fraction")) - 1
    DT::datatable(displayTable,
                  rownames=FALSE,
                  filter="top",
                  escape=FALSE,
                  options = list(autoWidth=TRUE,
                                 deferRender=TRUE,
                                 processing=TRUE,
                                 columnDefs=list(
                                   list(width = '250px', targets = as.list(wideCols))
                                   ),
                                 scrollX=TRUE,
                                 scrollY="80vh",
                                 scrollCollapse=TRUE,
                                 paging=TRUE,
                                 pageLength=100,
                                 search.caseInsensitive=TRUE,
                                 scroller=TRUE,
                                 preDrawCallback = JS(
                                   'function() { 
                                   Shiny.unbindAll(this.api().table().node()); 
                                   }'),
                                 drawCallback = JS(
                                   'function() { 
                                   Shiny.bindAll(this.api().table().node()); 
                                   } ')
                  ),
    )
  }
  
  output$dataFile <- DT::renderDataTable({
    req(xlTable())
    prepDT(xlTable())
  })

  output$dataFileSelected <- DT::renderDataTable({
    req(xlTableSelected())
    prepDT(xlTableSelected())
  })

  output$dataFileDecoy <- DT::renderDataTable({
    req(xlTableDecoy())
    prepDT(xlTableDecoy())
  })
  
  fdr <- reactiveVal()
  observe({
    if (input$separateFDRs) {
      req(input$svmThresholdIntra, input$svmThresholdInter)
      fdr(calculateSeparateFDRs(tabLevelFiltered(),
                                thresholdIntra = input$svmThresholdIntra,
                                thresholdInter = input$svmThresholdInter,
                                scalingFactor = input$scalingFactor)
      )
    } else {
      req(input$svmThreshold)
      fdr(calculateFDR(tabLevelFiltered(),
                       threshold = input$svmThreshold,
                       scalingFactor = input$scalingFactor)
      )
    }
  })
  
  output$FDR <- renderText({
    req(tabLevel())
    str_c("FDR: ", as.character(round(100 * fdr(), 2)), "%")
  })

  output$FDRplot <- renderPlot({
    req(tabLevelFiltered())
    if(input$separateFDRs) {
      par(mfrow=c(2, 1),
          mar=c(2.1, 4.1, 4.1, 2.1),
          oma=c(3.1, 0, 0, 0))
      plotMin <- mmin(min(tabLevelFiltered()$SVM.score), 1)
      plotMax <- mmax(max(tabLevelFiltered()$SVM.score), 1)
      fdrPlots(filter(tabLevelFiltered(), str_detect(xlinkClass, "intraProtein")), 
               scalingFactor = input$scalingFactor, cutoff=input$svmThresholdIntra,
               minValue = plotMin, maxValue = plotMax, title = "intraProtein FDR", 
               xlabel=NA, addLegend=F)
      fdrPlots(filter(tabLevelFiltered(), str_detect(xlinkClass, "interProtein")), 
               scalingFactor = input$scalingFactor, cutoff=input$svmThresholdInter,
               minValue = plotMin, maxValue = plotMax, title = "interProtein FDR")
    } else {
      fdrPlots(tabLevelFiltered(), 
               cutoff = input$svmThreshold,
               scalingFactor = input$scalingFactor, 
               title = "FDR plot")
    }
  })
  
  numHits <- reactiveVal()
  observe({
    if(input$separateFDRs) {
      numHits(generateErrorTableSeparate(tabLevelFiltered(), 
                                         scalingFactor = input$scalingFactor))
    } else {
      numHits(generateErrorTable(tabLevelFiltered(),
                                 scalingFactor = input$scalingFactor)
      )
    }
  })
  
  observeEvent(input$findThreshold, {
    req(tabLevel())
    if (input$separateFDRs) {
      threshold <- findSeparateThresholds(tabLevelFiltered(),
                                          targetER = input$targetFDR / 100,
                                          scalingFactor = input$scalingFactor)
      updateSliderInput(session, "svmThresholdIntra", value = threshold$intraThresh)
      updateSliderInput(session, "svmThresholdInter", value = threshold$interThresh)
    } else {
      threshold <- findThreshold(tabLevelFiltered(),
                                 targetER = input$targetFDR / 100,
                                 scalingFactor = input$scalingFactor)
      updateSliderInput(session, "svmThreshold", value = threshold$globalThresh)
    }
  })
  
  output$IIratio <- renderText({
    req(tabLevel())
    str_c("percent interProtein: ", round(100*classRatio(),1), "%")
  })
  
  output$thresholdPlot <- renderPlot({
    req(numHits())
    numHitsPlot(numHits(), fdr())
  })
  
  classedMassErrors <- reactive({
    subset(xlTable())$ppm
  })

  VR <- reactive({
    req(moduleFile(), pdbInfo())
    targetDists <- xlTable() %>% 
      filter(!is.na(distance)) %>% 
      pull(distance)
    sum(targetDists > input$distanceThreshold) / length(targetDists)
  })

  output$VR <- renderText({
    req(moduleFile, pdbInfo())
    str_c("Violation Rate: ", as.character(round(100 * VR(), 2)), "%"
    )
  })

  output$distancePlot <- renderPlot({
    req(moduleFile(), pdbInfo())
    distancePlot2(xlTable(), threshold = input$distanceThreshold)
  })

  output$massErrorPlot <- renderPlot({
    req(tabLevel())
    massErrorPlot(classedMassErrors(), 
                  lowPlotRange = floor(min(tabLevel()$ppm, na.rm=T)),
                  highPlotRange = ceiling(max(tabLevel()$ppm, na.rm=T)),
                  lowThresh = input$ms1MassError[1], 
                  highThresh = input$ms1MassError[2])
  })
  
  output$meanError <- renderText({
    req(tabLevel())
    str_c("mean error: ", round(mean(classedMassErrors(), na.rm=T),2), " ppm")
  })

  output$distanceSlider <- renderUI({
    req(moduleFile())
    sliderInput("distanceThreshold", "Violation Dist.", min = 0, max = 100,
                step = 1, value = 30) 
  })

  summaryProtData <- reactive({
    req(csmTab())
    summarizeProtData(xlTable())
  })
  
  summaryModulData <- reactive({
    req(moduleFile())
    summarizeModuleData(xlTable())
  })
  
  output$proteinPlot <- renderPlot({
    req(summaryProtData())
    summaryProtData() %>%
      ggplot(aes(Acc.1, Acc.2, size=counts, col=counts)) + 
      geom_point(na.rm=T) +
      scale_size_area(max_size = 10)+#range = c(-1, 10)) +
      scale_color_viridis_c(option="D") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, hjust=1),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            plot.title = element_text(face = "bold", size = 14)) +
      ggtitle("Protein Pairs")
  })

  xlTableSelected <- reactiveVal()
    protFilter <- reactiveVal()
    observeEvent(req(input$plot_click), {
      if (str_detect(input$plot_click$mapping$x, "Acc")) {
        protFilter(nearPoints(summaryProtData(), input$plot_click, threshold = 10) %>%
                     select(starts_with("Acc")) %>% 
                     slice(1) %>%
                     as.character
        )
        xlTableSelected(
          xlTable() %>% filter((Acc.1 == protFilter()[1] & Acc.2 == protFilter()[2]) |
                                 (Acc.1 == protFilter()[2] & Acc.2 == protFilter()[1]))
        )
      } else if (str_detect(input$plot_click$mapping$x, "Modul")) {
        protFilter(nearPoints(summaryModulData(), input$plot_click, threshold = 10) %>%
                     select(starts_with("Modul")) %>%
                     slice(1) %>%
                     as.character
        )
        xlTableSelected(
          xlTable() %>% filter((Module.1 == protFilter()[1] & Module.2 == protFilter()[2]) |
                                 (Module.1 == protFilter()[2] & Module.2 == protFilter()[1]))
        )
      }
      updateTabsetPanel(session, "navbar", selected = "Selected Crosslinks")
    })

  output$protHover <- renderText({
    req(xlTable())
    nearPoints(summaryProtData(), input$prot_hover) %>%
      select(starts_with("Acc")) %>%
      slice(1) %>%
      as.character
  })
  
  output$modHover <- renderText({
    req(xlTable())
    nearPoints(summaryModulData(), input$mod_hover) %>%
      select(starts_with("Modul")) %>%
      slice(1) %>%
      as.character
  })
  
  output$modulePlot <- renderPlot({
    req(summaryModulData)
    summaryModulData() %>% 
      ggplot(aes(Module.1, Module.2, size=counts, col=counts)) + 
      geom_point(na.rm=T) +
      scale_size_area(max_size = 10) + #range=c(-1,8)) +
      scale_color_viridis_c(option="D") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, hjust=1),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            plot.title = element_text(face = "bold", size = 14)) +
      ggtitle("Module Pairs")
  })

  observeEvent(input$viewXiNet, {
    req(xlTable())
    xiFile <- makeXiNetFile(xlTable())
    qs <- as.character(Sys.time()) %>% str_replace_all("[[\\s\\-\\:]]","")
    xiFileName <- str_c("xiFile", qs)
    xiFilePath <- str_c(pathToXiFile, "/", xiFileName, ".csv")
    write_csv(xiFile, xiFilePath)
    output$ui_open_tab <- renderUI({
      baseLink <- linkToXiView
      link <- str_c(baseLink, xiFileName, sep="?fileName=")
      tags$script(paste0("window.open('", link, "', '_blank')"))
    })
  })

  observeEvent(input$viewXiNetTwo, {
    req(xlTable())
    xiFile <- makeXiNetFile(xlTable())
    qs <- as.character(Sys.time()) %>% str_replace_all("[[\\s\\-\\:]]","")
    xiFileName <- str_c("xiFile", qs)
    xiFilePath <- str_c(pathToXiFile, "/", xiFileName, ".csv")
    write_csv(xiFile, xiFilePath)
    output$ui_open_tab_sel_two <- renderUI({
      baseLink <- linkToXiView
      link <- str_c(baseLink, xiFileName, sep="?fileName=")
      tags$script(paste0("window.open('", link, "', '_blank')"))
    })
  })
  
  observeEvent(input$viewXiNetSel, {
    req(xlTableSelected())
    xiFile <- makeXiNetFile(xlTableSelected())
    qs <- as.character(Sys.time()) %>% str_replace_all("[[\\s\\-\\:]]","")
    xiFileName <- str_c("xiFile", qs)
    xiFilePath <- str_c(pathToXiFile, "/", xiFileName, ".csv")
    write_csv(xiFile, xiFilePath)
    output$ui_open_tab_sel <- renderUI({
      baseLink <- linkToXiView
      link <- str_c(baseLink, xiFileName, sep="?fileName=")
      tags$script(paste0("window.open('", link, "', '_blank')"))
    })
  })
  
  observeEvent(input$scrapeMSP, {
    req(scResults())
    datTab <- scResults()
    require(future)
    require(furrr)
    future::plan(multicore)
    # The links should be generated in touchStone upon reading the SC file:
    msvFilePath <- parseFilePaths(exDir, input$clmsData)$datapath
    msvFiles <- list.dirs(dirname(msvFilePath), recursive = F)
    btName <- dirname(dirname(msvFilePath))
    btNameDir <- dirname(btName)
    btParamFile <- dir(btNameDir, str_c(basename(btName), ".xml"))
    if (length(btParamFile) > 0) {
      btParams <- readParamsFile(file.path(btNameDir, btParamFile))
    } else {
      btParams <- NA
    }
    if (!is.na(btParams)) {
      instrumentType <- btParams %>%
        xml_find_all("instrument_name") %>%
        xml_text()
    } else {
      instrumentType <- NA
    }
    if (input$experimentType == "ms3") {
    } else {
#      withProgress(message = "scraping MS-Product", value = 0, {
#        numPoints <- nrow(datTab)
        ms.product.info <-
          pmap_chr(list(msvFiles, datTab$Fraction, datTab$z, datTab$Peptide.1, 
                        datTab$Peptide.2, datTab$MSMS.Info, instrumentType,
                        outputType="Tab delimited text"), generateMSViewerLink) %>%
          future_map(function(msvLink) {
            spec.html <- xml2::read_html(msvLink)
            spec.node <- rvest::html_node(spec.html, xpath = '//*[@id="centerbody"]')
            spec.table <- read_tsv(rvest::html_text(spec.node))
#            incProgress(1/numPoints)
            return(spec.table)
          })
      }#)}
    percentMatched <- getPercentMatched(ms.product.info)
    datTab <- cbind(datTab, percentMatched)
    datTab <- buildClassifier(datTab, params.best, "SVM.new")
    scResults(datTab)
  })

  # These are test function to explore how to make the checkboxes work right.
  
  # output$test = renderPrint({
  #   req(xlTable())
  #   id <- xlTable()  %>% slice(1:100) %>% pull(id)
  #   values <- map_lgl(id, function(x) input[[x]])
  #   return(data.frame(id, values))
  # })
    
  # observeEvent(input$removeSel, {
  #   req(xlTable())
  #   unSelected <- map_lgl(xlTable()$id, function(x) input[[x]])
  #   removedXLtable <- xlTable() %>% filter(unSelected)
  #   xlTable(removedXLtable)
  # })
  # 
  
}
