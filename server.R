function(input, output, session) {
  options(shiny.maxRequestSize=1000*1024^2)
  
  consoleMessage <- reactiveVal("")
  output$consoleOut <- renderText({
    consoleMessage()
  })
  
  shinyFileChoose(input, "clmsData", roots=exDir, filetypes=c('', 'txt'))
  shinyFileChoose(input, "modules", roots=exDir, filetypes=c('', 'txt'))
  shinyFileChoose(input, "pdbID", roots=exDir, filetypes=c('', 'txt', 'pdb', 'cif'))
  shinyFileChoose(input, "chainmap", roots=exDir, filetypes=c('', 'txt'))

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
  
  output$pdbFileName <- renderPrint({
    if (is.integer(input$pdbID)) {
      cat("No file selected")
    } else {
      cat(parseFilePaths(exDir, input$pdbID)$name)
      consoleMessage("*** Reading Protein Structure File ***")
    }
  })
  
  output$chainmapFileName <- renderPrint({
    if (is.integer(input$chainmap)) {
      cat("No file selected")
    } else {
      cat(parseFilePaths(exDir, input$chainmap)$name)
      consoleMessage("*** Reading Chainmap File ***")
    }
  })
  
  scReader <- reactive({
    if (!is.integer(input$clmsData)) {
      inFile <- parseFilePaths(exDir, input$clmsData)
      consoleMessage("*** Loading Search Compare File ***")
      scTable <- readProspectorXLOutput(inFile$datapath)
      consoleMessage("*** Building SVM Classifier ***")
      scTable <- buildClassifier(scTable)
      return(scTable)
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
  
  pdbFile <- reactive({
    if (!is.integer(input$pdbID)) {
      inFile <- parseFilePaths(exDir, input$pdbID)
      return(inFile$datapath)
    }
  })
  
  chainmapFile <- reactive({
    if (!is.integer(input$chainmap)) {
      inFile <- parseFilePaths(exDir, input$chainmap)
      return(inFile$datapath)
    }
  })

  observeEvent(req(moduleFile()), {
    consoleMessage("*** Assigning Modules ***")
    scResults(assignModules(scResults(), moduleFile()))
  })

  observeEvent(req(chainmapFile(), pdbFile()), {
    consoleMessage("***calculating distances on PDB***")
    scResults(measureDistances(scResults(), parsePDB(pdbFile()), readChainMap(chainmapFile())))
  })

  csmTab <- reactive({
    req(scResults())
      datTab <- scResults()
      # The links should be generated in touchStone upon reading the SC file:
      msvFilePath <- parseFilePaths(exDir, input$clmsData)$datapath
      msvFiles <- system2("ls", c("-d", file.path(dirname(msvFilePath), "*/")), stdout=T)
      msvFiles <- str_replace(msvFiles, "\\/$", "")
      datTab <- datTab %>%
        mutate(link = pmap_chr(
          list(msvFiles, Fraction, RT, z, Peptide.1, Peptide.2, Spectrum), generateMSViewerLink))
      datTab <- generateCheckBoxes(datTab)
      minPPM = mfloor(min(datTab$ppm, na.rm=T))
      maxPPM = mfloor(max(datTab$ppm, na.rm=T))
      updateSliderInput(session, "ms1MassError", value = c(minPPM, maxPPM),
                        min = minPPM, max = maxPPM)
      minSVM = mfloor(min(datTab$dvals, na.rm=T), 1)
      maxSVM = mfloor(max(datTab$dvals, na.rm=T), 1)
      updateSliderInput(session, "svmThreshold", value = 0,
                        min = minSVM, max = maxSVM)
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
    tabLevelFiltered() %>% 
      filter(Decoy=="Target", dvals >= input$svmThreshold)
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
  
  output$dataFile <- DT::renderDataTable({
    req(csmTab())
    displayTable <- formatXLTable(xlTable())
    wideCols <- which(names(displayTable) %in% c("xlinkeResPair",
                                                 "DB.Peptide.1",
                                                 "DB.Peptide.2",
                                                 "Protein.1",
                                                 "Protein.2",
                                                 "Peptide.1",
                                                 "Peptide.2",
                                                 "Modul.1",
                                                 "Modul.2",
                                                 "Fraction"))
    DT::datatable(displayTable,
                  options = list(autoWidth=TRUE,
                                 deferRender=TRUE,
                                 processing=TRUE,
                                 columnDefs=list(list(width = '200px', targets = as.list(wideCols))),
                                 #                                 columnDefs=list(list(width = '200px', targets = c(3,7,8,17,21,26,27))),
                                 scrollX=TRUE,
                                 scrollY="80vh",
                                 scrollCollapse=TRUE,
                                 paging=TRUE,
                                 pageLength=100,
                                 search.caseInsensitive=TRUE,
                                 scroller=TRUE
                  ),
                  filter="top",
                  escape=FALSE
    )
  })
  
  output$dataFileSelected <- DT::renderDataTable({
    req(csmTab())
    displayTable <- formatXLTable(xlTableSelected())
    wideCols <- which(names(displayTable) %in% c("xlinkeResPair",
                                                 "DB.Peptide.1",
                                                 "DB.Peptide.2",
                                                 "Protein.1",
                                                 "Protein.2",
                                                 "Peptide.1",
                                                 "Peptide.2",
                                                 "Modul.1",
                                                 "Modul.2",
                                                 "Fraction"))
    DT::datatable(displayTable,
                  options = list(autoWidth=TRUE,
                                 deferRender=TRUE,
                                 processing=TRUE,
                                 columnDefs=list(list(width = '200px', targets = as.list(wideCols))),
                                 #                                 columnDefs=list(list(width = '200px', targets = c(3,7,8,17,21,26,27))),
                                 scrollX=TRUE,
                                 scrollY="80vh",
                                 scrollCollapse=TRUE,
                                 paging=TRUE,
                                 pageLength=100,
                                 search.caseInsensitive=TRUE,
                                 scroller=TRUE
                  ),
                  filter="top",
                  escape=FALSE
    )
  })

  fdr <- reactiveVal()
  observe({fdr(calculateFDR(tabLevelFiltered(), threshold = input$svmThreshold))})

  output$FDR <- renderText({
    req(tabLevel())
    str_c("FDR: ", as.character(round(100 * fdr(), 2)), "%")
  })

  output$FDRplot <- renderPlot({
    req(tabLevel())
    fdrPlots(tabLevelFiltered(), cutoff = input$svmThreshold)
  })

  numHits <- reactiveVal()
  observe({numHits(generateErrorTable(tabLevelFiltered()))})
  
  observeEvent(input$findThreshold, {
    req(tabLevel())
    threshold <- findThreshold(tabLevelFiltered(), targetER = input$targetFDR / 100)
    updateSliderInput(session, "svmThreshold", value = threshold[[1]])
  })

  output$IIratio <- renderText({
    str_c("percent interProtein: ", classRatio())
  })
  
  output$thresholdPlot <- renderPlot({
    req(numHits())
    numHitsPlot(numHits(), fdr())
  })
  
  randomDists <- reactive({
    req(pdbFile())
    consoleMessage("*** geting random Lys-Lys distances ***")
    getRandomCrosslinks(parsePDB(pdbFile()), 5000)
  })

  targetDists <- reactive({
    req(tabLevel())
    if (sum(is.na(tabLevel()$distance)) == nrow(tabLevel())) return(NULL)
    xlTable() %>% 
      pull(distance)
  })

  classedMassErrors <- reactive({
    subset(tabLevelFiltered(), dvals >= input$svmThreshold)$ppm
  })

  VR <- reactive({
    req(pdbFile())
    sum(targetDists() > input$distanceThreshold) / length(targetDists())
  })

  output$VR <- renderText({
    req(targetDists())
    str_c("Violation Rate: ", as.character(round(100 * VR(), 2)), "%"
    )
  })

  output$distancePlot <- renderPlot({
    req(pdbFile(), targetDists())
    distancePlot(targetDists(), randomDists(), threshold = input$distanceThreshold)
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
    req(pdbFile(), targetDists())
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
      scale_size(range = c(-1, 12)) +
      scale_color_viridis_c(option="D") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, hjust=1))
  })

    xlTableSelected <- reactive({
    req(xlTable(), input$plot_click)
    xlTable() %>% filter((Acc.1 == protFilter()[1] & Acc.2 == protFilter()[2]) |
                           (Acc.1 == protFilter()[2] & Acc.2 == protFilter()[1]),
    )
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
          xlTable() %>% filter((Modul.1 == protFilter()[1] & Modul.2 == protFilter()[2]) |
                                 (Modul.1 == protFilter()[2] & Modul.2 == protFilter()[1]))
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
      ggplot(aes(Modul.1, Modul.2, size=counts, col=counts)) + 
      geom_point(na.rm=T) +
      scale_size(range = c(-1, 12)) +
      scale_color_viridis_c(option="D") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, hjust=1))
  })
}
