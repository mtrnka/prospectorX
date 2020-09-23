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
          list(msvFiles, Fraction, RT, z, Peptide.1, Peptide.2), generateMSViewerLink))
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
    bestResPairHackDT(csmTab())
  })
  
  tabLevel <- reactive({
    req(csmTab())
    switch(input$summaryLevel,
           "CSMs" = csmTab(),
           "Unique Residue Pairs" = clTab()
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
  
  output$numberClassedLinks <- renderText({
    req(tabLevel())
      str_c("Number of ", input$summaryLevel, ": ", nrow(xlTable()))
  })
  
  output$dataFile <- DT::renderDataTable({
    req(csmTab())
    DT::datatable(formatXLTable(xlTable()), filter="top", escape=FALSE)
  })

  fdr <- reactiveVal()
  observe({fdr(calculateFDR(tabLevelFiltered(), threshold = input$svmThreshold))})

  output$FDR <- renderText({
    req(tabLevel())
    str_c("FDR: ", as.character(round(100 * fdr(), 2)), "%")
  })

  output$FDRplot <- renderPlot({
    req(tabLevel())
    numHitsPlot(tabLevelFiltered(), cutoff = fdr())
  })

  numHits <- reactiveVal()
  observe({numHits(generateErrorTable(tabLevelFiltered()))})
  
  observeEvent(input$findThreshold, {
    req(tabLevel())
    threshold <- findThreshold(tabLevelFiltered(), targetER = input$targetFDR / 100)
    updateSliderInput(session, "svmThreshold", value = threshold[[1]])
  })

  output$thresholdPlot <- renderPlot({
    req()
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
    tabLevelFiltered() %>% 
      filter(Decoy=="Target", dvals >= input$svmThreshold, !is.na(distance)) %>%
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
}
