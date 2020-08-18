function(input, output, session) {
  options(shiny.maxRequestSize=1000*1024^2)
  
  consoleMessage <- reactiveVal("")
  output$consoleOut <- renderText({
    consoleMessage()
  })
  
  exDir = c(wd= './DemoFiles')
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
      scResults <- readProspectorXLOutput(inFile$datapath)
      consoleMessage("*** Building SVM Classifier ***")
      scResults <- buildClassifier(scResults)
      return(scResults)
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
    scResults(populateModules(scResults(), readModuleFile(moduleFile())))
  })

  observeEvent(req(chainmapFile(), pdbFile()), {
    consoleMessage("***calculating distances on PDB***")
    scResults(measureDistances(scResults(), parsePDB(pdbFile()), readChainMap(chainmapFile())))
  })

  csmTab <- reactive({
    req(scResults())
      datTab <- scResults()
      # The links should be generated in touchStone upon reading the SC file:
      datTab <- datTab %>%
        mutate(link = pmap_chr(
          list(Fraction, RT, z, Peptide.1, Peptide.2), generateMSViewerLink))
      datTab <- generateCheckBoxes(datTab)
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
                          Len.Pep.1 >= input$peptideLengthFilter,
                          Len.Pep.2 >= input$peptideLengthFilter,
                          ppm >= input$ms1MassError[1],
                          ppm <= input$ms1MassError[2])
  })
  
  observeEvent(input$resetFilters, {
    updateSliderInput(session, "ms1MassError",
                      min = floor(min(tabLevel()$ppm, na.rm=T)),
                      max = ceiling(max(tabLevel()$ppm, na.rm=T))
    )
    updateSliderInput(session, "svmThreshold", 
                      min = floor(min(tabLevel()$dvals)),
                      max = ceiling(max(tabLevel()$dvals))
    )
  })
  
  xlTable <- reactive({
    req(csmTab())
    tabLevelFiltered() %>% 
      filter(Decoy=="Target", dvals >= input$svmThreshold)
  })
  
  output$numberClassedLinks <- renderText({
    req(tabLevel())
      paste0("Number of ", input$summaryLevel, ": ", nrow(xlTable()))
  })
  
  output$dataFile <- DT::renderDataTable({
    req(csmTab())
    DT::datatable(formatXLTable(xlTable()), filter="top", escape=FALSE)
  })
  
  output$FDR <- renderText({
    req(tabLevel())
    paste0("FDR: ", 
           as.character(round(100 * calculateFDR(tabLevelFiltered(), threshold = input$svmThreshold), 2)),
           "%"
    )
  })

  output$FDRplot <- renderPlot({
    req(tabLevel())
    fdrPlots(tabLevelFiltered(), cutoff = input$svmThreshold)
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
    paste0("Violation Rate: ", as.character(round(100 * VR(), 2)), "%"
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
    paste0("mean error: ", round(mean(classedMassErrors(), na.rm=T),2), " ppm")
  })    
}

