function(input, output, session) {
    options(shiny.maxRequestSize=1000*1024^2)

    exDir = c(wd= './transloconDemo')
    shinyFileChoose(input, "clmsData", roots=exDir, filetypes=c('', 'txt'))
    shinyFileChoose(input, "moduleFile", roots=exDir, filetypes=c('', 'txt'))
    shinyFileChoose(input, "pdbID", roots=exDir, filetypes=c('', 'txt', 'pdb', 'cif'))
    shinyFileChoose(input, "chainMapFile", roots=exDir, filetypes=c('', 'txt'))

    output$clmsDataFile <- renderPrint({
      if (is.integer(input$clmsData)) {
        cat("No file selected")
      } else {
        cat(parseFilePaths(exDir, input$clmsData)$name)
      }
    })

    output$modFile <- renderPrint({
      if (is.integer(input$moduleFile)) {
        cat("No file selected")
      } else {
        cat(parseFilePaths(exDir, input$moduleFile)$name)
        consoleMessage("*** Parsing Module File ***")
      }
    })

    output$pdbTabFile <- renderPrint({
      if (is.integer(input$pdbID)) {
        cat("No file selected")
      } else {
        cat(parseFilePaths(exDir, input$pdbID)$name)
        consoleMessage("*** Reading Protein Structure File ***")
      }
    })

    output$chainFile <- renderPrint({
      if (is.integer(input$chainMapFile)) {
        cat("No file selected")
      } else {
        cat(parseFilePaths(exDir, input$chainMapFile)$name)
        consoleMessage("*** Reading Chainmap File ***")
      }
    })

    pdbTab <- reactive({
      if (!is.integer(input$pdbID)) {
        inFile <- parseFilePaths(exDir, input$pdbID)
      return(parsePDB(inFile$datapath))
      }
    })

    chainTab <- reactive({
      if (!is.integer(input$chainMapFile)) {
        inFile <- parseFilePaths(exDir, input$chainMapFile)
        return(readChainMap(inFile$datapath))
      }
    })

    modTab <- reactive({
      if (!is.integer(input$moduleFile)) {
        inFile <- parseFilePaths(exDir, input$moduleFile)
        print("populate modules")
        return(readModuleFile(inFile$datapath))
      }
    })

    observe({
      if (!is.null(modTab())) {
        scRes(populateModules(scRes(), modTab()))
      }
    })

    scRead <- reactive({
      if (!is.integer(input$clmsData)) {
        inFile <- parseFilePaths(exDir, input$clmsData)
        consoleMessage("*** Loading Search Compare File ***")
        scRes <- readProspectorXLOutput(inFile$datapath)
        consoleMessage("*** Building SVM Classifier ***")
        scRes <- buildClassifier(scRes)
        return(scRes)
      }
    })

    scRes <- reactiveVal()
    observe({scRes(scRead())}) 
    
    
    consoleMessage <- reactiveVal("")
    output$consoleOut <- renderText({
      consoleMessage()
    })
    
    csmTab <- reactive({
      if (!is.null(scRes())) {
        datTab <- scRes()
        # The links should be generated in touchStone upon reading the SC file:
        datTab <- datTab %>%
          mutate(link = pmap_chr(
            list(Fraction, RT, z, Peptide.1, Peptide.2), generateMSViewerLink))
        datTab <- generateCheckBoxes(datTab)
        return(datTab)
      }
    })
    
    clTab <- reactive({
        consoleMessage("*** Calculating Residue-Pairs ***")
        bestResPairHackDT(csmTab())
    })

    tabLevel <- reactive({
      if (is.null(csmTab)) {return(NULL)}
      switch(input$summaryLevel,
             "CSMs" = csmTab(),
             "Unique Residue Pairs" = clTab()
      )
    })

    tabLevelFiltered <- reactive({
        #input$applyFilters, {
        if (is.null(csmTab)) {return(NULL)}
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
      if (is.null(csmTab())) return(NULL)
      tabLevelFiltered() %>% 
        filter(Decoy=="Target", dvals >= input$svmThreshold)
    })
    
    output$numberClassedLinks <- renderText({
      if (!is.null(tabLevel())) {
        paste0("Number of ", input$summaryLevel, ": ", nrow(xlTable()))
      }
    })
    
    output$dataFile <- DT::renderDataTable({
      if (is.null(csmTab()) ) return(NULL)
      DT::datatable(formatXLTable(xlTable()), filter="top", escape=FALSE)
    })
    
    output$FDR <- renderText({
      if (is.null(tabLevel())) return(NULL)
      paste0("FDR: ", 
             as.character(round(100 * calculateFDR(tabLevelFiltered(), threshold = input$svmThreshold), 2)),
             "%"
      )
    })

    output$FDRplot <- renderPlot({
      if (is.null(tabLevel())) return(NULL)
      fdrPlots(tabLevelFiltered(), cutoff = input$svmThreshold)
    })
    
    randomDists <- reactive({
      #        consoleMessage("*** geting random Lys-Lys distances ***")
      if (is.null(pdbTab())) return(NULL)
      getRandomCrosslinks(pdbTab(), 5000)
    })
    
    targetDists <- reactive({
      if (is.null(pdbTab())) return(NULL)
      tabLevelFiltered() %>% 
        filter(Decoy=="Target", dvals >= input$svmThreshold, !is.na(distance)) %>%
        pull(distance)
    })

    classedMassErrors <- reactive({
        subset(tabLevelFiltered(), dvals >= input$svmThreshold)$ppm
    })

    VR <- reactive({
      if (is.null(pdbTab())) return(NULL)
        sum(targetDists() > input$distanceThreshold) / length(targetDists())
    })

    output$VR <- renderText({
        if (is.null(targetDists())) return(NULL)
        paste0("Violation Rate: ", as.character(round(100 * VR(), 2)), "%"
        )
    })

    output$distancePlot <- renderPlot({
        if (is.null(targetDists())) return(NULL)
        distancePlot(targetDists(), randomDists(), threshold = input$distanceThreshold)
    })

    output$massErrorPlot <- renderPlot({
        if (is.null(tabLevel())) return(NULL)
        massErrorPlot(classedMassErrors(), 
                      lowPlotRange = floor(min(tabLevel()$ppm, na.rm=T)),
                      highPlotRange = ceiling(max(tabLevel()$ppm, na.rm=T)),
                      lowThresh = input$ms1MassError[1], 
                      highThresh = input$ms1MassError[2])
    })

    output$meanError <- renderText({
        if (is.null(tabLevel())) return(NULL)
        paste0("mean error: ", round(mean(classedMassErrors(), na.rm=T),2), " ppm")
    })    
}

