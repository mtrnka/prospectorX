function(input, output, session) {
    options(shiny.maxRequestSize=1000*1024^2)

    exDir = c(wd= './transloconDemo')
    shinyFileChoose(input, "moduleFile", roots=exDir, filetypes=c('', 'txt'))
    shinyFileChoose(input, "pdbID", roots=exDir, filetypes=c('', 'txt', 'pdb', 'cif'))
    shinyFileChoose(input, "chainMapFile", roots=exDir, filetypes=c('', 'txt'))
    shinyFileChoose(input, "clmsData", roots=exDir, filetypes=c('', 'txt'))

  
    output$etdMS2fileName <- renderTable({ input$etdMS2peakFile[, 1:2] })
    output$ms2fileName <- renderTable({ input$ms2peakFile[, 1:2] })
    output$ms3fileName <- renderTable({ input$ms3peakFile[, 1:2] })
    output$clReagentName <- renderPrint({ input$clReagent })
    output$ms1Tol <- renderPrint({ input$ms1Tol })
    output$ms2Tol <- renderPrint({ input$ms2Tol })
    output$ms3Tol <- renderPrint({ input$ms3Tol })
    output$secondaryDB <- renderUI({
        if (is.null(input$databaseType))
            return()
        switch(input$databaseType,
               "UniProtKB" = selectInput("species", h4("Species Code"),
                                         choices = c("HUMAN", "MOUSE", "YEAST",
                                                     "ECOLI", "etc"),
                                         width = "70%"),
               "SwissProt" = selectInput("species", h4("Species Code"),
                                         choices = c("HUMAN", "MOUSE", "YEAST",
                                                     "ECOLI", "etc"),
                                         width = "70%"),
               "fastaFile" = fileInput("fastaFile", h4("Fasta File"))
        )
    })
    output$tertiaryDB <- renderUI({
        if (is.null(input$databaseType))
            return()
        switch(input$databaseType,
               "UniProtKB" = textAreaInput("accNos", h4("List of Accession Numbers"),
                                           width = "70%"),
               "SwissProt" = textAreaInput("accNos", h4("List of Accession Numbers"),
                                           width = "70%"),
               "fastaFile" = ""
        )
    })
    
    output$ms3 <- reactive({
        str_detect(input$clStrategy, "MS3")
    })

    output$etdms2 <- reactive({
        str_detect(input$clStrategy, "ETD-MS2")
    })

    output$cidhcdms2 <- reactive({
        str_detect(input$clStrategy, "(HCD|CID)-MS2")
    })
    
    outputOptions(output, "ms3", suspendWhenHidden = FALSE)
    outputOptions(output, "etdms2", suspendWhenHidden = FALSE)
    outputOptions(output, "cidhcdms2", suspendWhenHidden = FALSE)
    
    # clStratMS2ETD <- reactive({
    #     str_detect(input$clStrategy, "ETD-MS2")
    # })
      
    # ToDo on GCE module.
    # Modularize the GCE code.
    # Re-think reactivity of the instance table.
    # New instances should load image from prospectorX disk image.

    pdbTab <- reactive({
      if (!is.integer(input$pdbID)) {
        inFile <- parseFilePaths(exDir, input$pdbID)
#        if (is.null(inFile)) return(NULL)
      return(parsePDB(inFile$datapath))
      }
    })

    chainTab <- reactive({
      if (!is.integer(input$chainMapFile)) {
        inFile <- parseFilePaths(exDir, input$chainMapFile)
        #        if (is.null(inFile)) return(NULL)
        return(readChainMap(inFile$datapath))
      }
    })

    csmTab <- reactive({
      if (!is.integer(input$clmsData) & !is.integer(input$moduleFile)) {
        inFile <- parseFilePaths(exDir, input$clmsData)
#        if (is.integer(inFile$file)) return(NULL)
        modFile <- parseFilePaths(exDir, input$moduleFile)
#        consoleMessage("*** Measuring PDB file crosslinks ***")
        datTab <- new(Class="PPsearchCompareXL",
                      dataFile=inFile$datapath,
                      modFile=modFile$datapath,
                      preProcessFunction=nameAccSwap,
                      chainMapFile=chainTab(),
                      pdbFile=pdbTab())
        # Renumber to account for N-terminal Flag and make residue numbers correspond
        # to Uniprot sequence
#        datTab <- renumberProtein(datTab,"Q9UM00ntf",-25)
        head(getSearchTable(datTab))
#        consoleMessage("*** Building SVM Classifier ***")
        datTab <- buildClassifier(datTab)
        datTab <- getSearchTable(datTab)
        datTab <- datTab %>%
            mutate(link = pmap_chr(
                list(Fraction, RT, z, Peptide.1, Peptide.2), generateMSViewerLink))
        datTab <- generateCheckBoxes(datTab)
        return(datTab)
      }
    })

#    consoleMessage <- reactiveVal("")
    
    # output$consoleOut <- renderText({
    #     consoleMessage()
    # })
    
    clTab <- reactive({
#        consoleMessage("*** Calculating Residue-Pairs ***")
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
        filter(Decoy=="Target",
               dvals >= input$svmThreshold)
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
        getRandomCrosslinks(pdbTab(), 5000)
    })

    targetDists <- reactive({
      if (is.null(csmTab())) return(NULL)
      tabLevelFiltered() %>% 
        filter(Decoy=="Target", dvals >= input$svmThreshold, !is.na(distance)) %>%
        pull(distance)
    })

    classedMassErrors <- reactive({
        subset(tabLevelFiltered(), dvals >= input$svmThreshold)$ppm
    })

    VR <- reactive({
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

