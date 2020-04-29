function(input, output, session) {
    options(shiny.maxRequestSize=1000*1024^2)
    output$etdMS2fileNamee <- renderTable({ input$etdMS2peakFile[, 1:2] })
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
    
    prospX <- reactiveVal()
    prospX(NULL)
    consoleStreamRegEx <- "(?<=span\\sid=.{1,20}\\>).+(?=\\<\\/span)"
    
    consoleRead <- reactivePoll(1000, session,
                                checkFunc = function() {
                                    if (file.exists(consoleFile))
                                        file.mtime(consoleFile)
                                    else
                                        ""
                                },
                                valueFunc = function() {
                                    runStream <- system2("tail", 
                                                         c("-n 1", consoleFile), 
                                                         stdout=T)
                                    str_extract_all(runStream, consoleStreamRegEx)
                                }
    )
    
    observe({
        if (req(input$navbar == "serverTab")) {
            projectName <- gce_get_project()$name
            output$projectStatus <- renderText({projectName})
            instances <- fetchInstanceList()
            selected <- as.integer(input$instanceNo)
            gceConnection <- prospX()
            positionInList <- str_which(gceConnection$name,
                                        instances$names[selected])
            if (!is.null(gceConnection)) {
                instances$status[positionInList] <- "CONNECT-PROSPX"
            }
            output$instanceTable <- renderTable({ instances })
        }        
    })
    
    observeEvent(input$instances, {
        projectName <- gce_get_project()$name
        output$projectStatus <- renderText({projectName})
        instances <- fetchInstanceList()
        selected <- as.integer(input$instanceNo)
        gceConnection <- prospX()
        positionInList <- str_which(gceConnection$name,
                                    instances$names[selected])
        if (!is.null(gceConnection)) {
            instances$status[positionInList] <- "CONNECT-PROSPX"
        }
        output$instanceTable <- renderTable({ instances })
        updateSelectInput(session, "instanceNo",
                          choices = instances$ID,
                          selected = head(instances$ID, 1)
        )
    })
    
    observeEvent(input$connect, {
        instances <- fetchInstanceList()
        selected <- as.integer(input$instanceNo)
        gceConnection <- gce_vm(instances$names[selected])
        positionInList <- str_which(instances$names,
                                    instances$names[selected])
        if (!is.null(gceConnection)) {
            instances$status[positionInList] <- "CONNECT-PROSPX"
        }
        output$instanceTable <- renderTable({ instances })
        prospX(gceConnection)
    })
    
    observeEvent(input$startVM, {
        instances <- fetchInstanceList()
        selected <- as.integer(input$instanceNo)
        gce_vm_start(instances$names[selected])
        gceConnection <- prospX()
        if (!is.null(gceConnection)) {
            instances$status[positionInList] <- "CONNECT-PROSPX"
            positionInList <- str_which(instances$names,
                                        gceConnection$name)
        }
    })
    
    observeEvent(input$stopVM, {
        instances <- fetchInstanceList()
        selected <- as.integer(input$instanceNo)
        gce_vm_stop(instances$names[selected])
        gceConnection <- prospX()
        if (!is.null(gceConnection)) {
            positionInList <- str_which(instances$names,
                                        gceConnection$name)
            instances$status[positionInList] <- "CONNECT-PROSPX"
        }
    })
    
    observeEvent(input$deleteVM, {
        instances <- fetchInstanceList()
        selected <- as.integer(input$instanceNo)
        gce_vm_delete(instances$names[selected])
        gceConnection <- prospX()
        if (!is.null(gceConnection)) {
            positionInList <- str_which(instances$names,
                                        gceConnection$name)
            instances$status[positionInList] <- "CONNECT-PROSPX"
        }
    })
    
    observeEvent(input$createNew, {
        # Needs to be connected to a prospector disk -
        gceConnection <- gce_vm_create(name = input$newInstanceName,
                                       predefined_type = input$gcpMachineType,
                                       image_family = "debian-9")
        prospX(gceConnection)
    })
    
    observeEvent(input$submitSearch, {
        gceConnection <- prospX()
        gce_ssh(gceConnection, demoScriptLocation, 
                capture_text = consoleFile,
                wait = FALSE)
    })
    
    output$consoleOutput <- renderUI({ h4(consoleRead()) })
    output$gcePricing <- renderTable({ gcePrices.OR })

    pdbTab <- reactive({
        inFile <- input$pdbID
        if (is.null(inFile)) return(NULL)
        return(parsePDB(inFile$datapath))
    })

    chainTab <- reactive({
        inFile <- input$chainMapFile
        if (is.null(inFile)) return(NULL)
        return(readChainMap(inFile$datapath))
    })

    csmTab <- reactive({
        inFile <- input$clmsData
        if (is.null(inFile)) return(NULL)
        modFile <- input$moduleFile
        consoleMessage("*** Measuring PDB file crosslinks ***")
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
        consoleMessage("*** Building SVM Classifier ***")
        datTab <- buildClassifier(datTab)
        datTab <- getSearchTable(datTab)
        datTab <- datTab %>% 
            mutate(link = pmap_chr(
                list(Fraction, RT, z, Peptide.1, Peptide.2), generateMSViewerLink))
        datTab <- generateCheckBoxes(datTab)
        return(datTab)
    })

    consoleMessage <- reactiveVal("")
    
    output$consoleOut <- renderText({
        consoleMessage()
    })
    
    clTab <- reactive({
        consoleMessage("*** Calculating Residue-Pairs ***")
        bestResPairHackDT(csmTab())
    })

    tabLevel <- reactive({
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
        if (is.null(tabLevel())) return(NULL)
        paste0("Number of ", input$summaryLevel, ": ", nrow(xlTable()))
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
        consoleMessage("*** geting random Lys-Lys distances ***")
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


