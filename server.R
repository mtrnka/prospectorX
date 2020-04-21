library(shiny)
library(googleComputeEngineR)
library(tidyverse)

source("R/gceUtils.R")

consoleFile <- "gceRunOutput.txt"
demoScriptLocation <- "bash ./runDemo.sh"

if (file.exists(consoleFile)) {
    system2("rm", consoleFile)
    system2("touch", consoleFile)
}

function(input, output, session) {
    options(shiny.maxRequestSize=1000*1024^2)
    output$ms1fileName <- renderTable({ input$ms1peakFile[, 1:2] })
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
    
    output$consoleOutput <- renderUI({ h5(consoleRead()) })
    
    
}

