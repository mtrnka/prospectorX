library(shiny)
# library(googleComputeEngineR)
# library(tidyverse)
#source("touchstone.R")
# source("helperFunctions.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output, server, session) {
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
    
    prospX <- reactiveVal()
    prospX(NULL)
    consoleFile <- "gceRunOutput.txt"
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
                                    str_extract_all(runStream,
                                                    "\\<span.+$")
                                }
    )
    

    observeEvent(input$instances, {
        projectName <- gce_get_project()$name
        output$projectStatus <- renderText({projectName})
        instances <- fetchInstanceList()
        gceConnection <- prospX()
        if (!is.null(gceConnection)) {
            if(gceConnection$status == "RUNNING") {
                instances$status[1] <- "CONNECT-PROSPX"
            }
        }
        output$instanceTable <- renderTable({ instances })
    })
    
    observeEvent(input$connect, {
        instances <- fetchInstanceList()
        gceConnection <- gce_vm(instances$names[1])
        if (!is.null(gceConnection)) {
            if(gceConnection$status == "RUNNING") {
                instances$status[1] <- "CONNECT-PROSPX"
            }
        }
        output$instanceTable <- renderTable({ instances })
        prospX(gceConnection)
    })
    
    observeEvent(input$createNew, {
        # Needs to be connected to a prospector disk -
        gceConnection <- gce_vm_create(name = input$newInstanceName,
                                       predefined_type = input$gcpMachineType,
                                       image_family = "debian-9")
        prospX(gceConnection)
    })
    
    demoScriptLocation <- "bash ./runDemo.sh"
    observeEvent(input$submitSearch, {
            gceConnection <- prospX()
            gce_ssh(gceConnection, demoScriptLocation, 
                    capture_text = consoleFile,
                    wait = FALSE)
    })
    
    output$consoleOutput <- renderUI({ consoleRead() })
    

})

