library(shiny)
# library(googleComputeEngineR)
# library(tidyverse)
#source("touchstone.R")
# source("helperFunctions.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output, server) {
    options(shiny.maxRequestSize=1000*1024^2)
    # output$ms1fileName <- renderPrint({ str(input$ms1peakFile) })
    # output$ms2fileName <- renderPrint({ str(input$ms2peakFile) })
    # output$ms3fileName <- renderPrint({ str(input$ms3peakFile) })
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
    
    prospX <- NULL
    makeReactiveBinding("prospX")
    
    observeEvent(input$instances, {
        projectName <- gce_get_project()$name
        output$projectStatus <- renderText({projectName})
        instances <- fetchInstanceList()
        if (!is.null(prospX$status)) {
            if(prospX$status == "RUNNING") {
                instances$status[1] <- "CONNECT-PROSPX"
            }
        }
        output$instanceTable <- renderTable({ instances })
    })
    
    observeEvent(input$connect, {
        instances <- fetchInstanceList()
        prospX <- gce_vm(instances$names[1])
        if (!is.null(prospX$status)) {
            if(prospX$status == "RUNNING") {
                instances$status[1] <- "CONNECT-PROSPX"
            }
        }
        output$instanceTable <- renderTable({ instances })
    })
    
    observeEvent(input$createNew, {
        # Needs to be connected to a prospector disk -
        prospX <- gce_vm_create(name = input$newInstanceName,
                                predefined_type = input$gcpMachineType,
                                image_family = "debian-9")
    })
    
})
