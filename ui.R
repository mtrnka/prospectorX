fluidPage(
    tags$head(tags$style(HTML(".shiny-text-output {background-color:#fff;}"))),
    
    titlePanel("ProspectorX Demo"),
    div(img(src="UCSF_logo_navy_RGB.png", height="75px"),
        styles="text-align: right"),
    br(),
    
    navbarPage("prospectorX-touchstone", 

               # Server Connection Tabset
               tabPanel("Server Connection", fluidRow(
                   column(5,
                          wellPanel(
                              h3("Server Information"),
                              helpText("This area contains parameter that enable \
                                            connections either to GCP instance or a users \
                                            local server (eg Wynton)"),
                              img(src="gcp.png", height= "100px"),
                              fluidRow(
                                  column(2, h4("Actions"),
                                         actionButton("instances", "View",
                                                      width = "100px"),
                                         actionButton("connect", "Connect",
                                                      width = "100px"),
                                         actionButton("startVM", "Start VM",
                                                      width = "100px"),
                                         actionButton("stopVM", "Stop VM",
                                                      width = "100px"),
                                         actionButton("deleteVM", "Delete VM",
                                                      width = "100px"),
                                         selectInput("instanceNo", h4("VM ID"),
                                                     choices = c(""),
                                                     width = "100px")
                                  ),
                                  column(10, h4("GCP Project Name"),
                                         verbatimTextOutput("projectStatus"),
                                         h4("Instances"),
                                         tableOutput("instanceTable")
                                  )
                              )
                          ),
                   ),
                   
                   column(7,
                          wellPanel(
                              h3("Create New VM"),
                              fluidRow(
                                  column(3,
                                         textInput("newInstanceName",
                                                   h4("Instance Name")),
                                         selectInput("gcpMachineType",
                                                     label = h4("Machine Type"),
                                                     choices = c("f1-micro",
                                                                 "n1-standard-1",
                                                                 "n1-standard-2",
                                                                 "n1-standard-4",
                                                                 "n1-standard-8",
                                                                 "n1-standard-16",
                                                                 "n1-standard-32",
                                                                 "n1-standard-64"),
                                                     selected = "f1-micro",
                                                     width = "100%"),
                                         actionButton("createNew",
                                                      "Create")
                                  ),
                                  
                                  column(9,
                                         h4("Compute Engine Pricing - Oregon (us-west1)"),
                                         tableOutput("gcePricing")
                                  )
                              )
                          )
                   )
               )),
               
               # Server Connection Tabset
               tabPanel("Job Submission", fluidRow(
                   # Biological input area
                   column(3,
                          wellPanel(#style = "background: #E3E7AF",
                              h3("Biological System Input"),
                              helpText("This are includes database info for the search \
                                      but also optional module and pdb files that will be used \
                                      in data display and downstream analysis"),
                              selectInput("databaseType", h4("Sequence Database Type"),
                                          choices = c("UniProtKB",
                                                      "SwissProt",
                                                      "fastaFile"),
                                          selected = "fastaFile",
                                          width = "70%"),
                              uiOutput("secondaryDB"),
                              uiOutput("tertiaryDB"),
                              fileInput("moduleFile",
                                        label = h4("Module File")),
                              fileInput("pdbID",
                                        label = h4("pdb ID")),
                              fileInput("chainMapFile", label = h4("Chainmap File"))
                          )
                   ),
                   
                   column(5,
                          wellPanel(#style = "background: #BFC0C0",
                              h3("Peaklist Input"),
                              helpText("This area contains specifications for the \
                                            peaklist files... could be expanded to include raw \
                                            files and peaklist generation"),
                              uiOutput("ms1info"),
                              uiOutput("ms2info"),
                              uiOutput("ms3info"),
                              fileInput("ms1peakFile", h4("ms1peakFiles"),
                                        multiple = T),
                              tableOutput("ms1fileName"),
                              fileInput("ms2peakFile", h4("ms2peakFiles"),
                                        multiple = T),
                              tableOutput("ms2fileName"),
                              fileInput("ms3peakFile", h4("ms3peakFiles"),
                                        multiple = T),
                              tableOutput("ms3fileName")
                          ),
                          column(12,
                                 fluidRow(
                                     wellPanel(
                                         h3("Submit"),
                                         actionButton("submitSearch", "Submit", icon("paper-plane"),
                                                      style="color: #fff;
                                                   background-color: #337ab7;
                                                   border-color: #2e6da4,"),
                                         br(),
                                         htmlOutput("consoleOutput")
                                     )
                                 )
                          )
                   ),
                   
                   column(4,
                          wellPanel(#style = "background: #CEE7E6",
                              h3("MS Experiment Input"),
                              selectInput("clStrategy", label = h4("MS Strategy"),
                                          choices = list("CID-MS2-CID-MS3",
                                                         "CID-MS3",
                                                         "CID-MS2-HCD-MS2",
                                                         "CID-MS2-CID-MS3-ETD-MS2",
                                                         "stHCD-MS2",
                                                         "ETD-MS2",
                                                         "ETD-MS2-HCD-MS2"),
                                          selected = 1,
                                          width = "100%"),
                              
                              fluidRow(
                                  column(6,
                                         selectInput("clReagent", label = h4("Crosslinker"),
                                                     choices = list("DSSO" = 1,
                                                                    "DSS/BS3" = 2,
                                                                    "EDC" = 3,
                                                                    "etc..." = 4),
                                                     selected = 1,
                                                     width = "100%"),
                                         
                                         selectInput("Protease", label = h4("Enzyme"),
                                                     choices = list("Trypsin",
                                                                    "Trypsin/GluC",
                                                                    "GluC",
                                                                    "LysC",
                                                                    "etc"),
                                                     selected = 1,
                                                     width = "100%"),
                                  ),
                                  
                                  column(6,
                                         selectInput("clQuench", label = h4("Quencher"),
                                                     choices = list("Ammonia" = 1,
                                                                    "Tris" = 2),
                                                     selected = 1,
                                                     width = "100%"),
                                         
                                         sliderInput("MissedCl", label = h4("Mis. Cleavages"),
                                                     min = 0, max = 3, value = 1, width = "50%"),
                                  )
                              ),
                              
                              fluidRow(
                                  column(6,
                                         sliderInput("ms1Tol", label = h4("MS1 tolerance (ppm)"),
                                                     min = 0, max = 50, value = 10, width = "70%"),
                                         #verbatimTextOutput("ms1Tol"),
                                         
                                         sliderInput("ms2Tol", label = h4("MS2 tolerance (ppm)"),
                                                     min = 0, max = 50, value = 20, width = "70%"),
                                         #verbatimTextOutput("ms2Tol"),
                                  ),
                                  
                                  column(6,
                                         sliderInput("ms3Tol", label = h4("MS3 tolerance (Da)"),
                                                     min = 0, max = 2.0, value = 0.7, step = 0.1, width = "70%")
                                         #verbatimTextOutput("ms3Tol")
                                  )
                              )
                          )
                   )
               )),
               
               # Server Connection Tabset
               tabPanel("Analysis", fluidRow())
    )
)