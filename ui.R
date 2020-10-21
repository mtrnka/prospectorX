fluidPage(
    #    tags$head(tags$style(HTML(".shiny-text-output {background-color:#fff;}
    #                              hr {border-top: 0.5px solid #000000;}"))),
    titlePanel("ProspectorX-Touchstone Demo"),
    div(img(src="UCSF_logo_navy_RGB.png", height="75px"),
        styles="text-align: right"),
    br(),
    
    sidebarLayout(
        sidebarPanel(width=3,
                     h3("Parameter Selection"),
                     h4("Search Compare XL Output"),
                     column(4, style='padding:0px; margin:0px',
                            shinyFilesButton("clmsData", "Browse...",
                                             "Select the search compare crosslink output",
                                             multiple = F, viewtype = "detail")
                     ),
                     column(8, style='padding:0px; margin:0px',
                            verbatimTextOutput("clmsDataFileName"),
                            tags$style(type = 'text/css', '#clmsDataFileName {white-space:pre-wrap; padding:7px;}')
                     ),
                     h4("Module File"),
                     column(4, style='padding:0px; margin:0px',
                            shinyFilesButton("modules", "Browse...",
                                             "Select the module file",
                                             multiple = F, viewtype = "detail")
                     ),
                     column(8, style='padding:0px; margin:0px',
                            verbatimTextOutput("moduleFileName"),
                            tags$style(type = 'text/css', '#moduleFileName {white-space:pre-wrap; padding:7px;}')
                     ),
                     h4("PDB ID"),
                     column(4, style='padding:0px; margin:0px',
                            shinyFilesButton("pdbID", "Browse...",
                                             "Select the pdb/cif file",
                                             multiple = F, viewtype = "detail")
                     ),
                     column(8, style='padding:0px; margin:0px',
                            verbatimTextOutput("pdbFileName"),
                            tags$style(type = 'text/css', '#pdbFileName {white-space:pre-wrap; padding:7px;')
                     ),
                     h4("Chainmap File"),
                     column(4, style='padding:0px; margin:0px',
                            shinyFilesButton("chainmap", "Browse...",
                                             "Select the chainmap file",
                                             multiple = F, viewtype = "detail")
                     ),
                     column(8, style='padding:0px; margin:0px',
                            verbatimTextOutput("chainmapFileName"),
                            tags$style(type = 'text/css', '#chainmapFileName {white-space:pre-wrap; padding:7px;}')
                     ),
                     selectInput("summaryLevel", label = h4("Summarization Level"),
                                 choices = list("CSMs",
                                                "Unique Residue Pairs",
                                                "Protein Pairs",
                                                "Module Pairs")
                     ),
                     tags$hr(),
                     fluidRow(
                         column(6,
                                sliderInput("targetFDR", "Target FDR", min = 0, max = 25,
                                            step = 0.5, value = 1.0),
                                actionButton("findThreshold", "Classify")
                         ),
                         column(6,
                                br(),br(),
                                downloadButton("saveClassified", label = "Save results"),
                                br(),br(),
                                actionButton("viewXiNet", label = "View in xiNet"),
                                uiOutput("ui_open_tab"),
                         )
                     ),
                     tags$hr(),
                     fluidRow(
                         column(6,
                                sliderInput("svmThreshold", "Min. SVM Score", min = -5, max = 10, 
                                            step = 0.1, value = 0),
                                br(),
                                sliderInput("scoreDiffThreshold", "Min. Score.Diff", min = 0, max = 30,
                                            step = 0.5, value = 5),
                                br(),
                                sliderInput("peptideLengthFilter", "Peptide Length", min = 3, max = 35,
                                            step = 1, value = c(4, 30))
                         ),
                         column(6,
                                sliderInput("ms1MassError", "Prec Mass Error", min = -20, max = 20,
                                            step = 0.5, value = c(-25, 25)),
                                br(),
                                uiOutput("distanceSlider")
                         )
                     ),
                     tags$hr()
        ),
        mainPanel(width=9,
                  navbarPage(
                      title="Touchstone", id="navbar",
                      tabPanel("Classification",
                               h4("Console Output"),
                               verbatimTextOutput("consoleOut"),
                               h4("Classification Plots"),
                               verbatimTextOutput("numberClassedLinks"),
                               fluidRow(
                                   column(4,
                                          plotOutput("FDRplot"),
                                          verbatimTextOutput("FDR")
                                   ),
                                   column(4,
                                          plotOutput("massErrorPlot"),
                                          verbatimTextOutput("meanError")
                                   ),
                                   column(4,
                                          plotOutput("distancePlot"),
                                          verbatimTextOutput("VR")
                                   )
                               ),
                               br(),
                               fluidRow(
                                   column(4,
                                          plotOutput("thresholdPlot"),
                                          verbatimTextOutput("IIratio")
                                   ),
                                   column(4,
                                          plotOutput("proteinPlot",
                                                     click = "plot_click",
                                                     hover = "prot_hover"
                                                     #                   brush = brushOpts(
                                                     #                      id = "protPlot_brush"
                                                     #                 )
                                          ),
                                          verbatimTextOutput("protHover")
                                   ),
                                   column(4,
                                          plotOutput("modulePlot",
                                                     click = "plot_click",
                                                     hover = "mod_hover"
                                                     # brush = brushOpts(
                                                     #     id = "modPlot_brush"
                                                     # )
                                          ),
                                          verbatimTextOutput("modHover"),
                                   )
                               )
                      ),
                      tabPanel("Crosslink Table",
                               fluidRow(
                                   actionButton("viewXiNetTwo", label = "View in xiNet"),
                                   uiOutput("ui_open_tab_two"),
                                   br(),br(),
                                   DT::dataTableOutput("dataFile")
                               )
                      ),
                      tabPanel("Selected Crosslinks",
                               fluidRow(
                                   actionButton("viewXiNetSel", label = "View in xiNet"),
                                   uiOutput("ui_open_tab_sel"),
                                   br(),br(),
                                   DT::dataTableOutput("dataFileSelected")
                               )
                      ),
                      tabPanel("Decoy Hits",
                               fluidRow(
                                   DT::dataTableOutput("dataFileDecoy")
                               )
                      )
                  )
        )
    )
)
