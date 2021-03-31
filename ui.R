fluidPage(
    tags$head(
        tags$style(HTML("hr {border-top: 0.75px solid #AAAAAA;}"))
    ),
    theme = "style.css",
    titlePanel("ProspectorX-Touchstone Demo"),
    div(img(src="UCSF_logo_navy_RGB.png", height="75px"),
        styles="text-align: right"),
    br(),
    
    sidebarLayout(
        sidebarPanel(width=3,
                     h3("Parameter Selection"),
                     fluidRow(
                         column(4, selectInput("experimentType", label = h4("Experiment"),
                                               choices = list("ms2", "ms3"))
                         ),
                         column(8, numericInput("scalingFactor", label = h4("DecoyDB Factor"),
                                            min = 1, max = 10, value = 10, step = 1)
                         )
                     ),
                     h4("Search Compare XL Output"),
                     fluidRow(
                         column(4, shinyFilesButton("clmsData", "Browse...",
                                                    "Select the search compare crosslink output",
                                                    multiple = F, viewtype = "detail")
                         ),
                         column(8, verbatimTextOutput("clmsDataFileName"),
                                tags$style(type = 'text/css', '#clmsDataFileName {white-space:pre-wrap; padding:7px;}')
                         )
                     ),
                     conditionalPanel(condition = "input.experimentType == 'ms3'",
                                      h4("MS2 Peaklists"),
                                      fluidRow(
                                          column(4, shinyFilesButton("ms2pkls", "Browse...",
                                                                     "Select the MS2 peaklist files",
                                                                     multiple = T, viewtype = "detail")
                                          ),
                                          column(8, verbatimTextOutput("ms2pklsName"),
                                                 tags$style(type = 'text/css', '#ms2pklsName {white-space:pre-wrap; padding:7px;}')
                                          )),
                                      h4("MS3 Peaklists"),
                                      fluidRow(
                                          column(4, shinyFilesButton("ms3pkls", "Browse...",
                                                                     "Select the MS3 peaklist files",
                                                                     multiple = T, viewtype = "detail")
                                          ),
                                          column(8, verbatimTextOutput("ms3pklsName"),
                                                 tags$style(type = 'text/css', '#ms3pklsName {white-space:pre-wrap; padding:7px;}')
                                          ))
                     ),
                     h4("Module File"),
                     fluidRow(
                         column(4, shinyFilesButton("modules", "Browse...",
                                                    "Select the module file",
                                                    multiple = F, viewtype = "detail")
                         ),
                         column(8, verbatimTextOutput("moduleFileName"),
                                tags$style(type = 'text/css', '#moduleFileName {white-space:pre-wrap; padding:7px;}')
                         )
                     ),
                     tags$hr(),
                     selectInput("summaryLevel", label = h4("Summarization Level"),
                                 choices = list("CSMs",
                                                "Unique Residue Pairs",
                                                "Protein Pairs")
                     ),
                     fluidRow(
                         column(6,
                                sliderInput("targetFDR", "Target FDR", min = 0, max = 25,
                                            step = 0.5, value = 1.0),
                                checkboxInput("separateFDRs", "Separate FDRs?",
                                              value = FALSE),
                                actionButton("findThreshold", "Classify")
                         ),
                         column(6,
                                br(),br(),
                                downloadButton("saveClassified", label = "Save results"),
                                br(),br(),
                                actionButton("viewXiNet", label = "View in xiNet"),
                                uiOutput("ui_open_tab"),
                                br(),br(),
                                actionButton("scrapeMSP", label = "Scrape MS Product")
                         )
                     ),
                     tags$hr(),
                     h3("Prefiltering"),
                     fluidRow(
                         column(6,
                                h4("min SVM Score"),
                                uiOutput("svmThreshSliders"),
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
