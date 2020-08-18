fluidPage(
    #    tags$head(tags$style(HTML(".shiny-text-output {background-color:#fff;}
    #                              hr {border-top: 0.5px solid #000000;}"))),
    titlePanel("ProspectorX Demo"),
    div(img(src="UCSF_logo_navy_RGB.png", height="75px"),
        styles="text-align: right"),
    br(),
   
    navbarPage(
        title=#div(img(src="UCSF_logo_navy_RGB.png", height="30px"), 
            "prospectorX-touchstone", 
        id="navbar",
        
        # Touchstone Analysis Tabset
        tabPanel("Analysis",
                 sidebarLayout(
                     sidebarPanel(width = 3,
                                  h3("File Selection"),
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
                                                             "Unique Residue Pairs")#,
                                              # "Protein Pairs",
                                              # "Domains",
                                              # "Modules")
                                  ),
                                  tags$hr(),
                                  fluidRow(
                                      column(6, 
                                             sliderInput("svmThreshold", "SVM Threshold", min = -5, max = 15, 
                                                         step = 0.1, value = 0),
                                             actionButton("applyFilters", h5("Apply Filters")),
                                             br(),br(),
                                             sliderInput("scoreDiffThreshold", "Score.Diff Threshold", min = 0, max = 30,
                                                         step = 0.5, value = 0),
                                             br(),
                                             sliderInput("peptideLengthFilter", "Min Peptide Length", min = 3, max = 7,
                                                         step = 1, value = 4)
                                      ),
                                      column(6,
                                             sliderInput("distanceThreshold", "Violation Distance", min = 0, max = 100,
                                                         step = 1, value = 30),
                                             actionButton("resetFilters", h5("Reset Filters")),
                                             br(),br(),
                                             sliderInput("ms1MassError", "Prec Mass Error", min = -20, max = 20,
                                                         step = 0.5, value = c(-5, 5))
                                      )
                                  )
                     ),
                     mainPanel(
                         h4("Console Output"),
                         verbatimTextOutput("consoleOut"),
                         h4("Classification Plots"),
                         verbatimTextOutput("numberClassedLinks"),
                         fluidRow(
                             column(4,
                                    verbatimTextOutput("FDR"),
                                    plotOutput("FDRplot")
                             ),
                             column(4,
                                    verbatimTextOutput("meanError"),
                                    plotOutput("massErrorPlot")
                             ),
                             column(4,
                                    verbatimTextOutput("VR"),
                                    plotOutput("distancePlot")
                             )
                         ),
                         fluidRow(
                             h3("Output Panel"),
                             DT::dataTableOutput("dataFile")
                         )
                     )
                 )
        )
    )
)