# script to run touchstone from command line.
source("./global.R")

readTstoneParams <- function(paramsFile) {
   tstoneParams <- read_json(paramsFile)$parameters
   return(tstoneParams)
}

processXLresult <- function(paramsSCresult) {
   scTable <- switch(paramsSCresult$ms_level,
                     "MS2" = processMS2result(paramsSCresult),
                     "MS3" = processMS3result(paramsSCresult),
                     stop("Invalid MS Experiment Type")
   )
   # Load Modules.
   if (!is.null(paramsSCresult$module_file)) {
      modFile <- readModuleFile(paramsSCresult$module_file)
      scTable <- tryCatch(assignModules(scTable, modFile),
                          error = function(e) {
                             message(str_c("Problem reading Module File-", e, sep="\n"))
                             return(scTable)
                          }
      )
   }
   # Do Distance Calculation.
   if (!is.null(modFile)) {
      pdbFiles <- gatherPDBs(modFile, paramsSCresult$pdb_directory)
      if (length(pdbFiles) > 0) {
         scTable <- tryCatch(measureCrosslinkDistance(scTable, pdbFiles),
                             error = function(e) {
                                message(str("Problem with distance measurements", e, sep="\n"))
                                return(scTable)
                             }
         )
      }
   }
   # Write the touchstone annotated file.
      
   # Select Reporting Level.
   # Select Filters.
   # Select Split or Unified.
   # Generate crosslink table for display.
   # Calculate Error Table.
   # Update Params File.
   # Write Output.
}

# csmTab <- scResults with checkboxes and links to MS-Product.  Sliders are updated.

#clTab <- residuePairs(csmTab)
#plTab <- protPairs(csmTab)
#modTab <- modPairs(csmtab)

# tabLevel <- reactive({
#    req(csmTab())
#    switch(input$summaryLevel,
#           "CSMs" = csmTab(),
#           "Unique Residue Pairs" = clTab(),
#           "Protein Pairs" = plTab(),
#           "Module Pairs" = modTab()
#    )
# })
# 
# tabLevelFiltered <- reactive({
#    req(csmTab())
#    tabLevel() %>% filter(Score.Diff >= input$scoreDiffThreshold,
#                          between(Len.Pep.1, input$peptideLengthFilter[1], input$peptideLengthFilter[2]),
#                          between(Len.Pep.2, input$peptideLengthFilter[1], input$peptideLengthFilter[2]),
#                          between(ppm, input$ms1MassError[1], input$ms1MassError[2]))
# })
# 
# xlTable <- reactive({
#    req(csmTab())
#    if (input$separateFDRs) {
#       req(input$svmThresholdIntra, input$svmThresholdInter)
#       tabLevelFiltered() %>% 
#          filter(Decoy=="Target",
#                 (str_detect(xlinkClass, "intraProtein") & dvals >= input$svmThresholdIntra) |
#                    (str_detect(xlinkClass, "interProtein") & dvals >= input$svmThresholdInter))
#    } else {
#       req(input$svmThreshold)
#       tabLevelFiltered() %>% 
#          filter(Decoy=="Target", dvals >= input$svmThreshold)
#    }
# })
# 
# xlTableDecoy <- reactive({
#    req(csmTab())
#    if (input$separateFDRs) {
#       req(input$svmThresholdIntra, input$svmThresholdInter)
#       tabLevelFiltered() %>% 
#          filter(Decoy!="Target",
#                 (str_detect(xlinkClass, "intraProtein") & dvals >= input$svmThresholdIntra) |
#                    (str_detect(xlinkClass, "interProtein") & dvals >= input$svmThresholdInter))
#    } else {
#       req(input$svmThreshold)
#       tabLevelFiltered() %>% 
#          filter(Decoy!="Target", dvals >= input$svmThreshold)
#    }
# })

processMS2result <- function(paramsSCresult) {
   scTable <- readProspectorXLOutput(paramsSCresult$scFile_name)
   scTable <- buildClassifier(scTable)
   return(scTable)
}

processMS3result <- function(paramsSCresult) {
   paramsName <- deparse(substitute(paramsSCresult))
   scan_info <- paramsSCresult$experiment_scan_info
   if(!is.null(scan_info$master_scan_file)) {
      masterScanFile <- read_tsv(file.path(scan_info$peaklist_dir,
                                           scan_info$master_scan_file))
   } else {
      peaklists <- unlist(scan_info$peaklists)
      ms2peaklists <- file.path(scan_info$peaklist_dir, names(peaklists))
      ms3peaklists <- file.path(scan_info$peaklist_dir, peaklists)
      masterScanFile <- createMasterScanFile(ms3peaklists, ms2peaklists)
      tryCatch(write_tsv(masterScanFile, file.path(scan_info$peaklist_dir, 
                                                   str_c(paramsName, "scanFile.txt"))),
               error = function(e) {
                  message(str_c("problem writing master scan file", e, "\n"))
               }
      )
   }
   searchCompareMS3 <- readMS3results(paramsSCresult$scFile_name)
   searchCompareMS3 <- addMasterScanInfo(searchCompareMS3, masterScanFile)
   scTable <- processMS3xlinkResults(searchCompareMS3)
   return(scTable)
}

tstone3 <- readTstoneParams("tstoneParamsMS3test.json")
testMS3 <- tstone3$searchCompare_results$scResults_1
testMS3$experiment_scan_info$master_scan_file <- "paramsSCresultscanFile.txt"
demo3 <- processXLresult(testMS3)

tstone2 <- readTstoneParams("tstoneParamsMS2test2.json")
testMS2 <- tstone2$searchCompare_results$scResults_1
demo2 <- processXLresult(testMS2)

newMod <- readModuleFile(testMS2$module_file)
#pdbFiles <- gatherPDBs(newMod, pdbFileDir = testMS2$pdb_directory)

#demo3 <- measureCrosslinkDistance(demo2, pdbFiles)
demo2 %>% filter(SVM.score > 0) %>%
   ggplot(aes(distance, fill=PDB)) + 
   geom_histogram(binwidth = 2.5, col="black", position="stack") +
   geom_vline(xintercept = 35) + 
   scale_fill_brewer(palette=4, type="seq") + 
   theme_bw()

