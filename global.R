# load modules and do forth here.
library(shiny)
library(googleComputeEngineR)
library(tidyverse)
library(DT)

source("R/gceUtils.R")
source("transloconDemo/touchStoneWorking.R")


consoleFile <- "gceRunOutput.txt"
demoScriptLocation <- "bash ./runDemo.sh"

if (file.exists(consoleFile)) {
   system2("rm", consoleFile)
   system2("touch", consoleFile)
}

nameAccSwap <- function(dataTable) {
   proteinName.1 <- dataTable$Protein.1
   proteinName.2 <- dataTable$Protein.2
   dataTable$Protein.1 <- dataTable$Acc.1
   dataTable$Protein.2 <- dataTable$Acc.2
   dataTable$Acc.1 <- proteinName.1
   dataTable$Acc.2 <- proteinName.2
   return(dataTable)
}




