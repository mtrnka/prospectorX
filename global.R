# load modules and do forth here.
library(shiny)
library(shinyFiles)
library(tidyverse)
library(DT)
library(urltools)

source("transloconDemo/touchStoneWorking.R")

nameAccSwap <- function(dataTable) {
   proteinName.1 <- dataTable$Protein.1
   proteinName.2 <- dataTable$Protein.2
   dataTable$Protein.1 <- dataTable$Acc.1
   dataTable$Protein.2 <- dataTable$Acc.2
   dataTable$Acc.1 <- proteinName.1
   dataTable$Acc.2 <- proteinName.2
   return(dataTable)
}




