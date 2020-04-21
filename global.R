# load modules and do forth here.
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

