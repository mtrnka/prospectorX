# load modules and do forth here.
library(shiny)
library(shinyFiles)
library(tidyverse)
library(DT)
library(urltools)
library(e1071)

source("./touchStone.R")

mfloor <- function(x, base=5) {
   if (x < 0) {
      floored <- base * floor(x/base)
   } else if (x > 0) {
      floored <- base * ceiling(x/base)
   } else {
      floored <- 0
   } 
   return(floored)
}
