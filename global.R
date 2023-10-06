# load modules and so forth here.
library(shiny)
library(shinyFiles)
library(tidyverse)
library(DT)
library(urltools)
library(e1071)
library(xml2)
library(jsonlite)
library(bio3d)
library(future)
library(furrr)
library(minpack.lm)
future::plan(multicore)
options(future.globals.maxSize = 8000 * 1024^2)
source("./touchStone.R")

mmax <- function(x, base=5) {
   base * ceiling(x/base)
}

mmin <- function(x, base=5) {
   base * floor(x/base)
}

readParamsFile <- function(paramsFile) {
   tryCatch(read_xml(paramsFile),
            error = function(e) {
               message(e)
               return(NA)
            }
   )
}
#SVM score parameters:
params.default <- c("Score.Diff","z","Score","numCSM","massError", "Rk.2","Rk.1")
params.best <- c("Score.Diff", "percMatched", "massError", "z", "numURP", "numCSM", "xlinkClass")
params.best.nop <- c("Score.Diff", "percMatched", "massError", "z", "numCSM")
   
proton <- 1.007276
H2O <- 18.01002
SOH2 <- 49.98209
sulfur <- 31.97152
   
# Path to MS-Viewer Data
#exDir <- c(wd= './DemoFiles')
exDir <- c(wd= '/mnt/pipeline/projects')
#pathToXiFile <- "DemoFiles/xinetDemo/xiDemo.csv"
pathToXiFile <- "/var/www/html/crosslink-viewer/demo/data"
linkToXiView <- 'http://lanhuang2.physiology.uci.edu/crosslink-viewer/demo/Demo2.html'

# For generating MS-Viewer Links Correctly
queryTemplate <- "
http://lanhuang2.physiology.uci.edu/prospector/cgi-bin/mssearch.cgi?
search_name=msproduct&
output_type=HTML&
report_title=MS-Product&version=6.2.29&
data_source=Data%20From%20File&
data_filename=%2Fmnt%2Fpipeline%2Fprojects%2FDSSOstar_rRibo%2Fsc%2FrRibo_DSSO_star_HCD%2FtstoneMS2.1%2FDSSOstar_rRibo%2FZ20200519-63_FTMSms2hcd.mgf&
use_instrument_ion_types=1&
msms_min_precursor_mass=0&
instrument_name=ESI-EThcD-high-res&display_graph=1&
msms_parent_mass_tolerance=10&
msms_parent_mass_tolerance_units=ppm&
fragment_masses_tolerance=30&
fragment_masses_tolerance_units=ppm&
msms_pk_filter=Max%20MSMS%20Pks&
msms_max_peaks=100&
scan_number=1000&
max_charge=4&
msms_precursor_charge=4&
sequence=SQK%28%2BDSG%29AIQDEIR&
s=1&
sequence2=Q%28Gln-%3Epyro-Glu%29QLPLPYEQLK%28%2BDSG%29HFYR&
s2=1&
count_pos_z=Ignore%20Basic%20AA&
link_search_type=No%20Link&
"

queryTemplate <- unlist(str_split(str_replace_all(queryTemplate, "\\n", ""), "&"))
queryTemplate <- str_split(queryTemplate, "=")
templateKeys <- map_chr(queryTemplate, function(x) {x[1]})
templateVals <- map_chr(queryTemplate, function(x) {x[2]})
templateVals <- url_decode(templateVals)
names(templateVals) <- templateKeys
templateVals <- templateVals[!is.na(templateVals)]

queryTemplate.ms3 <- "
http://lanhuang2.physiology.uci.edu/prospector/cgi-bin/mssearch.cgi?
search_name=msproduct&
output_type=HTML&
report_title=MS-Product&
version=6.2.29&
data_source=Data%20From%20File&
data_filename=
%2Fmnt%2Fpipeline%2Fprojects%2FDSSOms3_rRibo%2Fsc%2FrRibo_DSSO_ms3%2FtstoneMS3.1%2FDSSOms3_rRibo%2FZ20200519-49_ITMSms3cid.mgf&
use_instrument_ion_types=1&
msms_min_precursor_mass=0&
instrument_name=ESI-ION-TRAP-low-res&
display_graph=1&
msms_parent_mass_tolerance=10&
msms_parent_mass_tolerance_units=ppm&
fragment_masses_tolerance=0.7&
fragment_masses_tolerance_units=Da&
msms_pk_filter=Max%20MSMS%20Pks&
msms_max_peaks=40&
scan_number=37130&
max_charge=5&
msms_precursor_charge=5&
sequence=DHASIQMNVAEVDK%28XL:A-Alkene%29VTGR
count_pos_z=Ignore%20Basic%20AA&
s=1
"

queryTemplate.ms3 <- unlist(str_split(str_replace_all(queryTemplate.ms3, "\\n", ""), "&"))
queryTemplate.ms3 <- str_split(queryTemplate.ms3, "=")
templateKeys.ms3 <- map_chr(queryTemplate.ms3, function(x) {x[1]})
templateVals.ms3 <- map_chr(queryTemplate.ms3, function(x) {x[2]})
templateVals.ms3 <- url_decode(templateVals.ms3)
names(templateVals.ms3) <- templateKeys.ms3
templateVals.ms3 <- templateVals.ms3[!is.na(templateVals.ms3)]
