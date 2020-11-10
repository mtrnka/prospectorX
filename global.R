# load modules and do forth here.
library(shiny)
library(shinyFiles)
library(tidyverse)
library(DT)
library(urltools)
library(e1071)
source("./touchStone.R")

mmax <- function(x, base=5) {
   base * ceiling(x/base)
}

mmin <- function(x, base=5) {
   base * floor(x/base)
}

# Path to MS-Viewer Data
#exDir <- c(wd= './DemoFiles')
exDir <- c(wd= '/var/lib/prospector/seqdb/web/results/msviewer/')
#pathToXiFile <- "DemoFiles/xinetDemo/xiDemo.csv"
pathToXiFile <- "/var/www/crosslink-viewer/demo/data/"

# For generating MS-Viewer Links Correctly
queryTemplate <- "
http://rodin05.ucsf.edu/prospector/cgi-bin/mssearch.cgi?
search_name=msproduct&
output_type=HTML&
report_title=MS-Product&version=6.2.29&
data_source=Data%20From%20File&
data_filename=%2Fvar%2Flib%2Fprospector%2Fseqdb%2Fweb%2Fresults%2Fmsviewer%2F6%2Fk%2F6k6yrj8tuo%2FZ20200703_ethcd%2FZ20200703-05_FTMSms2ethcd.mgf&
use_instrument_ion_types=1&
msms_min_precursor_mass=0&
instrument_name=ESI-EThcD-high-res&display_graph=1&
msms_parent_mass_tolerance=10&
msms_parent_mass_tolerance_units=ppm&
fragment_masses_tolerance=20&
fragment_masses_tolerance_units=ppm&
msms_pk_filter=Max%20MSMS%20Pks&
msms_max_peaks=100&
fraction=1&
spot_number=59.210&
run=1&
spectrum_number=1&
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

queryTemplate.ms3 <- "
http://rodin05.ucsf.edu/prospector/cgi-bin/mssearch.cgi?
search_name=msproduct&
output_type=HTML&
report_title=MS-Product&
version=6.2.29&
data_source=Data%20From%20File&
data_filename=%2fvar%2flib%2fprospector%2fseqdb%2fweb%2fresults%2fmsviewer%2f3%2fk%2f3kobffhhss%2fdssoMS3%2fZ20200519-39_ITMSms3cid.mgf&
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
fraction=1&
spot_number=12.704&
run=1&
spectrum_number=1&
max_charge=5&
msms_precursor_charge=5&
sequence=ATGDETGAK%28Xlink:DSSO_s_fragment%29VER&
count_pos_z=Ignore%20Basic%20AA&
s=1
"

queryTemplate.ms3 <- unlist(str_split(str_replace_all(queryTemplate.ms3, "\\n", ""), "&"))
queryTemplate.ms3 <- str_split(queryTemplate.ms3, "=")
templateKeys.ms3 <- map_chr(queryTemplate.ms3, function(x) {x[1]})
templateVals.ms3 <- map_chr(queryTemplate.ms3, function(x) {x[2]})
templateVals.ms3 <- url_decode(templateVals.ms3)
