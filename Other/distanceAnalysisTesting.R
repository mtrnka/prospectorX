library(tidyverse)
source("R/touchStone.R")

nameAccSwap <- function(dataTable) {
   proteinName.1 <- dataTable$Protein.1
   proteinName.2 <- dataTable$Protein.2
   dataTable$Protein.1 <- dataTable$Acc.1
   dataTable$Protein.2 <- dataTable$Acc.2
   dataTable$Acc.1 <- proteinName.1
   dataTable$Acc.2 <- proteinName.2
   return(dataTable)
}

pdb <- parsePDB("transloconDemo/4ug0.cif")
chains <- readChainMap("transloconDemo/chain_map.txt")
modFile <- readModuleFile("transloconDemo/modules.txt")

datTab <- new(Class="PPsearchCompareXL",
              dataFile="transloconDemo/etd200422a_res.txt",
              modFile="transloconDemo/modules.txt",
              pdbFile=pdb,
              chainMapFile=chains,
              preProcessFunction=nameAccSwap)

datTab <- buildClassifier(datTab)
dt <- getSearchTable(datTab)
dt <- bestResPairHackDT(dt)
dt <- dt %>% 
   filter(Decoy=="Target", dvals >= 0, !is.na(distance))

randomDistances <- getRandomCrosslinks(pdb, 5000)

distancePlot(pull(dt, distance), randomDistances, 35)

dtt <- dt %>% filter(dvals >= 2) %>% mutate(SVM.score = dvals)
test <- generateXVisOutput(dtt)
#write_csv(test, "~/Sites/xVis/docs/xvisTest.csv")

dttt <- dtt %>% filter(Modul.1 == "40S Ribosome", Modul.2 == "40S Ribosome")
modfile40S <- modFile[which(unlist(lapply(modFile, function(x) {x$Module=="40S Ribosome"})))]
.pairPlot(dttt, modfile40S)

################################################################################
## Making MS-Viewer Links -

library(urltools)

# vals[6] <- str_replace(vals[6], "(?<=\\/)[[A-Z]][[0-9]]+\\-.+?\\.mgf", dt[1,] %>% pull(Fraction))
# vals[18] <- dt[1,] %>% pull(RT)
# vals[21] <- dt[1,] %>% pull(z)
# vals[22] <- dt[1,] %>% pull(z)
# vals[23] <- dt[1,] %>% pull(Peptide.1)
# vals[25] <- dt[1,] %>% pull(Peptide.2)
# 
# test4 <- str_c(keys[1:28], vals[1:28], sep="=", collapse="&")
# test5 <- URLencode(test4)

temp <- "http://msviewer.ucsf.edu/prospector/cgi-bin/mssearch.cgi?search_name=msproduct&output_type=HTML&report_title=MS-Product&version=6.2.1%20Basic&data_source=Data%20From%20File&data_filename=%2Fvar%2Flib%2Fprospector%2Fweb%2Fresults%2Fmsviewer%2Fv%2Fd%2Fvdibnsypj7%2FZ190207_filt2%2FZ20190207-09etd.mgf&use_instrument_ion_types=1&msms_min_precursor_mass=0&instrument_name=ESI-ETD-high-res&display_graph=1&msms_parent_mass_tolerance=15&msms_parent_mass_tolerance_units=ppm&fragment_masses_tolerance=25&fragment_masses_tolerance_units=ppm&msms_pk_filter=Max%20MSMS%20Pks&msms_max_peaks=80&fraction=1&spot_number=86.533&run=1&spectrum_number=1&max_charge=5&msms_precursor_charge=5&sequence=HPGSFDVVHVK%28%2BDSS%29DANGNSFATR&s=1&sequence2=ASTSK%28%2BDSS%29SESSQK&s2=1&count_pos_z=Ignore%20Basic%20AA&link_search_type=DSS&"
temp2 <- unlist(str_split(temp, "&"))
temp3 <- str_split(temp2, "=")
keys <- map_chr(temp3, function(x) {x[1]})
vals <- map_chr(temp3, function(x) {x[2]})
vals <- url_decode(vals)

generateMSViewerLink <- function(templateKeys, templateVals, fraction, rt, z, peptide.1, peptide.2) {
   templateVals[6] <- str_replace(vals[6], "(?<=\\/)[[A-Z]][[0-9]]+\\-.+?\\.mgf", fraction)
   templateVals[18] <- rt
   templateVals[21] <- z
   templateVals[22] <- z
   templateVals[23] <- peptide.1
   templateVals[25] <- peptide.2
   templateVals <- url_encode(templateVals)
   zipped <- str_c(templateKeys[1:28], templateVals[1:28], sep="=", collapse="&")
   str_c('a href=\"', zipped, '\">link</a>')
}


