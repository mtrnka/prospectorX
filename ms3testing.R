source("global.R")
scResults2 <- "DemoFiles/dssoMS3/MS3_DSSO_2.txt"
ms3files2 <- c("DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-31_ITMSms3cid.txt",
              "DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-39_ITMSms3cid.txt",
              "DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-49_ITMSms3cid.txt",
              "DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-59_ITMSms3cid.txt")

ms2files2 <- c("DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-31_FTMSms2cid.txt",
              "DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-39_FTMSms2cid.txt",
              "DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-49_FTMSms2cid.txt",
              "DemoFiles/dssoMS3/rRiboDSSOms3/Z20200519-59_FTMSms2cid.txt")
demo <- processMS3xlinkResultsMultiFile(scResults2, ms3files2, ms2files2)

searchCompareMS3 <- readMS3results(scResults2)
searchCompareMS3 <- searchCompareMS3 %>% group_by(Fraction.ms3) %>% nest()
baseFiles <- str_replace(searchCompareMS3$Fraction.ms3, "_[[A-Z]]+MSms[[0-9]][[a-z]]+", "")
p1 <- searchCompareMS3 %>% filter(str_detect(Fraction.ms3, baseFiles[1])) %>% unnest(cols=c(data))
p2 <- ms3files2[str_which(ms3files2, baseFiles[1])]
p3 <- ms2files2[str_which(ms2files2, baseFiles[1])]

parseProteinMods <- function(ms3tab) {
   mods <- str_split(ms3tab$Protein.Mods.ms3, ";")
   ms3tab$thiolMods <- map_int(mods, function(x) sum(str_count(x, "XL:A-Thiol\\(Unsaturated\\)")))
   ms3tab$alkeneMods <- map_int(mods, function(x) sum(str_count(x, "XL:A-Alkene")))
   ms3tab$sulfenicMods <- map_int(mods, function(x) sum(str_count(x, "XL:A-Sulfenic")))
   #ms3tab$thiolPos <- map_chr(mod, )
   return(ms3tab)
}

p1a <- parseProteinMods(p1)
test <- p1a %>% filter(ms3ScanNo %in% c(20539, 10075))
test <- test$Protein.Mods.ms3
test <- unlist(str_split(test, ";"))
test2 <- test %>% str_split("@") %>% map_chr(function(x) x[2])

gregexpr("\\((?>[^()]|(?R))*\\)", test2, perl=T)
gregexpr("\\((?>[^()]|(?R))*\\)", 
         substr(test2, 2, 48), perl=T)
str_sub(str_sub(test2, 2, 48), c(1, 17, 33), c(1+14, 17+14, 33+14))
#processMS3xlinkResultsSingleFile(p1, p2, p3)

testFunction <- function(parantheticalThing) {
   print(parantheticalThing)
   text <- str_sub(parantheticalThing, 2, nchar(parantheticalThing) - 1)
   matches <- gregexpr("\\((?>[^()]|(?R))*\\)", text, perl=T)
   if (do.call(sum, matches) <= -1) {return(text)}
   matched <- str_sub(text, matches[[1]], 
                      matches[[1]] + attr(matches[[1]], "match.length") - 1)
   testFunction(matched)
}

testNewScore <- function(datTab, scores, maxFDR = 1, datTab2=NULL) {
   #score is vector of score names
   score.test <- map_dfr(scores, function(x) {
      e.table <- generateErrorTable(datTab, classifier = x)
      e.table$classifier <- x
      e.table <- e.table %>% filter(!is.infinite(fdr), total >= 0.05 * max(total))
      return(e.table)
   })
   if (!is.null(datTab2)) {
      score.test2 <- map_dfr(scores, function(x) {
         e.table <- generateErrorTable(datTab2, x)
         e.table$classifier <- str_c(x, "2")
         e.table <- e.table %>% filter(!is.infinite(fdr), total >= 0.05 * max(total))
         return(e.table)
      })
      score.test <- rbind(score.test, score.test2)
   }
   score.plot <- score.test %>% 
      filter(fdr <= maxFDR) %>%
      ggplot(aes(x=fdr, y=total, col=classifier)) + 
      geom_path()
   print(score.plot)
}

params.best <- c("Score.Diff", "percMatched", "massError", "z", "numPPSM", "numCSM", "xlinkClass")

vh <- read_tsv("DemoFiles/vh7gv76czk/vh.test.1.txt")
for (i in 2:20) {
   dump <- read_tsv(str_c("DemoFiles/vh7gv76czk/vh.test.", i, ".txt"))
   vh <- rbind(vh, dump)
}

vh$Decoy2 <- as.factor(vh$Decoy2)
vh <- buildClassifier(vh)
vh <- buildClassifierExperimental(vh, params.best, "SVM.best")
testNewScore(vh, c("Score.Diff", "dvals", "SVM.best"), 0.2)

vh$pep2.percRun <- vh$pep2.longest.run / vh$Len.Pep.2
vh$pep1.percRun <- vh$pep1.longest.run / vh$Len.Pep.1

#vh$mScore <- vh$pep2.percRun * vh$percInt.pep2
#vh$mScore2 <- log(vh$Rk.1) + log(vh$Rk.2)
#vh <- buildClassifierExperimental(vh, c("Score.Diff", "mScore", "mScore2", "z", "massError"), "SVM.mScore")

vh <- buildClassifierExperimental(vh, c(params.best, "pep2.longest.run"), "SVM.p2lr")
vh <- buildClassifierExperimental(vh, c(params.best, "pep2.percRun"), "SVM.p2pr")
testNewScore(vh, c("Score.Diff", "dvals", "SVM.best", "SVM.p2pr"))
#testNewScore(vh, c("Score.Diff", "dvals", "SVM.best", "SVM.p2pr", "SVM.mScore"))
vh.rp <- bestResPair(vh)
testNewScore(vh.rp, c("Score.Diff", "dvals", "SVM.best", "SVM.p2pr"), 0.2)

vh %>% ggplot(aes(percDiag.pep2 / percMatched, fill=Decoy)) + geom_histogram(col="white", binwidth=0.025) +
   facet_grid(rows="Decoy", scales="free_y")

vh %>% ggplot(aes(percInt.pep2, fill=Decoy)) + geom_histogram(col="white", binwidth=0.025) +
   facet_grid(rows="Decoy", scales="free_y")

vh %>% ggplot(aes(z, fill=Decoy)) + geom_histogram(col="white", binwidth=1) +
   facet_grid(rows="Decoy", scales="free_y")

vh %>% filter(between(dvals, -0.5, 2.5)) %>%
   ggplot(aes(pep2.percRun * percInt.pep2, fill=Decoy)) + geom_histogram(col="white", binwidth=0.01) +
   facet_grid(rows="Decoy", scales="free_y")

vh %>% filter(between(dvals, -0.5, 2.5)) %>%
   ggplot(aes(log(Rk.1 * Rk.2), fill=Decoy)) + geom_histogram(col="white", binwidth=0.5) +
   facet_grid(rows="Decoy", scales="free_y")

vh %>% filter(between(dvals, -0.5, 2.5), Decoy!="DoubleDecoy") %>% 
   ggplot(aes((percDiag.pep1  + percDiag.pep2), group=Decoy)) + 
   geom_histogram(aes(fill=Decoy), binwidth=0.05, position=position_dodge2(preserve="single"), col="white") + 
   facet_grid(rows="Decoy", scales="free_y")

vh %>% filter(between(dvals, -0.5, 2.5), Decoy!="DoubleDecoy") %>% 
   ggplot(aes(percMatched, group=Decoy)) + 
   geom_histogram(aes(fill=Decoy), binwidth=0.1, position=position_dodge2(preserve="single"), col="white") #+
#   facet_grid(rows="z", scales="free_y")

vh %>% filter(between(dvals, -0.5, 2.5), Decoy!="DoubleDecoy") %>% 
   ggplot(aes(x=jitter(pep2.longest.run), y=percMatched, group=Decoy)) + 
   geom_point(aes(col=Decoy)) +
   facet_grid(rows="pep2.longest.run", scales="fixed")



vh8 <- read_tsv("DemoFiles/5155555555/vh_max8.test.1.txt")
for (i in 2:20) {
   dump <- read_tsv(str_c("DemoFiles/5155555555/vh_max8.test.", i, ".txt"))
   vh8 <- rbind(vh8, dump)
}
vh8$Decoy2 <- as.factor(vh8$Decoy2)
vh8 <- buildClassifier(vh8)
vh8 <- buildClassifierExperimental(vh8, params.best, "SVM.best")
vh8.rp <- bestResPair(vh8)
testNewScore(vh8, c("Score.Diff", "dvals", "SVM.best"), 0.2)
testNewScore(vh, c("Score.Diff", "dvals", "SVM.best"), 0.1, vh8)
testNewScore(vh.rp, c("Score.Diff", "dvals", "SVM.best"), 0.2, vh8.rp)


################################################################################
##  Peter asked me to compare where extra CSMs were coming from when he tried
##  restricting the protein-pairs for crosslinking search...
################################################################################

vh51 <- readProspectorXLOutput("DemoFiles/5155555555/tstone2.txt")
vh51 <- buildClassifier(vh51)

vh52 <- readProspectorXLOutput("DemoFiles/5255555555/tstone2.txt")
vh52 <- buildClassifier(vh52)

vh51.thresh <- vh51 %>% 
   lengthFilter(4, 30) %>% 
   filter(Score.Diff > 5) %>%
   findThreshold(targetER = 0.01)
vh51.thresh
vh52.thresh <- vh52 %>% 
   lengthFilter(4, 30) %>% 
   filter(Score.Diff > 5) %>%
   findThreshold(targetER = 0.01)
vh52.thresh
testNewScore(vh51, "dvals", datTab2=vh52)

vh52.1csm <- vh52 %>%
   lengthFilter(4, 30) %>%
   filter(Score.Diff > 5,
          Decoy == "Target",
          dvals >= vh52.thresh[[1]]
          )
vh52.1csm.rp <- bestResPair(vh52.1csm)

vh51.1csm <- vh51 %>%
   lengthFilter(4, 30) %>%
   filter(Score.Diff > 5,
          Decoy == "Target",
          dvals >= vh51.thresh[[1]]
   )
vh51.1csm.rp <- bestResPair(vh51.1csm)

summaryPlot2 <- function(summarizedData) {
   accPlot <- summarizedData %>% 
      ggplot(aes(Acc.1, Acc.2, size=abs(counts), col=counts)) + 
      geom_point(na.rm=T) +
      scale_size(range = c(0, 12)) +
      scale_color_viridis_c(option="D") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, hjust=1))
   print(accPlot)
}

vh51.prot <- summarizeProtData(vh51.1csm)
vh52.prot <- summarizeProtData(vh52.1csm)
vh51.prot$counts[is.na(vh51.prot$counts)] <- 0
vh52.prot$counts[is.na(vh52.prot$counts)] <- 0
vh53.prot <- full_join(vh51.prot, vh52.prot, by=c("Acc.1", "Acc.2"))
vh53.prot$counts <- vh53.prot$counts.y - vh53.prot$counts.x
vh53.prot$counts[vh53.prot$counts == 0] <- NA
summaryPlot2(vh53.prot)
vh53.prot %>% filter(!is.na(counts)) %>% ggplot(aes(counts)) + geom_histogram(binwidth=5)
which(vh53.prot$counts < -20)

###############################################################################
# Residue pair level analysis

vh51.1csm.rp$xlinkedResPair <- as.character(vh51.1csm.rp$xlinkedResPair)
vh51.1csm.rp$source <- "vh51"
vh52.1csm.rp$xlinkedResPair <- as.character(vh52.1csm.rp$xlinkedResPair)
vh52.1csm.rp$source <- "vh52"
vh53.1csm.rp <- rbind(vh51.1csm.rp, vh52.1csm.rp)
vh53.1csm.rp %>% group_by(xlinkedResPair, source) %>% 
   summarize(counts = sum(numCSM),
             pep1.len = mean(Len.Pep.1),
             pep2.len = mean(Len.Pep.2)) %>% 
   ungroup() %>%
   pivot_wider(names_from = source, values_from = counts) %>%
   mutate(vh51 = replace_na(vh51, replace = 0),
          vh52 = replace_na(vh52, replace = 0)) %>%
   mutate(diff = vh52 - vh51) %>% arrange(desc(diff))
   
   
################################################################################

rm1 <- readProspectorXLOutput("DemoFiles/rm11111111/tstone2.txt")
rm1 <- rm1 %>% filter(str_detect(Fraction, "hSAX_f28t30_rep"))
rm1 <- buildClassifier(rm1)
rm1 <- lengthFilter(rm1, 4, 30)
rm1.test <- rm1 %>% filter(Score.Diff > 5)
rm1.test <- bestResPair(rm1.test)

rm2 <- readProspectorXLOutput("DemoFiles/rm22222222/tstone2.txt")
rm2 <- rm2 %>% filter(str_detect(Fraction, "hSAX_f28t30_rep"))# %>% filter(Score.Diff > 5)
rm2 <- buildClassifier(rm2)
rm2 <- lengthFilter(rm2, 4, 30)
rm2.test <- rm2 %>% filter(Score.Diff > 5)
rm2.test <- bestResPair(rm2.test)

rm1.thresh <- findThreshold(rm1.test, targetER = 0.01, minThreshold = -5, classifier = "dvals", errorFUN = calculateFDR, scalingFactor = 10)
rm2.thresh <- findThreshold(rm2.test, targetER = 0.01, minThreshold = -5, classifier = "dvals", errorFUN = calculateFDR, scalingFactor = 1)
calculateHits(removeDecoys(rm1.test), threshold = rm1.thresh[[1]])
calculateHits(removeDecoys(rm2.test), threshold = rm2.thresh[[1]])

rm1.thresh.sep <- findSeparateThresholds(rm1.test, targetER = 0.01, minThreshold = -5, classifier = "dvals", errorFUN = calculateFDR, scalingFactor = 10)
rm2.thresh.sep <- findSeparateThresholds(rm2.test, targetER = 0.01, minThreshold = -5, classifier = "dvals", errorFUN = calculateFDR, scalingFactor = 1)
calculateSeparateHits(removeDecoys(rm1.test), thresholds = rm1.thresh.sep, classifier="dvals") 
calculateSeparateHits(removeDecoys(rm2.test), thresholds = rm2.thresh.sep, classifier="dvals") 

generateErrorTable(rm1, scalingFactor=10) %>% filter(fdr <= 0.5) %>% numHitsPlot
generateErrorTable(rm2, scalingFactor=1) %>% filter(fdr <= 0.5) %>% numHitsPlot

test <- seq(0, 0.25, by=0.01) %>% 
   map(function(x) findSeparateThresholds(rm, targetER = x, classifier = "dvals"))
test2 <- test %>% 
   map_dfr(function(x) calculateSeparateHits(removeDecoys(rm), x, classifier = "dvals")) %>%
   mutate(total = intra + inter)
test3 <- test %>% 
   map_dfr(function(x) calculateSeparateHits(filter(rm, Decoy == "Decoy"), x, classifier = "dvals")) %>%
   transmute(totalTD = (intra + inter) / 10)
test4 <- test %>% 
   map_dfr(function(x) calculateSeparateHits(filter(rm, Decoy == "DoubleDecoy"), x, classifier = "dvals")) %>%
   transmute(totalDD = (intra + inter) / 100)
test5 <- cbind(test2, test3, test4) %>% 
   mutate(fdr = (totalTD + totalDD)/total) %>%
   select(-totalDD, -totalTD)
numHitsPlot(test5, 0.01)

rm.et.inter <- generateErrorTable(filter(rm, str_detect(xlinkClass, "interProtein")))
rm.et.intra <- generateErrorTable(filter(rm, str_detect(xlinkClass, "intraProtein")))
max.fdr <- mmax(max(rm.et.inter$fdr, rm.et.intra$fdr), 0.05)
fdr.spacing <- max.fdr / 50
fdr.seq <- seq(0, max.fdr, fdr.spacing)

sepFDR <- fdr.seq %>% map_dfr(function(i) {
   intra <- rm.et.intra %>% slice(which.min(abs(fdr - i))) %>%
      pull(intra)
  inter <- rm.et.inter %>% slice(which.min(abs(fdr - i))) %>%
     pull(inter)
  tableLine <- tibble(thresh = NA, intra, inter, fdr = i, total = intra + inter)
})

################################################################################
# Mod File Creation for Uniprot Rabbit Ribosome
###############################################################################

old.chain <- read_tsv("DemoFiles/chain.txt", col_names = F)
old.mod <- read_tsv("DemoFiles/modules.txt")
new.mod <- read_tsv("DemoFiles/dssoMS3_unprotRabbit_full.txt")
new.mod$module <- ifelse(str_detect(new.mod$`Protein names`, "40S"), 
                         "40S Ribosome", 
                         ifelse(str_detect(new.mod$`Protein names`, "60S"),
                                "60S Ribosome", "other"))
mods <- left_join(old.mod, new.mod, by=c("Module"="module", "Range_high"="Length"))
mods %>% group_by(Subunit) %>% count() %>% arrange(desc(n))

#fast <- read_file("DemoFiles/dssoMS3_uniprotRabbit.fasta")
fast <- read_file("DemoFiles/UniProtKB.2020.12.01.RABIT.fasta")
fast1 <- str_split(fast, ">")[[1]][2:45455]
fast1.names <- str_extract(fast1, "^(tr|sp)\\|.+\\n")
for (i in seq_along(fast1)) {
   fast1[i] <- str_replace(fast1[i], fixed(fast1.names[i]), "")
}
fast1.names <- str_replace_all(fast1.names, fixed("\n"), "")
fast1 <- str_replace_all(fast1, fixed("\n"), "")
names(fast1) <- fast1.names

rRibo <- read_file("DemoFiles/PA.rRibo.fasta")
rRibo <- str_split(rRibo, ">")[[1]][2:78]
rRibo.names <- str_extract(rRibo, "^[[A-Za-z]].+?\\n")
for (i in seq_along(rRibo)) {
   rRibo[i] <- str_replace(rRibo[i], fixed(rRibo.names[i]), "")
}
rRibo.names <- str_replace_all(rRibo.names, fixed("\n"), "")
rRibo <- str_replace_all(rRibo, fixed("\n"), "")
names(rRibo) <- rRibo.names

table(rRibo %in% fast1)
which(!rRibo %in% fast1)
#test <- rRibo %>% map(function(x) str_which(fast1, x))
#table(map(test, is.na))





###############################################################################
##                         Compare to Clinton MS3                            ##
###############################################################################

clin <- read_tsv("DemoFiles/dssoMS3/clintonsSearch.txt")
clin$missed <- str_count(clin$Sequence, "[[KR]]") - str_count(clin$Sequence, "[[KR]][[P]]") - str_count(clin$Sequence, "[[KR]]$")
me <- read_tsv("DemoFiles/dssoMS3/ms3data.txt", skip=2)
me$missed <- str_count(me$`DB Peptide`, "[[KR]]") - str_count(me$`DB Peptide`, "[[KR]][[P]]") - str_count(me$`DB Peptide`, "[[KR]]$")
me$id <- str_replace(
   str_c(
      str_replace(
         me$Fraction, "_ITMSms3cid", ""
         ), me$`MSMS Info`, sep = ":"
      ), "-", "_"
   )


msf1 <- createMasterScanFile(ms3files2[1], ms2files2[1]) %>% 
   mutate(Fraction = str_extract(ms3files2[1], "Z[[0-9]]{8}\\-[[0-9]]{2}"))
msf2 <- createMasterScanFile(ms3files2[2], ms2files2[2]) %>% 
   mutate(Fraction = str_extract(ms3files2[2], "Z[[0-9]]{8}\\-[[0-9]]{2}"))
msf3 <- createMasterScanFile(ms3files2[3], ms2files2[3]) %>% 
   mutate(Fraction = str_extract(ms3files2[3], "Z[[0-9]]{8}\\-[[0-9]]{2}"))
msf4 <- createMasterScanFile(ms3files2[4], ms2files2[4]) %>% 
   mutate(Fraction = str_extract(ms3files2[4], "Z[[0-9]]{8}\\-[[0-9]]{2}"))
msf <- rbind(msf1,msf2,msf3,msf4)

me2 <- me %>% mutate(Fraction = str_replace(Fraction, "_ITMSms3cid", ""))
names(me2)[which(names(me2)=="MSMS Info")] <- "ms3ScanNo"
me2 <- left_join(me2, msf[,c(2,3,9)], by=c("Fraction", "ms3ScanNo"))
me2$id <- str_c(me2$Fraction, me2$ms2cidScanNo, sep=":") %>% str_replace("-", "_")

length(intersect(me2$id, clin$`Search:Scan`))
length(setdiff(me2$id, clin$`Search:Scan`))
length(setdiff(clin$`Search:Scan`, me2$id))

clinUnique <- filter(clin, `Search:Scan` %in% setdiff(clin$`Search:Scan`, me2$id))
meUnique <- filter(me2, id %in% setdiff(me2$id, clin$`Search:Scan`))
clinUnique %>% 
   filter(missed <= 2, EV <= 0.25, abs(PPM) <= 15) %>% 
   mutate(sulf = str_detect(`Protein Mods`, "Sulfenic")) %>% 
   #group_by(sulf) %>% 
   #count()
   filter(!sulf) %>% View

#write_tsv(clinUnique, "clin_unique.txt")
#write_tsv(meUnique, "me_unique.txt")

clinXL <- read_tsv("clin_xl.txt")
clinXL$`Search:Scan` <- str_replace(clinXL$`Search:Scan`, "_", "-")

names(clinXL) <- names(clinXL) %>% 
   str_replace_all("[[:space:]]", ".") %>% 
   str_replace_all("A$", "1") %>% 
   str_replace("B$", "2")

#test <- clinXL %>% separate(XL.1, c("XL.1.a", "XL.1.b", "XL.1.c"), ";")

demo$id <- str_c(str_replace(demo$Fraction, "_ITMSms3cid", ""), demo$ms2cidScanNo, sep=":")
agree <- intersect(unique(demo$id), unique(clinXL$`Search:Scan`))
length(agree)
length(unique(demo$id))
length(unique(clinXL$`Search:Scan`))
clinXLunique <- clinXL %>% filter(!`Search:Scan` %in% agree)
#   ggplot(aes(`XL PPM`)) + geom_histogram(binwidth=1, col="white")


clinXLunique %>% filter((str_detect(`Protein.Mods.1`, "Alkene@[[0-9\\|]]+") &
                            str_detect(`Protein.Mods.2`, "Thiol\\(Unsaturated\\)@[[0-9\\|]]+") |
                            str_detect(`Protein.Mods.2`, "Sulfenic@[[0-9\\|]]+")) |
                           (str_detect(`Protein.Mods.2`, "Alkene@[[0-9\\|]]+") &
                               str_detect(`Protein.Mods.1`, "Thiol\\(Unsaturated\\)@[[0-9\\|]]+") |
                               str_detect(`Protein.Mods.1`, "Sulfenic@[[0-9\\|]]+"))
) %>% write_tsv("clin_unique.txt")




scanGroupTibble$alkene <- str_count(scanGroupTibble$Peptide.ms3, fixed("XL:A-Alkene"))
scanGroupTibble$thiol <- str_count(scanGroupTibble$Peptide.ms3, fixed("XL:A-Thiol(Unsaturated)"))
scanGroupTibble$sulfenic <- str_count(scanGroupTibble$Peptide.ms3, fixed("XL:A-Sulfenic"))






###########################################
### monocle testing
##########################################

pk1 <- tibble(fdr = rep(c(1,5),4), 
              incorrectMono = rep(c(rep(T, 2), rep(F, 2)), 2),
              algorithm = c(rep("PAVA", 4), rep("MONOCLE", 4)),
              numCSM = c(4641, 5365, 3546, 4030, 4947, 5814, 4502, 5221),
              intra = c(3361, 3361, 2690, 2690, 3595, 3595, 3342, 3342),
              inter = c(1280, 2004, 856, 1340, 1352, 2219, 1160, 1879),
              vr = c(3.86, 6.71, 3.32, 5.75, 3.68, 7.21, 3.57, 6.68)
)

#pdf("pk1.pdf", useDingbats = F)
pk1 %>% 
   ggplot(aes(x=incorrectMono, y=numCSM, fill=algorithm)) + 
   geom_col(position=position_dodge2()) + 
   geom_text(aes(label=numCSM), position = position_dodge(width=0.9), vjust=-0.25) + 
   facet_grid("fdr") + 
   scale_fill_discrete(type=c("lightseagreen", "aquamarine")) + 
   theme_bw()
#dev.off()

pk1.xl <- tibble(fdr = rep(c(1,5),4), 
              incorrectMono = rep(c(rep(T, 2), rep(F, 2)), 2),
              algorithm = c(rep("PAVA", 4), rep("MONOCLE", 4)),
              numXL = c(686, 879, 621, 780, 705, 893, 712, 870),
              intra = c(504, 562, 470, 536, 532, 564, 529, 556),
              inter = c(182, 317, 151, 244, 173, 329, 183, 314),
              vr = c(6.99, 11.39, 6.37, 10.3, 6.05, 11.74, 7.53, 11.66)
)

#pdf("pk1.xl.pdf", useDingbats = F)
pk1.xl %>% 
   ggplot(aes(x=incorrectMono, y=numXL, fill=algorithm)) + 
   geom_col(position=position_dodge2()) + 
   geom_text(aes(label=numXL), position = position_dodge(width=0.9), vjust=-0.25) + 
   facet_grid("fdr") + 
   scale_fill_discrete(type=c("lightseagreen", "aquamarine")) + 
   theme_bw()
#dev.off()

test <- read_tsv("~/Projects-Collaboration/BurlingameLab/Jason/monocle_test/tstone2-2.txt")
test$incor <- str_detect(test$Peptide.1, fixed("Incorrect Mono"))
table(test$incor)
test %>% ggplot(aes(ppm, fill=incor)) + geom_histogram(color="white", binwidth = 1)
test$Da.error <- test$`m/z` * test$ppm * 1e-6 * test$z
#pdf("massErrorDa.pdf", useDingbats = F)
test %>% ggplot(aes(Da.error, fill=incor)) + geom_histogram(color="white", binwidth = 0.002)
#dev.off()
test %>% group_by(incor) %>% summarize(mean=mean(Da.error, na.rm=T))

test <- read_tsv("~/Projects-Collaboration/BurlingameLab/Jason/monocle_test/linearIncorrect2.txt",skip=2)
test$incor <- str_detect(test$`Variable Mods`, fixed("Incorrect Mono"))
test$mmod <- str_detect(test$`Variable Mods`, "[[0-9]]{3,4}\\.[[0-9]]+")
table(test$mmod, test$incor)
test %>% filter(!mmod) %>%
   ggplot(aes(ppm, fill=incor)) + geom_histogram(color="white", binwidth = 1)
test$Da.error <- test$`m/z` * test$ppm * 1e-6 * test$z
#pdf("massErrorDaLinear.pdf", useDingbats = F)
test %>% filter(!mmod) %>%
   ggplot(aes(Da.error, fill=incor)) + 
   geom_histogram(color="white", binwidth = 0.002, position=position_dodge())
#dev.off()
test %>% filter(!mmod) %>%
   group_by(incor) %>% summarize(mean=mean(Da.error, na.rm=T))

test %>% filter(!mmod) %>% View

