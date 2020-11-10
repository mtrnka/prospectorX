source("touchStone.R")
scResults2 <- "DemoFiles/dssoMS3/ms3data.txt"
ms3files2 <- c("DemoFiles/dssoMS3/peaklists/Z20200519-31_ITMSms3cid.txt",
              "DemoFiles/dssoMS3/peaklists/Z20200519-39_ITMSms3cid.txt",
              "DemoFiles/dssoMS3/peaklists/Z20200519-49_ITMSms3cid.txt",
              "DemoFiles/dssoMS3/peaklists/Z20200519-59_ITMSms3cid.txt")

ms2files2 <- c("DemoFiles/dssoMS3/peaklists/Z20200519-31_FTMSms2cid.txt",
              "DemoFiles/dssoMS3/peaklists/Z20200519-39_FTMSms2cid.txt",
              "DemoFiles/dssoMS3/peaklists/Z20200519-49_FTMSms2cid.txt",
              "DemoFiles/dssoMS3/peaklists/Z20200519-59_FTMSms2cid.txt")
demo <- processMS3xlinkResultsMultiFile(scResults2, ms3files2, ms2files2)

searchCompareMS3 <- readMS3results(scResults2)
searchCompareMS3 <- searchCompareMS3 %>% group_by(Fraction.ms3) %>% nest()
baseFiles <- str_replace(searchCompareMS3$Fraction.ms3, "_[[A-Z]]+MSms[[0-9]][[a-z]]+", "")
p1 <- searchCompareMS3 %>% filter(str_detect(Fraction.ms3, baseFiles[2])) %>% unnest(cols=c(data))
p2 <- ms3files2[str_which(ms3files2, baseFiles[2])]
p3 <- ms2files2[str_which(ms2files2, baseFiles[2])]

# processMS3xlinkResultsSingleFile(p1, p2, p3)


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
      geom_line()
   print(score.plot)
}

params.best <- c("Score.Diff", "percMatched", "massError", "z", "numPPSM", "numCSM", "xlinkClass")
